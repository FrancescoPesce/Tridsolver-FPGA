#include "thomas_accel.h"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <iostream>

namespace {

// TDMA_batch arguments (ThomasVsPcr/Thomas/TDMA_solver.cpp)
enum TDMA_Arg { ARG_A = 0, ARG_B, ARG_C, ARG_D, ARG_U, ARG_M, ARG_N, ARG_BATCH, ARG_ITERS };

const int kMaxLength = 128;            // N_MAX of the kernel
// Systems per call: thomas_interleave computes (bat << 5) on an ap_uint<12>, which wraps after
// 127 groups of 256 systems and makes the kernel read past its input (a hang on hardware).
const size_t kMaxRows = 127 * 256;
const size_t kBankBytes = 256ul << 20; // TDMA.ini connects every argument to one 256 MB HBM bank

size_t round_up( size_t value, size_t multiple ) { return ( value + multiple - 1 ) / multiple * multiple; }

double seconds_since( std::chrono::steady_clock::time_point start )
{ return std::chrono::duration<double>( std::chrono::steady_clock::now() - start ).count(); }

bool same( const Thomas_Accelerator::Coefficients& x, const Thomas_Accelerator::Coefficients& y )
{
	return x.off_diagonal == y.off_diagonal && x.diagonal == y.diagonal
	       && x.end_diagonal == y.end_diagonal && x.point_diagonal == y.point_diagonal;
}

} // namespace

void Thomas_Accelerator::attach( cl::Context* context, cl::CommandQueue* queue, const std::vector<cl::Kernel>& compute_units )
{
	release_buffers();
	context_ = context;
	queue_ = queue;
	cus_ = compute_units;
}

void Thomas_Accelerator::release( void )
{
	if( solves_ > 0 )
	{
		std::cout << "Thomas FPGA summary: " << solves_ << " sweeps on TDMA_batch (" << seconds_
		          << " s incl. transfers, " << kernel_seconds_ << " s kernel)" << std::endl;
	}
	release_buffers();
	cus_.clear();
	context_ = nullptr;
	queue_ = nullptr;
}

void Thomas_Accelerator::release_buffers( void )
{
	if( queue_ != nullptr )
	{
		cl_int err;
		bool unmapped = false;
		for( Sweep& sweep : sweeps_ )
		{
			for( Chunk& chunk : sweep.chunks )
			{
				if( chunk.d_ptr != nullptr )
				{ OCL_CHECK( err, err = queue_->enqueueUnmapMemObject( chunk.d, chunk.d_ptr ) ); unmapped = true; }
				if( chunk.u_ptr != nullptr )
				{ OCL_CHECK( err, err = queue_->enqueueUnmapMemObject( chunk.u, chunk.u_ptr ) ); unmapped = true; }
			}
		}
		if( unmapped )
		{ OCL_CHECK( err, err = queue_->finish() ); }
	}
	for( Sweep& sweep : sweeps_ )
	{ sweep = Sweep(); }
	configured_ = false;
}

cl::Buffer Thomas_Accelerator::create_buffer( int cu, int arg, size_t bytes )
{
	// XRT places the buffer in the memory bank connected to `arg` of this compute unit
	// (HBM[0..4] for TDMA_batch_1 with TDMA.ini).
	cl_mem_ext_ptr_t ext;
	ext.argidx = arg;
	ext.host_ptr_ = nullptr;
	ext.kernel = cus_[cu]();

	cl_int err;
	OCL_CHECK( err, cl::Buffer buffer( *context_, CL_MEM_READ_WRITE | CL_MEM_EXT_PTR_XILINX, bytes, &ext, &err ) );
	return buffer;
}

bool Thomas_Accelerator::configure( int nx, int ny, int nz, const Coefficients& coefficients, std::string& reason )
{
	if( configured_ && nx == nx_ && ny == ny_ && nz == nz_ && same( coefficients, coefficients_ ) )
	{ return true; }

	release_buffers();
	const int num_substrates = static_cast<int>( coefficients.diagonal.size() );
	if( cus_.empty() || context_ == nullptr || queue_ == nullptr )
	{
		reason = "no TDMA_batch compute unit attached";
		return false;
	}
	if( num_substrates < 1 )
	{
		reason = "no substrates";
		return false;
	}
	if( std::max( nx, std::max( ny, nz ) ) > kMaxLength )
	{
		reason = "grid " + std::to_string( nx ) + "x" + std::to_string( ny ) + "x" + std::to_string( nz )
		       + " exceeds the TDMA_batch limit of " + std::to_string( kMaxLength ) + " voxels per line";
		return false;
	}

	const int n[3] = { nx, ny, nz };
	const size_t lines[3] = { (size_t) ny * nz, (size_t) nx * nz, (size_t) nx * ny };
	const size_t num_cus = cus_.size();

	// Layout: rows of every direction split over the compute units; all limits checked first.
	std::vector<size_t> bytes_per_cu( num_cus, 0 ); // per argument, the three sweeps share the banks
	for( int dir = 0; dir < 3; dir++ )
	{
		Sweep& sweep = sweeps_[dir];
		sweep.n = n[dir];
		sweep.M = static_cast<int>( round_up( n[dir], 16 ) );
		sweep.lines = lines[dir];
		const size_t rows = sweep.lines * num_substrates;
		// calls of at most kMaxRows rows, the same number on every compute unit
		const size_t calls_per_cu = ( rows + num_cus * kMaxRows - 1 ) / ( num_cus * kMaxRows );
		const size_t num_chunks = num_cus * calls_per_cu;
		sweep.rows_per_chunk = ( rows + num_chunks - 1 ) / num_chunks;

		for( size_t c = 0; c * sweep.rows_per_chunk < rows; c++ )
		{
			Chunk chunk;
			chunk.cu = static_cast<int>( c % num_cus );
			chunk.first_row = c * sweep.rows_per_chunk;
			chunk.rows = std::min( sweep.rows_per_chunk, rows - chunk.first_row );
			chunk.N = 8;
			chunk.B = static_cast<int>( round_up( chunk.rows, 8 ) / 8 );
			bytes_per_cu[chunk.cu] += (size_t) chunk.N * chunk.B * sweep.M * sizeof( float );
			sweep.chunks.push_back( chunk );
		}
	}
	for( size_t cu = 0; cu < num_cus; cu++ )
	{
		if( bytes_per_cu[cu] > kBankBytes )
		{
			reason = "the sweeps need " + std::to_string( bytes_per_cu[cu] >> 20 ) + " MB per argument, more than a "
			       + std::to_string( kBankBytes >> 20 ) + " MB HBM bank";
			release_buffers();
			return false;
		}
	}

	cl_int err;
	std::vector<cl::Memory> to_device;
	for( int dir = 0; dir < 3; dir++ )
	{
		Sweep& sweep = sweeps_[dir];
		for( Chunk& chunk : sweep.chunks )
		{
			const size_t rows = (size_t) chunk.N * chunk.B;
			const size_t bytes = rows * sweep.M * sizeof( float );
			chunk.a = create_buffer( chunk.cu, ARG_A, bytes );
			chunk.b = create_buffer( chunk.cu, ARG_B, bytes );
			chunk.c = create_buffer( chunk.cu, ARG_C, bytes );
			chunk.d = create_buffer( chunk.cu, ARG_D, bytes );
			chunk.u = create_buffer( chunk.cu, ARG_U, bytes );

			float* a;
			float* b;
			float* c;
			OCL_CHECK( err, a = (float*) queue_->enqueueMapBuffer( chunk.a, CL_TRUE, CL_MAP_WRITE_INVALIDATE_REGION, 0, bytes, nullptr, nullptr, &err ) );
			OCL_CHECK( err, b = (float*) queue_->enqueueMapBuffer( chunk.b, CL_TRUE, CL_MAP_WRITE_INVALIDATE_REGION, 0, bytes, nullptr, nullptr, &err ) );
			OCL_CHECK( err, c = (float*) queue_->enqueueMapBuffer( chunk.c, CL_TRUE, CL_MAP_WRITE_INVALIDATE_REGION, 0, bytes, nullptr, nullptr, &err ) );
			for( size_t r = 0; r < rows; r++ )
			{
				const size_t row = r * sweep.M;
				const bool padding_row = r >= chunk.rows;
				const int s = padding_row ? 0 : static_cast<int>( ( chunk.first_row + r ) / sweep.lines );
				for( int p = 0; p < sweep.M; p++ )
				{
					float lower = 0.0f, diag = 1.0f, upper = 0.0f; // decoupled identity row
					if( !padding_row && p < sweep.n )
					{
						if( sweep.n == 1 )
						{ diag = coefficients.point_diagonal[s]; }
						else
						{
							const bool first = ( p == 0 ), last = ( p == sweep.n - 1 );
							lower = first ? 0.0f : coefficients.off_diagonal[s];
							upper = last ? 0.0f : coefficients.off_diagonal[s];
							diag = ( first || last ) ? coefficients.end_diagonal[s] : coefficients.diagonal[s];
						}
					}
					a[row + p] = lower;
					b[row + p] = diag;
					c[row + p] = upper;
				}
			}
			OCL_CHECK( err, err = queue_->enqueueUnmapMemObject( chunk.a, a ) );
			OCL_CHECK( err, err = queue_->enqueueUnmapMemObject( chunk.b, b ) );
			OCL_CHECK( err, err = queue_->enqueueUnmapMemObject( chunk.c, c ) );

			// d and u stay mapped: padding entries of d must remain zero, only grid points are written.
			OCL_CHECK( err, chunk.d_ptr = (float*) queue_->enqueueMapBuffer( chunk.d, CL_TRUE, CL_MAP_WRITE, 0, bytes, nullptr, nullptr, &err ) );
			std::memset( chunk.d_ptr, 0, bytes );
			OCL_CHECK( err, chunk.u_ptr = (float*) queue_->enqueueMapBuffer( chunk.u, CL_TRUE, CL_MAP_READ, 0, bytes, nullptr, nullptr, &err ) );

			to_device.push_back( chunk.a );
			to_device.push_back( chunk.b );
			to_device.push_back( chunk.c );
		}
	}
	OCL_CHECK( err, err = queue_->enqueueMigrateMemObjects( to_device, 0 ) );
	OCL_CHECK( err, err = queue_->finish() );

	nx_ = nx;
	ny_ = ny;
	nz_ = nz;
	coefficients_ = coefficients;
	configured_ = true;
	std::cout << "Thomas FPGA: TDMA_batch with " << num_cus << " compute unit(s), grid " << nx << "x" << ny << "x" << nz
	          << ", " << num_substrates << " substrates; x/y/z sweeps: rows of " << sweeps_[X].M << "/" << sweeps_[Y].M
	          << "/" << sweeps_[Z].M << " elements, " << sweeps_[X].chunks.size() << "/" << sweeps_[Y].chunks.size()
	          << "/" << sweeps_[Z].chunks.size() << " kernel calls" << std::endl;
	return true;
}

Thomas_Accelerator::Chunk& Thomas_Accelerator::chunk_of( int direction, int substrate, size_t line, size_t& offset )
{
	Sweep& sweep = sweeps_[direction];
	const size_t row = substrate * sweep.lines + line;
	Chunk& chunk = sweep.chunks[row / sweep.rows_per_chunk];
	offset = ( row - chunk.first_row ) * sweep.M;
	return chunk;
}

float* Thomas_Accelerator::rhs( int direction, int substrate, size_t line )
{
	size_t offset;
	Chunk& chunk = chunk_of( direction, substrate, line, offset );
	return chunk.d_ptr + offset;
}

const float* Thomas_Accelerator::solution( int direction, int substrate, size_t line )
{
	size_t offset;
	Chunk& chunk = chunk_of( direction, substrate, line, offset );
	return chunk.u_ptr + offset;
}

void Thomas_Accelerator::solve( int direction )
{
	const auto start = std::chrono::steady_clock::now();
	Sweep& sweep = sweeps_[direction];
	cl_int err;

	std::vector<cl::Memory> inputs, outputs;
	for( Chunk& chunk : sweep.chunks )
	{
		inputs.push_back( chunk.d );
		outputs.push_back( chunk.u );
	}

	// Arguments are captured at enqueue time: calls sharing a compute unit run one after the other.
	std::vector<cl::Event> h2d( 1 ), kernels( sweep.chunks.size() );
	cl::Event d2h;
	OCL_CHECK( err, err = queue_->enqueueMigrateMemObjects( inputs, 0, nullptr, &h2d[0] ) );
	for( size_t c = 0; c < sweep.chunks.size(); c++ )
	{
		Chunk& chunk = sweep.chunks[c];
		cl::Kernel& k = cus_[chunk.cu];
		OCL_CHECK( err, err = k.setArg( ARG_A, chunk.a ) );
		OCL_CHECK( err, err = k.setArg( ARG_B, chunk.b ) );
		OCL_CHECK( err, err = k.setArg( ARG_C, chunk.c ) );
		OCL_CHECK( err, err = k.setArg( ARG_D, chunk.d ) );
		OCL_CHECK( err, err = k.setArg( ARG_U, chunk.u ) );
		OCL_CHECK( err, err = k.setArg( ARG_M, sweep.M ) );
		OCL_CHECK( err, err = k.setArg( ARG_N, chunk.N ) );
		OCL_CHECK( err, err = k.setArg( ARG_BATCH, chunk.B ) );
		OCL_CHECK( err, err = k.setArg( ARG_ITERS, 1 ) );
		OCL_CHECK( err, err = queue_->enqueueTask( k, &h2d, &kernels[c] ) );
	}
	OCL_CHECK( err, err = queue_->enqueueMigrateMemObjects( outputs, CL_MIGRATE_MEM_OBJECT_HOST, &kernels, &d2h ) );
	OCL_CHECK( err, err = d2h.wait() );

	cl_ulong first_start = ~0ul, last_end = 0;
	bool have_profile = true;
	for( cl::Event& e : kernels )
	{
		cl_ulong t_start = 0, t_end = 0;
		have_profile = have_profile && e.getProfilingInfo( CL_PROFILING_COMMAND_START, &t_start ) == CL_SUCCESS
		                            && e.getProfilingInfo( CL_PROFILING_COMMAND_END, &t_end ) == CL_SUCCESS;
		first_start = std::min( first_start, t_start );
		last_end = std::max( last_end, t_end );
	}
	if( have_profile && last_end > first_start )
	{ kernel_seconds_ += ( last_end - first_start ) * 1e-9; }

	seconds_ += seconds_since( start );
	solves_++;
}
