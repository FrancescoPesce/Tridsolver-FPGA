#ifndef __THOMAS_ACCEL_H__
#define __THOMAS_ACCEL_H__

/*
    Host driver for the TDMA_batch kernel of FPGA/XIlinx/ThomasVsPcr/Thomas (TDMA_solver.cpp),
    a batched Thomas solver with per-element coefficients:

        TDMA_batch(a, b, c, d, u, M, N, B, iters)

    solves the N*B float32 systems  a[i] u[i-1] + b[i] u[i] + c[i] u[i+1] = d[i]  (i = 0..M-1,
    a[0] and c[M-1] ignored) stored as consecutive rows of M elements. It only solves along the
    contiguous dimension and `iters` just repeats the same solve. Constraints (otherwise wrong
    results or a hang): M % 16 == 0, M <= 128, N*B % 8 == 0 and N*B <= 32512 (127 groups of 256
    systems; beyond that a 12-bit counter in thomas_interleave wraps).

    Here it solves the three sweeps of BioFVM's 3-D LOD step. For every direction each
    (substrate, grid line) pair is one row, padded to M with decoupled identity rows
    (a = c = 0, b = 1, d = 0). The coefficients only depend on the grid and the substrates and are
    uploaded once; each sweep uploads d and downloads u. The rows are split into kernel calls of at
    most 32512 systems, spread over the compute units, and the buffers of all three sweeps stay on
    the device (TDMA.ini: one 256 MB HBM bank per argument).
*/

#include "vitis_common.h"

#include <cstddef>
#include <string>
#include <vector>

class Thomas_Accelerator
{
 public:
	enum { X = 0, Y = 1, Z = 2 };

	// Per-substrate coefficients of BioFVM's LOD systems (identical for the three directions)
	struct Coefficients
	{
		std::vector<double> off_diagonal;   // a = c = -D dt / dx^2              (thomas_constant1a)
		std::vector<double> diagonal;       // 1 + 2 D dt / dx^2 + lambda dt / 3 (thomas_constant3)
		std::vector<double> end_diagonal;   // 1 + D dt / dx^2 + lambda dt / 3   (thomas_constant3a)
		std::vector<double> point_diagonal; // 1 + lambda dt / 3, line of a single voxel
	};

	void attach( cl::Context* context, cl::CommandQueue* queue, const std::vector<cl::Kernel>& compute_units );
	// Releases all OpenCL objects and prints a summary. Call before the OpenCL context is destroyed.
	void release( void );

	// Prepares the sweeps of an nx x ny x nz grid (no-op when nothing changed). Returns false,
	// with a reason, when the kernel cannot solve them.
	bool configure( int nx, int ny, int nz, const Coefficients& coefficients, std::string& reason );

	// Row of `substrate` along `direction` through grid line `line` (X: j + ny*k, Y: i + nx*k,
	// Z: i + nx*j). rhs() is filled before solve(), solution() is read after it.
	float* rhs( int direction, int substrate, size_t line );
	void solve( int direction );
	const float* solution( int direction, int substrate, size_t line );

 private:
	// Rows [first_row, first_row + rows) of one direction, solved by one kernel call on compute unit `cu`
	struct Chunk
	{
		int cu = 0;
		size_t first_row = 0, rows = 0;
		int N = 0, B = 0; // kernel arguments, N*B >= rows
		cl::Buffer a, b, c, d, u;
		float* d_ptr = nullptr;
		float* u_ptr = nullptr;
	};
	struct Sweep
	{
		int n = 0;        // grid points along the direction
		int M = 0;        // row length
		size_t lines = 0; // grid lines per substrate
		size_t rows_per_chunk = 0;
		std::vector<Chunk> chunks;
	};

	void release_buffers( void );
	cl::Buffer create_buffer( int cu, int arg, size_t bytes );
	Chunk& chunk_of( int direction, int substrate, size_t line, size_t& offset );

	cl::Context* context_ = nullptr;
	cl::CommandQueue* queue_ = nullptr;
	std::vector<cl::Kernel> cus_;

	bool configured_ = false;
	int nx_ = 0, ny_ = 0, nz_ = 0;
	Coefficients coefficients_;
	Sweep sweeps_[3];

	unsigned long solves_ = 0;
	double seconds_ = 0.0, kernel_seconds_ = 0.0;
};

#endif
