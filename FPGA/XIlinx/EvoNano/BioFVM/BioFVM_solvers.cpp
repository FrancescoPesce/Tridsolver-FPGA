/*
#############################################################################
# If you use BioFVM in your project, please cite BioFVM and the version     #
# number, such as below:                                                    #
#                                                                           #
# We solved the diffusion equations using BioFVM (Version 1.1.7) [1]        #
#                                                                           #
# [1] A. Ghaffarizadeh, S.H. Friedman, and P. Macklin, BioFVM: an efficient #
#    parallelized diffusive transport solver for 3-D biological simulations,#
#    Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730 #
#                                                                           #
#############################################################################
#                                                                           #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)   #
#                                                                           #
# Copyright (c) 2015-2017, Paul Macklin and the BioFVM Project              #
# All rights reserved.                                                      #
#                                                                           #
# Redistribution and use in source and binary forms, with or without        #
# modification, are permitted provided that the following conditions are    #
# met:                                                                      #
#                                                                           #
# 1. Redistributions of source code must retain the above copyright notice, #
# this list of conditions and the following disclaimer.                     #
#                                                                           #
# 2. Redistributions in binary form must reproduce the above copyright      #
# notice, this list of conditions and the following disclaimer in the       #
# documentation and/or other materials provided with the distribution.      #
#                                                                           #
# 3. Neither the name of the copyright holder nor the names of its          #
# contributors may be used to endorse or promote products derived from this #
# software without specific prior written permission.                       #
#                                                                           #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       #
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED #
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           #
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER #
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  #
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       #
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        #
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    #
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      #
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        #
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              #
#                                                                           #
#############################################################################
*/

#include "BioFVM_solvers.h" 
#include "BioFVM_agent_container.h"
#include "BioFVM_microenvironment.h"
#include "BioFVM_vector.h" 

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <string>

// Per-step timing, printed like the krnl_compute version (Ours/BioFVM/BioFVM_solvers.cpp)
#define DEBUG_PRINT 1

namespace BioFVM{

// do I even need this? 
void diffusion_decay_solver__constant_coefficients_explicit( Microenvironment& M, double dt )
{
	static bool precomputations_and_constants_done = false; 
	if( !precomputations_and_constants_done )
	{
		std::cout	<< std::endl << "Using solver: " << __FUNCTION__ << std::endl 
					<< "     (constant diffusion coefficient with explicit stepping, implicit decay) ... " << std::endl << std::endl;  

		if( M.mesh.uniform_mesh == true )
		{
			std::cout << "Uniform mesh detected! Consider switching to a more efficient method, such as " << std::endl  
			<< "     diffusion_decay_solver__constant_coefficients_explicit_uniform_mesh" << std::endl  
			<< std::endl; 
		}

		precomputations_and_constants_done = true; 
	}

	return; 
}

void diffusion_decay_solver__constant_coefficients_explicit_uniform_mesh( Microenvironment& M, double dt )
{
	static bool precomputations_and_constants_done = false; 
	if( !precomputations_and_constants_done )
	{
		std::cout	<< std::endl << "Using solver: " << __FUNCTION__ << std::endl 
					<< "     (constant diffusion coefficient with explicit stepping, implicit decay, uniform mesh) ... " << std::endl << std::endl;  

		if( M.mesh.regular_mesh == false )
		{ std::cout << "Error. This code is only supported for regular meshes." << std::endl; }

		precomputations_and_constants_done = true; 
	}

	return; 
}

int compute_voxel_idx(int i, int j, int k, int y_coord_size, int x_coord_size) {
    return ( k * y_coord_size + j ) * x_coord_size + i;
}

void print_timestamp() {
	using namespace std::chrono;

    auto now = system_clock::now();
    auto now_time = system_clock::to_time_t(now);

    auto ms = duration_cast<microseconds>(now.time_since_epoch()) % seconds(1);

    std::tm tm = *std::localtime(&now_time);

    std::cout << std::put_time(&tm, "%Y-%m-%dT%H:%M:%S")
              << "." << std::setw(6) << std::setfill('0') << ms.count()
              << '\n';
}

namespace {

typedef std::chrono::high_resolution_clock::duration Duration;

long long to_us( Duration d )
{ return std::chrono::duration_cast<std::chrono::microseconds>( d ).count(); }

// Runs `work`; with DEBUG_PRINT its duration is added to `total`
template <class Work>
void timed( Duration& total, Work work )
{
#if DEBUG_PRINT
	const auto start = std::chrono::high_resolution_clock::now();
	work();
	total += std::chrono::high_resolution_clock::now() - start;
#else
	work();
#endif
}

// Grid lines along a direction: `points` voxels per line, `stride` between them and the first
// voxel of each line (BioFVM voxel index (k*ny + j)*nx + i)
struct Sweep_Lines
{
	int direction;
	size_t nx, ny, points, lines, stride;

	Sweep_Lines( Cartesian_Mesh& mesh, int dir ) : direction( dir )
	{
		nx = mesh.x_coordinates.size();
		ny = mesh.y_coordinates.size();
		const size_t nz = mesh.z_coordinates.size();
		switch( direction )
		{
			case Thomas_Accelerator::X: points = nx; lines = ny * nz; stride = 1; break;
			case Thomas_Accelerator::Y: points = ny; lines = nx * nz; stride = nx; break;
			default: points = nz; lines = nx * ny; stride = nx * ny; break;
		}
	}

	size_t first_voxel( size_t line ) const
	{
		switch( direction )
		{
			case Thomas_Accelerator::X: return line * nx; // line = j + ny*k
			case Thomas_Accelerator::Y: return ( line / nx ) * nx * ny + line % nx; // line = i + nx*k
			default: return line; // line = i + nx*j
		}
	}
};

} // namespace

/*
    3-D LOD step (BioFVM's algorithm) with the tridiagonal sweeps on the TDMA_batch kernel of
    FPGA/XIlinx/ThomasVsPcr/Thomas. The kernel solves float32 systems along rows only: every sweep
    gathers the grid lines of all substrates into rows and scatters the solution back, and the
    Dirichlet nodes are applied on the host between the sweeps, as in BioFVM. Grids the kernel
    cannot handle (more than 128 voxels along a direction) use BioFVM's CPU sweeps.

    With DEBUG_PRINT the step is timed like the krnl_compute version (Ours/BioFVM/BioFVM_solvers.cpp):
    std::chrono around each phase, with the device work of a phase finished before its timer stops,
    printed in microseconds, plus wall-clock timestamps at the start and end of the step and of every
    kernel execution (one per sweep). Phases repeated by the three sweeps are summed over the step:
      Setup time                    coefficient vectors
      Thomas setup time             configure(): buffers and coefficient rows when the grid changes
      Host Dirichlet time           the four apply_dirichlet_conditions() (done by the kernel in Ours)
      Density setup time            gather of the grid lines into rows
      Host-to-Device transfer time  right-hand sides (and new coefficients)
      Kernel time                   kernel calls
      Device-to-Host transfer time  solutions
      Density teardown time         scatter of the solutions into the grid
      Compute time                  whole step
    EVONANO_THOMAS_VERIFY=1 also runs the CPU step on a copy, outside the timed region, and reports
    the difference.
*/
void Microenvironment::compute( void )
{
	static const char* verify_env = std::getenv( "EVONANO_THOMAS_VERIFY" ); 
	static const bool verify = verify_env != NULL && std::string( verify_env ) != "" && std::string( verify_env ) != "0"; 
	std::vector< std::vector<double> > reference; 
	if( verify )
	{ reference = *p_density_vectors; }

#if DEBUG_PRINT
	auto start = std::chrono::high_resolution_clock::now();
	std::cout << "Start compute: ";
	print_timestamp();
#endif

	const int num_substrates = number_of_densities(); 

	Thomas_Accelerator::Coefficients coefficients; 
	coefficients.off_diagonal = thomas_constant1a; 
	coefficients.diagonal = thomas_constant3; 
	coefficients.end_diagonal = thomas_constant3a; 
	coefficients.point_diagonal = thomas_constant2; 
	for( int s = 0; s < num_substrates; s++ )
	{ coefficients.point_diagonal[s] += 1.0; }

#if DEBUG_PRINT
	auto end_setup = std::chrono::high_resolution_clock::now();
	auto duration_setup = std::chrono::duration_cast<std::chrono::microseconds>(end_setup - start);
	std::cout << "Setup time: " << duration_setup.count() << " us" << std::endl;

	auto start_thomas = std::chrono::high_resolution_clock::now();
#endif

	static std::string cpu_reason; 
	std::string reason; 
	const bool on_fpga = thomas_fpga.configure( mesh.x_coordinates.size(), mesh.y_coordinates.size(), mesh.z_coordinates.size(), coefficients, reason ); 

#if DEBUG_PRINT
	auto end_thomas = std::chrono::high_resolution_clock::now();
	auto duration_thomas = std::chrono::duration_cast<std::chrono::microseconds>(end_thomas - start_thomas);
	std::cout << "Thomas setup time: " << duration_thomas.count() << " us" << std::endl;
#endif

	if( !on_fpga )
	{
		if( reason != cpu_reason )
		{
			std::cout << "Thomas FPGA: " << reason << "; using BioFVM's CPU sweeps" << std::endl; 
			cpu_reason = reason; 
		}
		lod_3d_cpu(); 
	}
	else
	{
		cpu_reason.clear(); 
		Duration dirichlet( 0 ), density( 0 ), h2d( 0 ), kernel( 0 ), d2h( 0 ), teardown( 0 ); 

		for( int direction = Thomas_Accelerator::X; direction <= Thomas_Accelerator::Z; direction++ )
		{
			timed( dirichlet, [&]{ apply_dirichlet_conditions(); } ); 
			timed( density, [&]{ thomas_gather( direction ); } ); 
			timed( h2d, [&]{ thomas_fpga.upload( direction ); } ); 

#if DEBUG_PRINT
			auto start_kernel = std::chrono::high_resolution_clock::now();
			std::cout << "Start kernel: ";
			print_timestamp();
#endif
			thomas_fpga.run( direction ); 
#if DEBUG_PRINT
			kernel += std::chrono::high_resolution_clock::now() - start_kernel; 
			std::cout << "End kernel: ";
			print_timestamp();
#endif

			timed( d2h, [&]{ thomas_fpga.download( direction ); } ); 
			timed( teardown, [&]{ thomas_scatter( direction ); } ); 
		}
		timed( dirichlet, [&]{ apply_dirichlet_conditions(); } ); 

#if DEBUG_PRINT
		std::cout << "Host Dirichlet time: " << to_us( dirichlet ) << " us" << std::endl;
		std::cout << "Density setup time: " << to_us( density ) << " us" << std::endl;
		std::cout << "Host-to-Device transfer time: " << to_us( h2d ) << " us" << std::endl;
		std::cout << "Kernel time: " << to_us( kernel ) << " us" << std::endl;
		std::cout << "Device-to-Host transfer time: " << to_us( d2h ) << " us" << std::endl;
		std::cout << "Density teardown time: " << to_us( teardown ) << " us" << std::endl;
#endif
	}

#if DEBUG_PRINT
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	std::cout << "Compute time: " << duration.count() << " us" << std::endl;
	std::cout << "End compute: ";
	print_timestamp();
#endif

	if( verify && on_fpga )
	{
		std::vector< std::vector<double> >* current = p_density_vectors; 
		p_density_vectors = &reference; 
		lod_3d_cpu(); 
		p_density_vectors = current; 

		// largest difference relative to the largest magnitude of each substrate
		std::vector<double> max_diff( num_substrates, 0.0 ), max_value( num_substrates, 0.0 ); 
		for( unsigned int n = 0; n < reference.size(); n++ )
		{
			for( int s = 0; s < num_substrates; s++ )
			{
				const double diff = std::fabs( (*p_density_vectors)[n][s] - reference[n][s] ); 
				if( !( diff <= max_diff[s] ) )
				{ max_diff[s] = diff; }
				max_value[s] = std::max( max_value[s], std::fabs( reference[n][s] ) ); 
			}
		}
		std::cout << "Thomas verify: max |fpga - cpu| / max |cpu| ="; 
		for( int s = 0; s < num_substrates; s++ )
		{ std::cout << " " << density_names[s] << ":" << ( max_value[s] > 0.0 ? max_diff[s] / max_value[s] : max_diff[s] ); }
		std::cout << std::endl; 
	}
}

// Grid lines along `direction` -> accelerator rows
void Microenvironment::thomas_gather( int direction )
{
	const Sweep_Lines grid( mesh, direction ); 
	const int num_substrates = number_of_densities(); 

	#pragma omp parallel
	{
		std::vector<float*> rows( num_substrates ); 
		#pragma omp for
		for( long long line = 0; line < (long long) grid.lines; line++ )
		{
			for( int s = 0; s < num_substrates; s++ )
			{ rows[s] = thomas_fpga.rhs( direction, s, line ); }
			size_t voxel = grid.first_voxel( line ); 
			for( size_t p = 0; p < grid.points; p++, voxel += grid.stride )
			{
				const std::vector<double>& density = (*p_density_vectors)[voxel]; 
				for( int s = 0; s < num_substrates; s++ )
				{ rows[s][p] = static_cast<float>( density[s] ); }
			}
		}
	}
}

// Accelerator solutions -> grid lines along `direction`
void Microenvironment::thomas_scatter( int direction )
{
	const Sweep_Lines grid( mesh, direction ); 
	const int num_substrates = number_of_densities(); 

	#pragma omp parallel
	{
		std::vector<const float*> rows( num_substrates ); 
		#pragma omp for
		for( long long line = 0; line < (long long) grid.lines; line++ )
		{
			for( int s = 0; s < num_substrates; s++ )
			{ rows[s] = thomas_fpga.solution( direction, s, line ); }
			size_t voxel = grid.first_voxel( line ); 
			for( size_t p = 0; p < grid.points; p++, voxel += grid.stride )
			{
				std::vector<double>& density = (*p_density_vectors)[voxel]; 
				for( int s = 0; s < num_substrates; s++ )
				{ density[s] = rows[s][p]; }
			}
		}
	}
}

// BioFVM's original 3-D LOD step (double precision), using the LOD_3D setup
void Microenvironment::lod_3d_cpu( void )
{
	const int nx = mesh.x_coordinates.size(); 
	const int ny = mesh.y_coordinates.size(); 
	const int nz = mesh.z_coordinates.size(); 

	// x-diffusion 
	apply_dirichlet_conditions();
	#pragma omp parallel for 
	for( int k = 0; k < nz; k++ )
	{
		for( int j = 0; j < ny; j++ )
		{
			int n = voxel_index( 0, j, k ); 
			(*p_density_vectors)[n] /= thomas_denomx[0]; 
			for( int i = 1; i < nx; i++ )
			{
				n = voxel_index( i, j, k ); 
				axpy( &(*p_density_vectors)[n], thomas_constant1, (*p_density_vectors)[n - thomas_i_jump] ); 
				(*p_density_vectors)[n] /= thomas_denomx[i]; 
			}
			for( int i = nx - 2; i >= 0; i-- )
			{
				n = voxel_index( i, j, k ); 
				naxpy( &(*p_density_vectors)[n], thomas_cx[i], (*p_density_vectors)[n + thomas_i_jump] ); 
			}
		}
	}

	// y-diffusion 
	apply_dirichlet_conditions();
	#pragma omp parallel for 
	for( int k = 0; k < nz; k++ )
	{
		for( int i = 0; i < nx; i++ )
		{
			int n = voxel_index( i, 0, k ); 
			(*p_density_vectors)[n] /= thomas_denomy[0]; 
			for( int j = 1; j < ny; j++ )
			{
				n = voxel_index( i, j, k ); 
				axpy( &(*p_density_vectors)[n], thomas_constant1, (*p_density_vectors)[n - thomas_j_jump] ); 
				(*p_density_vectors)[n] /= thomas_denomy[j]; 
			}
			for( int j = ny - 2; j >= 0; j-- )
			{
				n = voxel_index( i, j, k ); 
				naxpy( &(*p_density_vectors)[n], thomas_cy[j], (*p_density_vectors)[n + thomas_j_jump] ); 
			}
		}
	}

	// z-diffusion 
	apply_dirichlet_conditions();
	#pragma omp parallel for 
	for( int j = 0; j < ny; j++ )
	{
		for( int i = 0; i < nx; i++ )
		{
			int n = voxel_index( i, j, 0 ); 
			(*p_density_vectors)[n] /= thomas_denomz[0]; 
			for( int k = 1; k < nz; k++ )
			{
				n = voxel_index( i, j, k ); 
				axpy( &(*p_density_vectors)[n], thomas_constant1, (*p_density_vectors)[n - thomas_k_jump] ); 
				(*p_density_vectors)[n] /= thomas_denomz[k]; 
			}
			for( int k = nz - 2; k >= 0; k-- )
			{
				n = voxel_index( i, j, k ); 
				naxpy( &(*p_density_vectors)[n], thomas_cz[k], (*p_density_vectors)[n + thomas_k_jump] ); 
			}
		}
	}

	apply_dirichlet_conditions();
}

void diffusion_decay_solver__constant_coefficients_LOD_3D( Microenvironment& M, double dt )
{
	if( M.mesh.regular_mesh == false || M.mesh.Cartesian_mesh == false )
	{
		std::cout << "Error: This algorithm is written for regular Cartesian meshes. Try: other solvers!" << std::endl << std::endl; 
		return; 
	}

	// define constants and pre-computed quantities 
	
	if( !M.diffusion_solver_setup_done )
	{
		std::cout << std::endl << "Using method " << __FUNCTION__ << " (implicit 3-D LOD with Thomas Algorithm) ... " 
		<< std::endl << std::endl;  
		
		M.thomas_denomx.resize( M.mesh.x_coordinates.size() , M.zero );
		M.thomas_cx.resize( M.mesh.x_coordinates.size() , M.zero );

		M.thomas_denomy.resize( M.mesh.y_coordinates.size() , M.zero );
		M.thomas_cy.resize( M.mesh.y_coordinates.size() , M.zero );
		
		M.thomas_denomz.resize( M.mesh.z_coordinates.size() , M.zero );
		M.thomas_cz.resize( M.mesh.z_coordinates.size() , M.zero );

		M.thomas_i_jump = 1; 
		M.thomas_j_jump = M.mesh.x_coordinates.size(); 
		M.thomas_k_jump = M.thomas_j_jump * M.mesh.y_coordinates.size(); 

		M.thomas_constant1 =  M.diffusion_coefficients; // dt*D/dx^2 
		M.thomas_constant1a = M.zero; // -dt*D/dx^2; 
		M.thomas_constant2 =  M.decay_rates; // (1/3)* dt*lambda 
		M.thomas_constant3 = M.one; // 1 + 2*constant1 + constant2; 
		M.thomas_constant3a = M.one; // 1 + constant1 + constant2; 		
			
		M.thomas_constant1 *= dt; 
		M.thomas_constant1 /= M.mesh.dx; 
		M.thomas_constant1 /= M.mesh.dx; 

		M.thomas_constant1a = M.thomas_constant1; 
		M.thomas_constant1a *= -1.0; 

		M.thomas_constant2 *= dt; 
		M.thomas_constant2 /= 3.0; // for the LOD splitting of the source 

		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant2; 

		M.thomas_constant3a += M.thomas_constant1; 
		M.thomas_constant3a += M.thomas_constant2; 

		// Thomas solver coefficients 

		M.thomas_cx.assign( M.mesh.x_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomx.assign( M.mesh.x_coordinates.size()  , M.thomas_constant3 ); 
		M.thomas_denomx[0] = M.thomas_constant3a; 
		M.thomas_denomx[ M.mesh.x_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.x_coordinates.size() == 1 )
		{ M.thomas_denomx[0] = M.one; M.thomas_denomx[0] += M.thomas_constant2; } 

		M.thomas_cx[0] /= M.thomas_denomx[0]; 
		for( unsigned int i=1 ; i <= M.mesh.x_coordinates.size()-1 ; i++ )
		{ 
			axpy( &M.thomas_denomx[i] , M.thomas_constant1 , M.thomas_cx[i-1] ); 
			M.thomas_cx[i] /= M.thomas_denomx[i]; // the value at  size-1 is not actually used  
		}

		M.thomas_cy.assign( M.mesh.y_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomy.assign( M.mesh.y_coordinates.size()  , M.thomas_constant3 ); 
		M.thomas_denomy[0] = M.thomas_constant3a; 
		M.thomas_denomy[ M.mesh.y_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.y_coordinates.size() == 1 )
		{ M.thomas_denomy[0] = M.one; M.thomas_denomy[0] += M.thomas_constant2; } 

		M.thomas_cy[0] /= M.thomas_denomy[0]; 
		for( unsigned int i=1 ; i <= M.mesh.y_coordinates.size()-1 ; i++ )
		{ 
			axpy( &M.thomas_denomy[i] , M.thomas_constant1 , M.thomas_cy[i-1] ); 
			M.thomas_cy[i] /= M.thomas_denomy[i]; // the value at  size-1 is not actually used  
		}

		M.thomas_cz.assign( M.mesh.z_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomz.assign( M.mesh.z_coordinates.size()  , M.thomas_constant3 ); 
		M.thomas_denomz[0] = M.thomas_constant3a; 
		M.thomas_denomz[ M.mesh.z_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.z_coordinates.size() == 1 )
		{ M.thomas_denomz[0] = M.one; M.thomas_denomz[0] += M.thomas_constant2; } 

		M.thomas_cz[0] /= M.thomas_denomz[0]; 
		for( unsigned int i=1 ; i <= M.mesh.z_coordinates.size()-1 ; i++ )
		{ 
			axpy( &M.thomas_denomz[i] , M.thomas_constant1 , M.thomas_cz[i-1] ); 
			M.thomas_cz[i] /= M.thomas_denomz[i]; // the value at  size-1 is not actually used  
		}	

		M.diffusion_solver_setup_done = true; 
	}

	// x-diffusion 

	M.compute();
	
	// M.apply_dirichlet_conditions();
	// M.thomas_solver_x();

	// // y-diffusion 

	// M.apply_dirichlet_conditions();
	// M.thomas_solver_y();

    // // z-diffusion 

	// M.apply_dirichlet_conditions();
	// M.thomas_solver_z();
 
	// M.apply_dirichlet_conditions();
	
	// reset gradient vectors 
    //	M.reset_all_gradient_vectors(); 

	return; 
}

void diffusion_decay_solver__constant_coefficients_LOD_2D( Microenvironment& M, double dt )
{
	if( M.mesh.regular_mesh == false )
	{
		std::cout << "Error: This algorithm is written for regular Cartesian meshes. Try: something else." << std::endl << std::endl; 
		return; 
	}
	
	// constants for the linear solver (Thomas algorithm) 
	
	if( !M.diffusion_solver_setup_done )
	{
		std::cout << std::endl << "Using method " << __FUNCTION__ << " (2D LOD with Thomas Algorithm) ... " << std::endl << std::endl;  
		
		M.thomas_denomx.resize( M.mesh.x_coordinates.size() , M.zero );
		M.thomas_cx.resize( M.mesh.x_coordinates.size() , M.zero );

		M.thomas_denomy.resize( M.mesh.y_coordinates.size() , M.zero );
		M.thomas_cy.resize( M.mesh.y_coordinates.size() , M.zero );
		
		// define constants and pre-computed quantities 

		M.thomas_i_jump = 1; 
		M.thomas_j_jump = M.mesh.x_coordinates.size(); 

		M.thomas_constant1 =  M.diffusion_coefficients; //   dt*D/dx^2 
		M.thomas_constant1a = M.zero; // -dt*D/dx^2; 
		M.thomas_constant2 =  M.decay_rates; // (1/2)*dt*lambda 
		M.thomas_constant3 = M.one; // 1 + 2*constant1 + constant2; 
		M.thomas_constant3a = M.one; // 1 + constant1 + constant2; 
		
		M.thomas_constant1 *= dt; 
		M.thomas_constant1 /= M.mesh.dx; 
		M.thomas_constant1 /= M.mesh.dx; 

		M.thomas_constant1a = M.thomas_constant1; 
		M.thomas_constant1a *= -1.0; 

		M.thomas_constant2 *= dt; 
		M.thomas_constant2 *= 0.5; // for splitting via LOD

		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant2; 

		M.thomas_constant3a += M.thomas_constant1; 
		M.thomas_constant3a += M.thomas_constant2; 
		
		// Thomas solver coefficients 

		M.thomas_cx.assign( M.mesh.x_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomx.assign( M.mesh.x_coordinates.size()  , M.thomas_constant3 ); 
		M.thomas_denomx[0] = M.thomas_constant3a; 
		M.thomas_denomx[ M.mesh.x_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.x_coordinates.size() == 1 )
		{ M.thomas_denomx[0] = M.one; M.thomas_denomx[0] += M.thomas_constant2; } 

		M.thomas_cx[0] /= M.thomas_denomx[0]; 
		for( unsigned int i=1 ; i <= M.mesh.x_coordinates.size()-1 ; i++ )
		{ 
			axpy( &M.thomas_denomx[i] , M.thomas_constant1 , M.thomas_cx[i-1] ); 
			M.thomas_cx[i] /= M.thomas_denomx[i]; // the value at  size-1 is not actually used  
		}

		M.thomas_cy.assign( M.mesh.y_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomy.assign( M.mesh.y_coordinates.size()  , M.thomas_constant3 ); 
		M.thomas_denomy[0] = M.thomas_constant3a; 
		M.thomas_denomy[ M.mesh.y_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.y_coordinates.size() == 1 )
		{ M.thomas_denomy[0] = M.one; M.thomas_denomy[0] += M.thomas_constant2; } 

		M.thomas_cy[0] /= M.thomas_denomy[0]; 
		for( unsigned int i=1 ; i <= M.mesh.y_coordinates.size()-1 ; i++ )
		{ 
			axpy( &M.thomas_denomy[i] , M.thomas_constant1 , M.thomas_cy[i-1] ); 
			M.thomas_cy[i] /= M.thomas_denomy[i]; // the value at  size-1 is not actually used  
		}

		M.diffusion_solver_setup_done = true; 
	}

	// set the pointer
	
	M.apply_dirichlet_conditions();

	// x-diffusion 
	#pragma omp parallel for 
	for( unsigned int j=0; j < M.mesh.y_coordinates.size() ; j++ )
	{
		// Thomas solver, x-direction

		// remaining part of forward elimination, using pre-computed quantities 
		unsigned int n = M.voxel_index(0,j,0);
		(*M.p_density_vectors)[n] /= M.thomas_denomx[0]; 

		n += M.thomas_i_jump; 
		for( unsigned int i=1; i < M.mesh.x_coordinates.size() ; i++ )
		{
			axpy( &(*M.p_density_vectors)[n] , M.thomas_constant1 , (*M.p_density_vectors)[n-M.thomas_i_jump] ); 
			(*M.p_density_vectors)[n] /= M.thomas_denomx[i]; 
			n += M.thomas_i_jump; 
		}

		// back substitution 
		n = M.voxel_index( M.mesh.x_coordinates.size()-2 ,j,0); 

		for( int i = M.mesh.x_coordinates.size()-2 ; i >= 0 ; i-- )
		{
			naxpy( &(*M.p_density_vectors)[n] , M.thomas_cx[i] , (*M.p_density_vectors)[n+M.thomas_i_jump] ); 
			n -= M.thomas_i_jump; 
		}
	}

	// y-diffusion 

	M.apply_dirichlet_conditions();
	#pragma omp parallel for 
	for( unsigned int i=0; i < M.mesh.x_coordinates.size() ; i++ )
	{
		// Thomas solver, y-direction

		// remaining part of forward elimination, using pre-computed quantities 

		int n = M.voxel_index(i,0,0);
		(*M.p_density_vectors)[n] /= M.thomas_denomy[0]; 

		n += M.thomas_j_jump; 
		for( unsigned int j=1; j < M.mesh.y_coordinates.size() ; j++ )
		{
			axpy( &(*M.p_density_vectors)[n] , M.thomas_constant1 , (*M.p_density_vectors)[n-M.thomas_j_jump] ); 
			(*M.p_density_vectors)[n] /= M.thomas_denomy[j]; 
			n += M.thomas_j_jump; 
		}

		// back substitution 
		n = M.voxel_index( i,M.mesh.y_coordinates.size()-2, 0); 

		for( int j = M.mesh.y_coordinates.size()-2 ; j >= 0 ; j-- )
		{
			naxpy( &(*M.p_density_vectors)[n] , M.thomas_cy[j] , (*M.p_density_vectors)[n+M.thomas_j_jump] ); 
			n -= M.thomas_j_jump; 
		}
	}

	M.apply_dirichlet_conditions();
	
	// reset gradient vectors 
//	M.reset_all_gradient_vectors(); 
	
	return; 
}

void diffusion_decay_explicit_uniform_rates( Microenvironment& M, double dt )
{
	using std::vector; 
	using std::cout; 
	using std::endl; 

	// static int n_jump_i = 1; 
	// static int n_jump_j = M.mesh.x_coordinates.size(); 
	// static int n_jump_k = M.mesh.x_coordinates.size() * M.mesh.y_coordinates.size(); 

	if( !M.diffusion_solver_setup_done )
	{	
		M.thomas_i_jump = 1; 
		M.thomas_j_jump = M.mesh.x_coordinates.size(); 
		M.thomas_k_jump = M.thomas_j_jump * M.mesh.y_coordinates.size(); 
	
		M.diffusion_solver_setup_done = true; 
	}
	
	if( M.mesh.uniform_mesh == false )
	{
		cout << "Error: This algorithm is written for uniform Cartesian meshes. Try: something else" << endl << endl; 
		return; 
	}

	// double buffering to reduce memory copy / allocation overhead 

	static vector< vector<double> >* pNew = &(M.temporary_density_vectors1);
	static vector< vector<double> >* pOld = &(M.temporary_density_vectors2);

	// swap the buffers 

	vector< vector<double> >* pTemp = pNew; 
	pNew = pOld; 
	pOld = pTemp; 
	M.p_density_vectors = pNew; 

	// static bool reaction_diffusion_shortcuts_are_set = false; 

	static vector<double> constant1 = (1.0 / ( M.mesh.dx * M.mesh.dx )) * M.diffusion_coefficients; 
	static vector<double> constant2 = dt * constant1; 
	static vector<double> constant3 = M.one + dt * M.decay_rates;

	static vector<double> constant4 = M.one - dt * M.decay_rates;

	#pragma omp parallel for
	for( unsigned int i=0; i < (*(M.p_density_vectors)).size() ; i++ )
	{
		unsigned int number_of_neighbors = M.mesh.connected_voxel_indices[i].size(); 

		double d1 = -1.0 * number_of_neighbors; 

		(*pNew)[i] = (*pOld)[i];  
		(*pNew)[i] *= constant4; 

		for( unsigned int j=0; j < number_of_neighbors ; j++ )
		{
			axpy( &(*pNew)[i], constant2, (*pOld)[  M.mesh.connected_voxel_indices[i][j] ] ); 
		}
		vector<double> temp = constant2; 
		temp *= d1; 
		axpy( &(*pNew)[i] , temp , (*pOld)[i] ); 
	}
	
	// reset gradient vectors 
//	M.reset_all_gradient_vectors(); 

	return; 
}

void diffusion_decay_solver__constant_coefficients_LOD_1D( Microenvironment& M, double dt )
{
	if( M.mesh.regular_mesh == false )
	{
		std::cout << "Error: This algorithm is written for regular Cartesian meshes. Try: something else." << std::endl << std::endl; 
		return; 
	}
	
	// constants for the linear solver (Thomas algorithm) 
	
	if( !M.diffusion_solver_setup_done )
	{
		std::cout << std::endl << "Using method " << __FUNCTION__ << " (2D LOD with Thomas Algorithm) ... " << std::endl << std::endl;  
		
		M.thomas_denomx.resize( M.mesh.x_coordinates.size() , M.zero );
		M.thomas_cx.resize( M.mesh.x_coordinates.size() , M.zero );

		// define constants and pre-computed quantities 

		M.thomas_i_jump = 1; 
		M.thomas_j_jump = M.mesh.x_coordinates.size(); 

		M.thomas_constant1 =  M.diffusion_coefficients; //   dt*D/dx^2 
		M.thomas_constant1a = M.zero; // -dt*D/dx^2; 
		M.thomas_constant2 =  M.decay_rates; // (1/2)*dt*lambda 
		M.thomas_constant3 = M.one; // 1 + 2*constant1 + constant2; 
		M.thomas_constant3a = M.one; // 1 + constant1 + constant2; 
		
		M.thomas_constant1 *= dt; 
		M.thomas_constant1 /= M.mesh.dx; 
		M.thomas_constant1 /= M.mesh.dx; 

		M.thomas_constant1a = M.thomas_constant1; 
		M.thomas_constant1a *= -1.0; 

		M.thomas_constant2 *= dt; 
		M.thomas_constant2 *= 1; // no splitting via LOD

		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant1; 
		M.thomas_constant3 += M.thomas_constant2; 

		M.thomas_constant3a += M.thomas_constant1; 
		M.thomas_constant3a += M.thomas_constant2; 
		
		// Thomas solver coefficients 

		M.thomas_cx.assign( M.mesh.x_coordinates.size() , M.thomas_constant1a ); 
		M.thomas_denomx.assign( M.mesh.x_coordinates.size()  , M.thomas_constant3 ); 
		M.thomas_denomx[0] = M.thomas_constant3a; 
		M.thomas_denomx[ M.mesh.x_coordinates.size()-1 ] = M.thomas_constant3a; 
		if( M.mesh.x_coordinates.size() == 1 )
		{ M.thomas_denomx[0] = M.one; M.thomas_denomx[0] += M.thomas_constant2; } 

		M.thomas_cx[0] /= M.thomas_denomx[0]; 
		for( unsigned int i=1 ; i <= M.mesh.x_coordinates.size()-1 ; i++ )
		{ 
			axpy( &M.thomas_denomx[i] , M.thomas_constant1 , M.thomas_cx[i-1] ); 
			M.thomas_cx[i] /= M.thomas_denomx[i]; // the value at  size-1 is not actually used  
		}

		M.diffusion_solver_setup_done = true; 
	}

	// set the pointer
	
	M.apply_dirichlet_conditions();

	// x-diffusion 
	#pragma omp parallel for 
	for( unsigned int j=0; j < M.mesh.y_coordinates.size() ; j++ )
	{
		// Thomas solver, x-direction

		// remaining part of forward elimination, using pre-computed quantities 
		unsigned int n = M.voxel_index(0,j,0);
		(*M.p_density_vectors)[n] /= M.thomas_denomx[0]; 

		n += M.thomas_i_jump; 
		for( unsigned int i=1; i < M.mesh.x_coordinates.size() ; i++ )
		{
			axpy( &(*M.p_density_vectors)[n] , M.thomas_constant1 , (*M.p_density_vectors)[n-M.thomas_i_jump] ); 
			(*M.p_density_vectors)[n] /= M.thomas_denomx[i]; 
			n += M.thomas_i_jump; 
		}

		// back substitution 
		n = M.voxel_index( M.mesh.x_coordinates.size()-2 ,j,0); 

		for( int i = M.mesh.x_coordinates.size()-2 ; i >= 0 ; i-- )
		{
			naxpy( &(*M.p_density_vectors)[n] , M.thomas_cx[i] , (*M.p_density_vectors)[n+M.thomas_i_jump] ); 
			n -= M.thomas_i_jump; 
		}
	}

	M.apply_dirichlet_conditions();
	
	// reset gradient vectors 
//	M.reset_all_gradient_vectors(); 
	
	return; 
}


};
