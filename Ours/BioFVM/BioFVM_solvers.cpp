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
#include "../vitis_common.h"
#include <CL/cl_ext_xilinx.h>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <omp.h>

// TODO: Remove these def from here
#define NUM_PE 8
#define LOG_NUM_PE 3

#define GROUP_SIZE 32
#define LOG_GROUP_SIZE 5
#define NUM_PE_GROUP (NUM_PE * GROUP_SIZE)
#define LOG_NUM_PE_GROUP (LOG_NUM_PE + LOG_GROUP_SIZE)

#define DEBUG_PRINT 1
#define SW_EMU_HOTFIX 0

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

static unsigned int ddr_bank_flag(int cu) {
    switch (cu) {
        case 0: return XCL_MEM_DDR_BANK0;
        case 1: return XCL_MEM_DDR_BANK1;
        case 2: return XCL_MEM_DDR_BANK2;
        case 3: return XCL_MEM_DDR_BANK3;
        default: return XCL_MEM_DDR_BANK0;
    }
}

cl::Buffer Microenvironment::make_bank_buffer(
    int cu,
    cl_mem_flags flags,
    size_t size,
    cl_int* err
) {
    cl_mem_ext_ptr_t ext;
    ext.obj = nullptr;
    ext.param = 0;
    ext.flags = ddr_bank_flag(cu_to_slr_mapping[cu]);

    return cl::Buffer(
        *context,
        flags | CL_MEM_EXT_PTR_XILINX,
        size,
        &ext,
        err
    );
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

void Microenvironment::compute() {
#if DEBUG_PRINT
	auto start = std::chrono::high_resolution_clock::now();
	std::cout << "Start compute: ";
	print_timestamp();
#endif
	uint32_t reload_mask = 0;
	bool is_coord_size_changed = false;
	cl_int err = CL_SUCCESS;
	
	// Dynamic vector to add the buffers that need to be migrated to the FPGA
	std::vector<cl::Memory> migrate_to_fpga;

	// Grid dimensions
	int z_coord_size = mesh.z_coordinates.size();
	int y_coord_size = mesh.y_coordinates.size();
	int x_coord_size = mesh.x_coordinates.size();
	int num_components = thomas_constant1.size();
	int comp_offset[NUM_CU];

	// If the components have changed, then compute the new distribution across the kernels
	if (num_components != curr_num_components) {
		curr_num_components = num_components;

		int base = curr_num_components / NUM_CU;
		int remainder = curr_num_components % NUM_CU;

		for (int cu = 0; cu < NUM_CU; cu++) {
            cu_buf[cu].cu_num_components = base + (cu < remainder ? 1 : 0);

			OCL_CHECK(err, err = krnl_compute[cu].setArg(15, cu_buf[cu].cu_num_components));
        }
	}

	// Set the real X and Y coord size to each kernels
	int current_comp_idx = 0;
	for (int cu = 0; cu < NUM_CU; cu++) {
		OCL_CHECK(err, err = krnl_compute[cu].setArg(17, x_coord_size));
		OCL_CHECK(err, err = krnl_compute[cu].setArg(18, y_coord_size));

		comp_offset[cu] = current_comp_idx;
        current_comp_idx += cu_buf[cu].cu_num_components;
	}

	// z-axis
	if (z_coord_size != curr_z_coord_size) {
		is_coord_size_changed = true;
		curr_z_coord_size = z_coord_size;

		// Assign to each kernel the coordinate Z
		for (int cu = 0; cu < NUM_CU; cu++) {
			OCL_CHECK(err, err = krnl_compute[cu].setArg(12, curr_z_coord_size)); // Scalar
		}
	}

	// y-axis
	if (y_coord_size != curr_y_coord_size) {
		is_coord_size_changed = true;
		curr_y_coord_size = y_coord_size;

		// The total number of rows is a multiple of NUM_PE * GROUP_SIZE
		curr_y_coord_padded_size = ((y_coord_size + NUM_PE_GROUP - 1) / NUM_PE_GROUP) * NUM_PE_GROUP;

		// Assign to each kernel the PADDED coordinate Y
		for (int cu = 0; cu < NUM_CU; cu++) {
			OCL_CHECK(err, err = krnl_compute[cu].setArg(13, curr_y_coord_padded_size)); // Scalar
		}
	}

	// x-axis
	if (x_coord_size != curr_x_coord_size) {
		is_coord_size_changed = true;
		curr_x_coord_size = x_coord_size;

		// The total number of rows is a multiple of NUM_PE * GROUP_SIZE
		curr_x_coord_padded_size = ((x_coord_size + NUM_PE_GROUP - 1) / NUM_PE_GROUP) * NUM_PE_GROUP;

		// Assign to each kernel the PADDED coordinate X
		for (int cu = 0; cu < NUM_CU; cu++) {
			OCL_CHECK(err, err = krnl_compute[cu].setArg(14, curr_x_coord_padded_size)); // Scalar
		}
	}

	// init_c1 and init_c2 are equal across all the components, therefore we can simply them for each kernel
	// We need 4 buffers, 4 pointer and 4 setArgs
	if (is_coord_size_changed) {
		reload_mask |= RELOAD_INIT;
		
		size_t size_init = 2 * NUM_PE * GROUP_SIZE * sizeof(int);

		for (int cu = 0; cu < NUM_CU; cu++) {
			OCL_CHECK(err, cu_buf[cu].buffer_init_c1 = make_bank_buffer(cu, CL_MEM_READ_ONLY, size_init, &err));
            OCL_CHECK(err, cu_buf[cu].buffer_init_c2 = make_bank_buffer(cu, CL_MEM_READ_ONLY, size_init, &err));
			
			OCL_CHECK(err, err = krnl_compute[cu].setArg(19, cu_buf[cu].buffer_init_c1));
			OCL_CHECK(err, err = krnl_compute[cu].setArg(20, cu_buf[cu].buffer_init_c2));

			// MAP
			OCL_CHECK(err, cu_buf[cu].ptr_init_c1 = (int*)q->enqueueMapBuffer(cu_buf[cu].buffer_init_c1, CL_TRUE, CL_MAP_WRITE, 0, size_init, NULL, NULL, &err));
            OCL_CHECK(err, cu_buf[cu].ptr_init_c2 = (int*)q->enqueueMapBuffer(cu_buf[cu].buffer_init_c2, CL_TRUE, CL_MAP_WRITE, 0, size_init, NULL, NULL, &err));

			// Pack
			for (int pe = 0; pe < NUM_PE; pe++) {
				for (int g = 0; g < GROUP_SIZE; g++) {
					// Logical line
					int lane = g * NUM_PE + pe;

					// Where to save the value in the buffer
					int idx_x = 0 * NUM_PE * GROUP_SIZE + pe * GROUP_SIZE + g;
					int idx_yz = 1 * NUM_PE * GROUP_SIZE + pe * GROUP_SIZE + g;

					// Config 0: Sweep X.
					// The external row is identified by (y, z).
					cu_buf[cu].ptr_init_c1[idx_x] = lane % curr_y_coord_padded_size;
					cu_buf[cu].ptr_init_c2[idx_x] = lane / curr_y_coord_padded_size;

					// Config 1: Sweep Y/Z.
					// The external row starts from x, with second coordinate initially 0.
					cu_buf[cu].ptr_init_c1[idx_yz] = lane;
					cu_buf[cu].ptr_init_c2[idx_yz] = 0;
				}
			}

			// Unmap
			OCL_CHECK(err, err = q->enqueueUnmapMemObject(cu_buf[cu].buffer_init_c1, cu_buf[cu].ptr_init_c1));
            OCL_CHECK(err, err = q->enqueueUnmapMemObject(cu_buf[cu].buffer_init_c2, cu_buf[cu].ptr_init_c2));

            migrate_to_fpga.push_back(cu_buf[cu].buffer_init_c1);
            migrate_to_fpga.push_back(cu_buf[cu].buffer_init_c2);
		}
	}
    
	unsigned int num_padded_voxels = curr_x_coord_padded_size * curr_y_coord_padded_size * curr_z_coord_size;

#if DEBUG_PRINT
	auto end_setup = std::chrono::high_resolution_clock::now();
	auto duration_setup = std::chrono::duration_cast<std::chrono::microseconds>(end_setup - start);
	std::cout << "Setup time: " << duration_setup.count() << " us" << std::endl;
#endif

	// =========================================================
    // 1. START DIRICHLET
    // =========================================================
#if DEBUG_PRINT
	auto start_dirichlet = std::chrono::high_resolution_clock::now();
#endif
	bool update_dirichlet = false;
	unsigned int num_voxels = mesh.voxels.size();

	if (num_voxels != curr_num_voxels) {
		curr_num_voxels = num_voxels;
		update_dirichlet = true;

		for (int cu = 0; cu < NUM_CU; cu++) {
            OCL_CHECK(err, err = krnl_compute[cu].setArg(2, num_padded_voxels)); 
        }
	}

	// Dirichlet buffers
	size_t req_size_dir = sizeof(int) * curr_num_components;
    size_t req_size_dir_val = sizeof(real_t) * curr_num_components;

	if (req_size_dir > cap_size_dir_general || req_size_dir_val > cap_size_dir_val_general) {
        update_dirichlet = true;
        cap_size_dir_general = req_size_dir; 
        cap_size_dir_val_general = req_size_dir_val; 
    }

	// Generalized reload of all the CU
	if (update_dirichlet) {
		reload_mask |= RELOAD_DIRICHLET;
		int actual_value_vector_size = dirichlet_value_vectors[0].size();

		for (int cu = 0; cu < NUM_CU; cu++) {
			if (cu_buf[cu].cu_num_components == 0) continue; // Guard

			// Compute the dimension for the specific CU
			size_t cu_req_size_dir = sizeof(char) * cu_buf[cu].cu_num_components;
            size_t cu_req_size_dir_val = sizeof(real_t) * cu_buf[cu].cu_num_components;

			// Allocate the new buffers
			OCL_CHECK(err, cu_buf[cu].buffer_apply_dirichlet = make_bank_buffer(cu, CL_MEM_READ_ONLY, cu_req_size_dir, &err));
            OCL_CHECK(err, cu_buf[cu].buffer_dirichlet_value = make_bank_buffer(cu, CL_MEM_READ_ONLY, cu_req_size_dir_val, &err));

			// Set the arguments
			OCL_CHECK(err, err = krnl_compute[cu].setArg(0, cu_buf[cu].buffer_apply_dirichlet));
            OCL_CHECK(err, err = krnl_compute[cu].setArg(1, cu_buf[cu].buffer_dirichlet_value));

			// Map
			OCL_CHECK(err, cu_buf[cu].ptr_apply_dirichlet = (char*)q->enqueueMapBuffer(cu_buf[cu].buffer_apply_dirichlet, CL_TRUE, CL_MAP_WRITE, 0, cu_req_size_dir, NULL, NULL, &err));
            OCL_CHECK(err, cu_buf[cu].ptr_dirichlet_value = (real_t*)q->enqueueMapBuffer(cu_buf[cu].buffer_dirichlet_value, CL_TRUE, CL_MAP_WRITE, 0, cu_req_size_dir_val, NULL, NULL, &err));

			// Pack
			for (int local_j = 0; local_j < cu_buf[cu].cu_num_components; local_j++) {
                int global_j = comp_offset[cu] + local_j; // Mappatura sul componente logico globale

                if (global_j < actual_value_vector_size) {
                    cu_buf[cu].ptr_apply_dirichlet[local_j] = dirichlet_activation_vectors[0][global_j] ? 1 : 0;
                    cu_buf[cu].ptr_dirichlet_value[local_j] = (real_t)dirichlet_value_vectors[0][global_j];
                } else {
                    cu_buf[cu].ptr_apply_dirichlet[local_j] = 0;
                    cu_buf[cu].ptr_dirichlet_value[local_j] = (real_t)0.0;
                }
            }

			// UnMap
            OCL_CHECK(err, err = q->enqueueUnmapMemObject(cu_buf[cu].buffer_apply_dirichlet, cu_buf[cu].ptr_apply_dirichlet));
            OCL_CHECK(err, err = q->enqueueUnmapMemObject(cu_buf[cu].buffer_dirichlet_value, cu_buf[cu].ptr_dirichlet_value));

            // Inserimento nel vettore di migrazione globale
            migrate_to_fpga.push_back(cu_buf[cu].buffer_apply_dirichlet);
            migrate_to_fpga.push_back(cu_buf[cu].buffer_dirichlet_value);
		}
	}

#if DEBUG_PRINT
	auto end_dirichlet = std::chrono::high_resolution_clock::now();
	auto duration_dirichlet = std::chrono::duration_cast<std::chrono::microseconds>(end_dirichlet - start_dirichlet);
	std::cout << "Dirichlet setup time: " << duration_dirichlet.count() << " us" << std::endl;
#endif

	// =========================================================
    // 2. THOMAS CONSTANTS (Unified for X, Y, Z)
    // =========================================================
#if DEBUG_PRINT
	auto start_thomas = std::chrono::high_resolution_clock::now();
#endif
	bool update_const1 = false;
    size_t req_thomas_constant1_size = sizeof(real_t) * curr_num_components;

    if (req_thomas_constant1_size > cap_thomas_constant1_size_general) {
        update_const1 = true;
        cap_thomas_constant1_size_general = req_thomas_constant1_size;
    }

    if (update_const1) {
        reload_mask |= RELOAD_CONST1;
        
        for (int cu = 0; cu < NUM_CU; cu++) {
			if (cu_buf[cu].cu_num_components == 0) continue; // Guard

            size_t cu_req_thomas_constant1_size = sizeof(real_t) * cu_buf[cu].cu_num_components;

            OCL_CHECK(err, cu_buf[cu].buffer_thomas_constant1 = make_bank_buffer(cu, CL_MEM_READ_ONLY, cu_req_thomas_constant1_size, &err));
            
            // Single binding to Argument 9
            OCL_CHECK(err, err = krnl_compute[cu].setArg(9, cu_buf[cu].buffer_thomas_constant1));

            // Map and Pack
            OCL_CHECK(err, cu_buf[cu].ptr_thomas_constant1 = (real_t*)q->enqueueMapBuffer(cu_buf[cu].buffer_thomas_constant1, CL_TRUE, CL_MAP_WRITE, 0, cu_req_thomas_constant1_size, NULL, NULL, &err));

            for (int local_j = 0; local_j < cu_buf[cu].cu_num_components; local_j++) {
                int global_j = comp_offset[cu] + local_j;
                cu_buf[cu].ptr_thomas_constant1[local_j] = (real_t)thomas_constant1[global_j];
            }

            OCL_CHECK(err, err = q->enqueueUnmapMemObject(cu_buf[cu].buffer_thomas_constant1, cu_buf[cu].ptr_thomas_constant1));

            migrate_to_fpga.push_back(cu_buf[cu].buffer_thomas_constant1);
        }
    }

	// =========================================================
    // 3. START THOMAS X
    // =========================================================
	bool update_thomas_x = false;
	
	size_t req_flat_thomas_array_size_x = sizeof(real_t) * (curr_x_coord_padded_size * curr_num_components);

	if (req_flat_thomas_array_size_x > cap_flat_thomas_array_size_x_general) {
        update_thomas_x = true;
        cap_flat_thomas_array_size_x_general = req_flat_thomas_array_size_x;
    }

	// Generalized reload for the thomas_x arrays
	if (update_thomas_x) {
        reload_mask |= RELOAD_THOMAS_X;

		for (int cu = 0; cu < NUM_CU; cu++) {
			size_t cu_req_flat_thomas_size = sizeof(real_t) * (curr_x_coord_padded_size * cu_buf[cu].cu_num_components);
			
			// Allocate the buffers
			OCL_CHECK(err, cu_buf[cu].buffer_thomas_denom_x = make_bank_buffer(cu, CL_MEM_READ_ONLY, cu_req_flat_thomas_size, &err));
			OCL_CHECK(err, cu_buf[cu].buffer_thomas_c_x = make_bank_buffer(cu, CL_MEM_READ_ONLY, cu_req_flat_thomas_size, &err));

			// Set arguments
			OCL_CHECK(err, err = krnl_compute[cu].setArg(3, cu_buf[cu].buffer_thomas_denom_x));
			OCL_CHECK(err, err = krnl_compute[cu].setArg(6, cu_buf[cu].buffer_thomas_c_x));

			// Map
			OCL_CHECK(err, cu_buf[cu].ptr_thomas_denom_x = (real_t*)q->enqueueMapBuffer(cu_buf[cu].buffer_thomas_denom_x, CL_TRUE, CL_MAP_WRITE, 0, cu_req_flat_thomas_size, NULL, NULL, &err));
			OCL_CHECK(err, cu_buf[cu].ptr_thomas_c_x = (real_t*)q->enqueueMapBuffer(cu_buf[cu].buffer_thomas_c_x, CL_TRUE, CL_MAP_WRITE, 0, cu_req_flat_thomas_size, NULL, NULL, &err));

			// Pack: Substrate-major limited by the components of the specific CU
			for (int i = 0; i < curr_x_coord_padded_size; i++) {
				for (int local_j = 0; local_j < cu_buf[cu].cu_num_components; local_j++) {
					
					int global_j = comp_offset[cu] + local_j; // Componente effettivo di PhysiCell
					int flat_idx = local_j * curr_x_coord_padded_size + i; // Indice flatten locale per la CU

					if (i < curr_x_coord_size) {
						// Real values
						cu_buf[cu].ptr_thomas_denom_x[flat_idx] = (real_t)(1.0 / thomas_denomx[i][global_j]);
						cu_buf[cu].ptr_thomas_c_x[flat_idx] = (real_t)thomas_cx[i][global_j];
					} else {
						// Padded values
						cu_buf[cu].ptr_thomas_denom_x[flat_idx] = (real_t)0.0;
						cu_buf[cu].ptr_thomas_c_x[flat_idx] = (real_t)0.0;
					}
				}
			}

			// UnMap
			OCL_CHECK(err, err = q->enqueueUnmapMemObject(cu_buf[cu].buffer_thomas_denom_x, cu_buf[cu].ptr_thomas_denom_x));
			OCL_CHECK(err, err = q->enqueueUnmapMemObject(cu_buf[cu].buffer_thomas_c_x, cu_buf[cu].ptr_thomas_c_x));

			migrate_to_fpga.push_back(cu_buf[cu].buffer_thomas_denom_x);
			migrate_to_fpga.push_back(cu_buf[cu].buffer_thomas_c_x);
		}
	}

	// =========================================================
    // 4. START THOMAS Y
    // =========================================================
	bool update_thomas_y = false;

	size_t req_flat_thomas_array_size_y = sizeof(real_t) * (curr_y_coord_padded_size * curr_num_components);
	size_t req_thomas_constant1_size_y = sizeof(real_t) * (curr_num_components); 

	if (req_flat_thomas_array_size_y > cap_flat_thomas_array_size_y_general) {
        update_thomas_y = true;
        cap_flat_thomas_array_size_y_general = req_flat_thomas_array_size_y; 
    }

	if (update_thomas_y) {
        reload_mask |= RELOAD_THOMAS_Y;

		for (int cu = 0; cu < NUM_CU; cu++) {
			size_t cu_req_flat_thomas_size = sizeof(real_t) * (curr_y_coord_padded_size * cu_buf[cu].cu_num_components);

			// Allocate buffers
            OCL_CHECK(err, cu_buf[cu].buffer_thomas_denom_y = make_bank_buffer(cu, CL_MEM_READ_ONLY, cu_req_flat_thomas_size, &err));
            OCL_CHECK(err, cu_buf[cu].buffer_thomas_c_y = make_bank_buffer(cu, CL_MEM_READ_ONLY, cu_req_flat_thomas_size, &err));

			// Set argurments
            OCL_CHECK(err, err = krnl_compute[cu].setArg(4, cu_buf[cu].buffer_thomas_denom_y));
            OCL_CHECK(err, err = krnl_compute[cu].setArg(7, cu_buf[cu].buffer_thomas_c_y));

            // Map
            OCL_CHECK(err, cu_buf[cu].ptr_thomas_denom_y = (real_t*)q->enqueueMapBuffer(cu_buf[cu].buffer_thomas_denom_y, CL_TRUE, CL_MAP_WRITE, 0, cu_req_flat_thomas_size, NULL, NULL, &err));
            OCL_CHECK(err, cu_buf[cu].ptr_thomas_c_y = (real_t*)q->enqueueMapBuffer(cu_buf[cu].buffer_thomas_c_y, CL_TRUE, CL_MAP_WRITE, 0, cu_req_flat_thomas_size, NULL, NULL, &err));

			// Pack Substrate-major
            for (int i = 0; i < curr_y_coord_padded_size; i++) { 
                for (int local_j = 0; local_j < cu_buf[cu].cu_num_components; local_j++) {
                    int global_j = comp_offset[cu] + local_j;
                    int flat_idx = local_j * curr_y_coord_padded_size + i; 

                    if (i < curr_y_coord_size) {
                        cu_buf[cu].ptr_thomas_denom_y[flat_idx] = (real_t)(1.0 / thomas_denomy[i][global_j]); 
                        cu_buf[cu].ptr_thomas_c_y[flat_idx] = (real_t)thomas_cy[i][global_j]; 
                    } else {
                        // Prevent forward sweep leakage
                        cu_buf[cu].ptr_thomas_denom_y[flat_idx] = (real_t)0.0;
                        cu_buf[cu].ptr_thomas_c_y[flat_idx] = (real_t)0.0;
                    }
                }
            }

			// Unmap
            OCL_CHECK(err, err = q->enqueueUnmapMemObject(cu_buf[cu].buffer_thomas_denom_y, cu_buf[cu].ptr_thomas_denom_y));
            OCL_CHECK(err, err = q->enqueueUnmapMemObject(cu_buf[cu].buffer_thomas_c_y, cu_buf[cu].ptr_thomas_c_y));

            migrate_to_fpga.push_back(cu_buf[cu].buffer_thomas_denom_y);
            migrate_to_fpga.push_back(cu_buf[cu].buffer_thomas_c_y);
		}
	}

	// =========================================================
    // 5. START THOMAS Z
    // =========================================================
	bool update_thomas_z = false;
	
	size_t req_flat_thomas_array_size_z = sizeof(real_t) * (curr_z_coord_size * curr_num_components);
    size_t req_thomas_constant1_size_z = sizeof(real_t) * curr_num_components;

	if (req_flat_thomas_array_size_z > cap_flat_thomas_array_size_z_general) {
        update_thomas_z = true;
        cap_flat_thomas_array_size_z_general = req_flat_thomas_array_size_z; 
    }

	if (update_thomas_z) {
        reload_mask |= RELOAD_THOMAS_Z;

        for (int cu = 0; cu < NUM_CU; cu++) {
            size_t cu_req_flat_thomas_size = sizeof(real_t) * (curr_z_coord_size * cu_buf[cu].cu_num_components);

            // Allocate buffers
            OCL_CHECK(err, cu_buf[cu].buffer_thomas_denom_z = make_bank_buffer(cu, CL_MEM_READ_ONLY, cu_req_flat_thomas_size, &err));
            OCL_CHECK(err, cu_buf[cu].buffer_thomas_c_z = make_bank_buffer(cu, CL_MEM_READ_ONLY, cu_req_flat_thomas_size, &err));

            OCL_CHECK(err, err = krnl_compute[cu].setArg(5, cu_buf[cu].buffer_thomas_denom_z));
            OCL_CHECK(err, err = krnl_compute[cu].setArg(8, cu_buf[cu].buffer_thomas_c_z));

            // Map
            OCL_CHECK(err, cu_buf[cu].ptr_thomas_denom_z = (real_t*)q->enqueueMapBuffer(cu_buf[cu].buffer_thomas_denom_z, CL_TRUE, CL_MAP_WRITE, 0, cu_req_flat_thomas_size, NULL, NULL, &err));
            OCL_CHECK(err, cu_buf[cu].ptr_thomas_c_z = (real_t*)q->enqueueMapBuffer(cu_buf[cu].buffer_thomas_c_z, CL_TRUE, CL_MAP_WRITE, 0, cu_req_flat_thomas_size, NULL, NULL, &err));

            // Pack Substrate-major (no padding needed for Z)
            for (int i = 0; i < curr_z_coord_size; i++) { 
                for (int local_j = 0; local_j < cu_buf[cu].cu_num_components; local_j++) {
                    int global_j = comp_offset[cu] + local_j;
                    int flat_idx = local_j * curr_z_coord_size + i; 

                    cu_buf[cu].ptr_thomas_denom_z[flat_idx] = (real_t)(1.0 / thomas_denomz[i][global_j]); 
                    cu_buf[cu].ptr_thomas_c_z[flat_idx] = (real_t)thomas_cz[i][global_j]; 
                }
            }

            // Unmap
            OCL_CHECK(err, err = q->enqueueUnmapMemObject(cu_buf[cu].buffer_thomas_denom_z, cu_buf[cu].ptr_thomas_denom_z));
            OCL_CHECK(err, err = q->enqueueUnmapMemObject(cu_buf[cu].buffer_thomas_c_z, cu_buf[cu].ptr_thomas_c_z));

            migrate_to_fpga.push_back(cu_buf[cu].buffer_thomas_denom_z);
            migrate_to_fpga.push_back(cu_buf[cu].buffer_thomas_c_z);
        }
    }

	#if DEBUG_PRINT
		auto end_thomas = std::chrono::high_resolution_clock::now();
		auto duration_thomas = std::chrono::duration_cast<std::chrono::microseconds>(end_thomas - start_thomas);
		std::cout << "Thomas setup time: " << duration_thomas.count() << " us" << std::endl;
	#endif

	// =========================================================
	// 5. START DENSITY (We always need to read it and send it to the FPGA)
	// =========================================================
	#if DEBUG_PRINT
		auto start_density = std::chrono::high_resolution_clock::now();
	#endif

	// For each CU, pack only its component slice directly into the CU-local FPGA buffer.
#pragma omp parallel for
	for (int cu = 0; cu < NUM_CU; cu++) {
		#if SW_EMU_HOTFIX
		reload_mask |= RELOAD_DIRICHLET;
		reload_mask |= RELOAD_THOMAS_X;
		reload_mask |= RELOAD_THOMAS_Y;
		reload_mask |= RELOAD_THOMAS_Z;
		reload_mask |= RELOAD_CONST1;
		// reload_mask |= RELOAD_INIT;
		#endif

		OCL_CHECK(err, err = krnl_compute[cu].setArg(16, reload_mask));

		if (cu_buf[cu].cu_num_components == 0) {
			continue;
		}

		size_t cu_density_size_bytes =
			sizeof(real_t) *
			static_cast<size_t>(num_padded_voxels) *
			static_cast<size_t>(cu_buf[cu].cu_num_components);

		// If the size has changed, re-allocate CU buffers.
		if (cu_density_size_bytes > cu_buf[cu].cap_flat_density_size) {
			cu_buf[cu].cap_flat_density_size =
				cu_density_size_bytes +
				cu_density_size_bytes / 2 +
				64;

			OCL_CHECK(err,
				cu_buf[cu].buffer_density_PING =
					make_bank_buffer(
						cu,
						CL_MEM_READ_WRITE,
						cu_buf[cu].cap_flat_density_size,
						&err
					)
			);

			OCL_CHECK(err,
				cu_buf[cu].buffer_density_PONG =
					make_bank_buffer(
						cu,
						CL_MEM_READ_WRITE,
						cu_buf[cu].cap_flat_density_size,
						&err
					)
			);

			OCL_CHECK(err,
				err = krnl_compute[cu].setArg(
					10,
					cu_buf[cu].buffer_density_PING
				)
			);

			OCL_CHECK(err,
				err = krnl_compute[cu].setArg(
					11,
					cu_buf[cu].buffer_density_PONG
				)
			);
		}

		OCL_CHECK(err,
			cu_buf[cu].ptr_density =
				(real_t*)q->enqueueMapBuffer(
					cu_buf[cu].buffer_density_PING,
					CL_TRUE,
					CL_MAP_WRITE,
					0,
					cu_density_size_bytes,
					NULL,
					NULL,
					&err
				)
		);

		pack_density_cu_slice(
			cu_buf[cu].ptr_density,
			curr_x_coord_size,
			curr_x_coord_padded_size,
			curr_y_coord_size,
			curr_y_coord_padded_size,
			curr_z_coord_size,
			static_cast<unsigned int>(comp_offset[cu]),
			static_cast<unsigned int>(cu_buf[cu].cu_num_components)
		);

		OCL_CHECK(err,
			err = q->enqueueUnmapMemObject(
				cu_buf[cu].buffer_density_PING,
				cu_buf[cu].ptr_density
			)
		);
	}

	for (int cu = 0; cu < NUM_CU; cu++) {
		migrate_to_fpga.push_back(cu_buf[cu].buffer_density_PING);
	}

	OCL_CHECK(err, q->finish());

	#if DEBUG_PRINT
		auto end_density = std::chrono::high_resolution_clock::now();
		auto duration_density = std::chrono::duration_cast<std::chrono::microseconds>(end_density - start_density);
		std::cout << "Density setup time: " << duration_density.count() << " us" << std::endl;
	#endif

	#if DEBUG_PRINT
		auto start_h2d = std::chrono::high_resolution_clock::now();
	#endif

	// Single migrate call for all buffers across all CUs.
	OCL_CHECK(err, err = q->enqueueMigrateMemObjects(migrate_to_fpga, 0));
	OCL_CHECK(err, q->finish());

	#if DEBUG_PRINT
		auto end_h2d = std::chrono::high_resolution_clock::now();
		auto duration_h2d = std::chrono::duration_cast<std::chrono::microseconds>(end_h2d - start_h2d);
		std::cout << "Host-to-Device transfer time: " << duration_h2d.count() << " us" << std::endl;
	#endif

	#if DEBUG_PRINT
		auto start_kernel = std::chrono::high_resolution_clock::now();
		std::cout << "Start kernel: ";
		print_timestamp();
	#endif

	// Parallel Kernel Execution.
	for (int cu = 0; cu < NUM_CU; cu++) {
		OCL_CHECK(err, err = q->enqueueTask(krnl_compute[cu]));
	}

	OCL_CHECK(err, q->finish());

	#if DEBUG_PRINT
		auto end_kernel = std::chrono::high_resolution_clock::now();
		auto duration_kernel = std::chrono::duration_cast<std::chrono::microseconds>(end_kernel - start_kernel);
		std::cout << "Kernel time: " << duration_kernel.count() << " us" << std::endl;
		std::cout << "End kernel: ";
		print_timestamp();
	#endif

	// =========================================================
	// READ BACK RESULTS FROM FPGA
	// =========================================================
	std::vector<cl::Memory> migrate_back_to_host;

	for (int cu = 0; cu < NUM_CU; cu++) {
		if (cu_buf[cu].cu_num_components == 0) {
			continue;
		}

		migrate_back_to_host.push_back(cu_buf[cu].buffer_density_PONG);
	}

	#if DEBUG_PRINT
		auto start_d2h = std::chrono::high_resolution_clock::now();
	#endif

	OCL_CHECK(err,
		err = q->enqueueMigrateMemObjects(
			migrate_back_to_host,
			CL_MIGRATE_MEM_OBJECT_HOST
		)
	);

	OCL_CHECK(err, q->finish());

	#if DEBUG_PRINT
		auto end_d2h = std::chrono::high_resolution_clock::now();
		auto duration_d2h = std::chrono::duration_cast<std::chrono::microseconds>(end_d2h - start_d2h);
		std::cout << "Device-to-Host transfer time: " << duration_d2h.count() << " us" << std::endl;
	#endif

	#if DEBUG_PRINT
		auto start_density_2 = std::chrono::high_resolution_clock::now();
	#endif

#pragma omp parallel for
	for (int cu = 0; cu < NUM_CU; cu++) {
		if (cu_buf[cu].cu_num_components == 0) {
			continue;
		}

		size_t cu_density_size_bytes =
			sizeof(real_t) *
			static_cast<size_t>(num_padded_voxels) *
			static_cast<size_t>(cu_buf[cu].cu_num_components);

		OCL_CHECK(err,
			cu_buf[cu].ptr_density =
				(real_t*)q->enqueueMapBuffer(
					cu_buf[cu].buffer_density_PONG,
					CL_TRUE,
					CL_MAP_READ,
					0,
					cu_density_size_bytes,
					NULL,
					NULL,
					&err
				)
		);

		unpack_density_cu_slice(
			cu_buf[cu].ptr_density,
			curr_x_coord_size,
			curr_x_coord_padded_size,
			curr_y_coord_size,
			curr_y_coord_padded_size,
			curr_z_coord_size,
			static_cast<unsigned int>(comp_offset[cu]),
			static_cast<unsigned int>(cu_buf[cu].cu_num_components)
		);

		OCL_CHECK(err,
			err = q->enqueueUnmapMemObject(
				cu_buf[cu].buffer_density_PONG,
				cu_buf[cu].ptr_density
			)
		);
	}

	OCL_CHECK(err, q->finish());

	#if DEBUG_PRINT
		auto end_density_2 = std::chrono::high_resolution_clock::now();
		auto duration_density_2 = std::chrono::duration_cast<std::chrono::microseconds>(end_density_2 - start_density_2);
		std::cout << "Density teardown time: " << duration_density_2.count() << " us" << std::endl;
	#endif

	#if DEBUG_PRINT
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
		std::cout << "Compute time: " << duration.count() << " us" << std::endl;
		std::cout << "End compute: ";
		print_timestamp();
	#endif
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
