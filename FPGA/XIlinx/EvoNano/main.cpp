/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

// 
//  EVO-NANO PHYSICELL 
// 
//  Author: Milos Savic (svc@dmi.uns.ac.rs)
//

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>
#include <chrono>
#include "vitis_common.h"

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h"

#include "./custom_modules/evo_nano_3D.h"
#include "./custom_modules/spexp.h"
#include "./custom_modules/csvlog.h"

using namespace BioFVM;
using namespace PhysiCell;


int count_cancer_cells() {
    int counter = 0;

    for (int i = 0; i < all_cells->size(); i++) {
        Cell* c = (*all_cells)[i];
        int type = c->type;
        bool dead = c->phenotype.death.dead;

        if (!dead && (type == CANCER_CELL_TYPE || type == CANCER_STEM_CELL_TYPE)) {
            counter++;
        }
    }    

    return counter;
}


void angiogen_report() {
    int angiogen_index = microenvironment.find_density_index("angiogen");
    if (angiogen_index == -1)
        return;
                
    double max_ang = -1000, min_ang = 1000;
    int max_ang_voxel[3], min_ang_voxel[3];
    int max_multiplicity, min_multiplicity;

    for (unsigned int i = 0; i < microenvironment.mesh.x_coordinates.size(); i++) {
        for (unsigned int j = 0; j < microenvironment.mesh.y_coordinates.size(); j++) {
            for (unsigned int k = 0; k < microenvironment.mesh.z_coordinates.size(); k++) {
                std::vector<double> den = microenvironment(i, j, k);
                double ang = den[angiogen_index];

                if (ang > max_ang) {
                    max_multiplicity = 1;
                    max_ang = ang;
                    max_ang_voxel[0] = i; max_ang_voxel[1] = j; max_ang_voxel[2] = k;
                } else if (ang == max_ang) {
                    max_multiplicity++;
                }

                if (ang < min_ang) {
                    min_multiplicity = 1;
                    min_ang = ang;
                    min_ang_voxel[0] = i; min_ang_voxel[1] = j; min_ang_voxel[2] = k;
                } else if (ang == min_ang) {
                    min_multiplicity++;
                }
            }
        }
    }

    std::cout << "Max ang = " << max_ang << ", multiplicity " << max_multiplicity 
              << " [" << max_ang_voxel[0] << "," << max_ang_voxel[1] << "," << max_ang_voxel[2] << "]" << std::endl;

    std::cout << "Min ang = " << min_ang << ", multiplicity " << min_multiplicity 
              << " [" << min_ang_voxel[0] << "," << min_ang_voxel[1] << "," << min_ang_voxel[2] << "]" << std::endl;

    /*
    for (int i = 0; i < all_cells->size(); i++) {
        Cell* c = (*all_cells)[i];
        double ang_inside = c->phenotype.molecular.internalized_total_substrates[angiogen_index];
        if (c->type == VASCULAR_CELL_TYPE) {
            int numexp_index = c->custom_data.find_variable_index("num_expansions");
            int num_expansions = c->custom_data[numexp_index];
            
            double near_ang = microenvironment.nearest_density_vector(c->position)[angiogen_index];
            std::cout << c->type_name << ", angiogen near = " << near_ang << ", num expansions = " << num_expansions << std::endl ;
        }
    }
    */
}

// Funzione custom per il test FPGA
void dump_simulation_state(Microenvironment& M, const std::string& filename) {
    std::ofstream outfile(filename);
    
    if (!outfile.is_open()) {
        std::cout << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    int total_voxels = M.mesh.voxels.size();
    
    for (int i = 0; i < total_voxels; i++) {
        int actual_size = M.density_vector(i).size();
        for (int j = 0; j < actual_size; j++) {
            outfile << M.density_vector(i)[j];
            if (j < actual_size - 1) outfile << ",";
        }
        outfile << "\n";
    }
    
    outfile.close();
    std::cout << ">>> STATO SIMULAZIONE SALVATO IN: " << filename << " <<<" << std::endl;
}

int main(int argc, char* argv[])
{
	/*
	=========================
	=== START BOILERPLATE ===
	=========================
	*/

	// TARGET_DEVICE macro needs to be passed from gcc command line
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " <xclbin> <xml (optional)>" << std::endl;
		return EXIT_FAILURE;
	}

	std::string xclbinFilename = argv[1];

	// Creates a vector of DATA_SIZE elements with an initial value of 10 and 32
	// using customized allocator for getting buffer alignment to 4k boundary

	std::vector<cl::Device> devices;
	cl_int err;
	std::vector<cl::Platform> platforms;
	bool found_device = false;

    // Create the global variables that will be used later on in the solver
    cl::Context context_object;
    cl::CommandQueue q_object;

    // kernels global variables: one cl::Kernel per TDMA_batch compute unit (ThomasVsPcr/Thomas xclbin)
    std::vector<cl::Kernel> krnl_tdma_objects;

    cl::Program program_object;

	// traversing all Platforms To find Xilinx Platform and targeted
	// Device in Xilinx Platform
	cl::Platform::get(&platforms);
	for (size_t i = 0; (i < platforms.size()) & (found_device == false); i++) {
		cl::Platform platform = platforms[i];
		std::string platformName = platform.getInfo<CL_PLATFORM_NAME>();
		if (platformName == "Xilinx") {
			devices.clear();
			platform.getDevices(CL_DEVICE_TYPE_ACCELERATOR, &devices);
			if (devices.size()) {
				found_device = true;
				break;
			}
		}
	}
	if (found_device == false) {
		std::cout << "Error: Unable to find Target Device " << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "INFO: Reading " << xclbinFilename << std::endl;
	FILE* fp;
	if ((fp = fopen(xclbinFilename.c_str(), "r")) == nullptr) {
		printf("ERROR: %s xclbin not available please build\n", xclbinFilename.c_str());
		exit(EXIT_FAILURE);
	}
	// Load xclbin
	std::cout << "Loading: '" << xclbinFilename << "'\n";
	std::ifstream bin_file(xclbinFilename, std::ifstream::binary);
	bin_file.seekg(0, bin_file.end);
	unsigned nb = bin_file.tellg();
	bin_file.seekg(0, bin_file.beg);
	char* buf = new char[nb];
	bin_file.read(buf, nb);

	// Creating Program from Binary File
	cl::Program::Binaries bins;
	bins.push_back({buf, nb});
	bool valid_device = false;
	for (unsigned int i = 0; i < devices.size(); i++) {
		auto device = devices[i];
		// Creating Context and Command Queue for selected Device
		OCL_CHECK(err, context_object = cl::Context(device, nullptr, nullptr, nullptr, &err));
		OCL_CHECK(err, q_object = cl::CommandQueue(context_object, device, CL_QUEUE_PROFILING_ENABLE | CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &err));
		std::cout << "Trying to program device[" << i << "]: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
		program_object = cl::Program(context_object, {device}, bins, nullptr, &err);
		if (err != CL_SUCCESS) {
			std::cout << "Failed to program device[" << i << "] with xclbin file!\n";
		} else {
			std::cout << "Device[" << i << "]: program successful!\n";

			// ThomasVsPcr/Thomas (TDMA_solver.cpp): TDMA_batch(a, b, c, d, u, M, N, B, iters),
			// CUs TDMA_batch_1..TDMA_batch_n. ADI3D_F32 also has a TDMA_batch kernel, with 13 arguments.
			cl_uint num_args = 0;
			cl::Kernel tdma_all_cus(program_object, "TDMA_batch", &err);
			if (err == CL_SUCCESS) {
				OCL_CHECK(err, err = tdma_all_cus.getInfo(CL_KERNEL_NUM_ARGS, &num_args));
			}
			if (num_args != 9) {
				std::cout << "Error: " << xclbinFilename << " does not contain the TDMA_batch kernel of FPGA/XIlinx/ThomasVsPcr/Thomas" << std::endl;
				exit(EXIT_FAILURE);
			}
			cl_uint num_cus = 0;
			OCL_CHECK(err, err = clGetKernelInfo(tdma_all_cus(), CL_KERNEL_COMPUTE_UNIT_COUNT, sizeof(num_cus), &num_cus, nullptr));
			for (cl_uint cu = 1; cu <= num_cus; cu++) {
				std::string cu_name = "TDMA_batch:{TDMA_batch_" + std::to_string(cu) + "}";
				cl::Kernel krnl;
				OCL_CHECK(err, krnl = cl::Kernel(program_object, cu_name.c_str(), &err));
				krnl_tdma_objects.push_back(krnl);
			}
			std::cout << "Found " << num_cus << " TDMA_batch compute units" << std::endl;

			valid_device = true;
			break; // we break because we found a valid device
		}
	}
	if (!valid_device) {
		std::cout << "Failed to program any device found, exit!\n";
		exit(EXIT_FAILURE);
	}

	/*
	=========================
	==== END BOILERPLATE ====
	=========================
	*/

    std::cout << "Physicell Evo-Nano experimental version by svc" << std::endl << std::endl;

    bool XML_status = false;
    if( argc > 2 ) {
        XML_status = load_PhysiCell_config_file( argv[2] );
    } else {
        // Proviamo il percorso standard
        XML_status = load_PhysiCell_config_file( "./config/PhysiCell_settings.xml" );
    }

    if( !XML_status ) {
        std::cerr << "ERRORE FATALE: Impossibile trovare il file di configurazione XML!" << std::endl;
        exit(-1);
    }

    // OpenMP setup
    omp_set_num_threads(PhysiCell_settings.omp_num_threads);

    // PNRG setup
    SeedRandom( 1UL );

    // time setup
    std::string time_units = "min";

    /* Microenvironment setup */
    setup_microenvironment(); // modify this in the custom code
    microenvironment.setup_opencl(&context_object, &q_object, krnl_tdma_objects);

    /* PhysiCell setup */

    // set mechanics voxel size, and match the data structure to BioFVM
    double mechanics_voxel_size = 30;
    Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );

    /* Users typically start modifying here. START USERMODS */
    create_cell_types();

    setup_tissue();
    /* Users typically stop modifying here. END USERMODS */

    // set MultiCellDS save options
    set_save_biofvm_mesh_as_matlab( true );
    set_save_biofvm_data_as_matlab( true );
    set_save_biofvm_cell_data( true );
    set_save_biofvm_cell_data_as_custom_matlab( true );

    // save a simulation snapshot

    char filename[1024];
    sprintf( filename , "%s/initial" , PhysiCell_settings.folder.c_str() );
    save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time );

    
    // save a quick SVG cross section through z = 0, after setting its
    // length bar to 200 microns

    PhysiCell_SVG_options.length_bar = 200;

    // for simplicity, set a pathology coloring function
    std::vector<std::string> (*cell_coloring_function)(Cell*) = my_coloring_function;

    sprintf( filename , "%s/initial.svg" , PhysiCell_settings.folder.c_str() );
    SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );

    display_citations();

    // set the performance timers

    BioFVM::RUNTIME_TIC();
    BioFVM::TIC();

    std::ofstream report_file;
    if( PhysiCell_settings.enable_legacy_saves == true )
    {
        sprintf( filename , "%s/simulation_report.txt" , PhysiCell_settings.folder.c_str() );

        report_file.open(filename); 	// create the data log file
        report_file<<"simulated time\tnum cells\tnum division\tnum death\twall time"<<std::endl;
    }

    CSV_log* csvlog = new CSV_log("tabela.csv"); 
    bool print_ang_report = parameters.bools("print_angiogen_report");

    std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

    
    int num_expansions = 0;
    bool expand_vascularity = parameters.bools("expand_vascularity_when_expanding_space");
    double expansion_scale_factor = 1;
    if (expand_vascularity) {
        expansion_scale_factor = parameters.doubles("expansion_scale_factor");
    }
    
    int cancer_cells_to_activate_NP = 0;
    if (NP_simul) {
        cancer_cells_to_activate_NP = parameters.ints("cancer_cells_to_activate_NP");
    }


    // main loop
    try
    {
        // TODO: REMOVE AFTER TESTING
        int test_iteration_count = 0; // Counter for testing

        while( PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1*diffusion_dt )
        {
            // save data if it's time.
            if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_full_save_time ) < 0.01 * diffusion_dt )
            {
                std::cout << "Saving data" << std::endl;
                microenvironment.dump_density_to_file(PhysiCell_settings.folder + "/density_dump_" + std::to_string(PhysiCell_globals.full_output_index) + ".txt");
                csvlog->log();

                display_simulation_status( std::cout );
                cells_report();
                if (NP_simul) {
                    std::cout << "NP_CC_deaths = " << NP_CC_deaths << ", ";
                    std::cout << "NP_HC_deaths = " << NP_HC_deaths << std::endl;
                }
                
                if (cytokine_caf_simul) {
                    std::cout << "Cytokine-based CC-->CSC conversions: " << cytokine_conversions << std::endl; 
                }
                
                if (prostaglandin_simul) {
                    std::cout << "Prostaglandin-based CC-->CSC conversions: " << prostaglandin_conversions << std::endl; 
                }
                
                if (print_ang_report) angiogen_report();
                std::cout << std::endl;

                if( PhysiCell_settings.enable_legacy_saves == true )
                {
                    log_output( PhysiCell_globals.current_time , PhysiCell_globals.full_output_index, microenvironment, report_file);
                }

                if( PhysiCell_settings.enable_full_saves == true )
                {
                    sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index );

                    save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time );
                }

                PhysiCell_globals.full_output_index++;
                PhysiCell_globals.next_full_save_time += PhysiCell_settings.full_save_interval;
            }

            // save SVG plot if it's time
            if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_SVG_save_time  ) < 0.01 * diffusion_dt )
            {
                std::cout << "Saving SVG" << std::endl;
                if( PhysiCell_settings.enable_SVG_saves == true )
                {
                    sprintf( filename , "%s/snapshot%08u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index );
                    SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );

                    PhysiCell_globals.SVG_output_index++;
                    PhysiCell_globals.next_SVG_save_time  += PhysiCell_settings.SVG_save_interval;
                }
            }

            // update the microenvironment
            test_iteration_count++;
            //std::cout << "Simulating microenvironment for the " << test_iteration_count << "th time" << std::endl;
            microenvironment.simulate_diffusion_decay(diffusion_dt);

            /*
            UPDATED TO SAVE THE STATE OF THE MODEL

            TODO: REMOVE AFTER TESTING
            */
            // Save only the first 10 iterations
            // test_iteration_count++;
            // if (test_iteration_count <= 10) {
            //     // Create a dynamic filename: "output_iter_1.csv", "output_iter_2.csv", etc.
            //     std::string filename = "test_iter_" + std::to_string(test_iteration_count) + ".csv";
                
            //     // Choose the correct prefix based on which project folder you are in
            //     // dump_simulation_state(microenvironment, "original_" + filename);
            //     dump_simulation_state(microenvironment, "fpga_" + filename);
            // }

            // // Exit after 10 iterations to save time
            // if (test_iteration_count == 10) {
            //     std::cout << std::endl << ">>> Test phase completed: 10 iterations saved. Exiting now. <<<" << std::endl;
            //     exit(0); 
            // }

            
            // run PhysiCell
            evonano_cell_division((Cell_Container *) microenvironment.agent_container, PhysiCell_globals.current_time);
            evonano_CSC_deatachment((Cell_Container *) microenvironment.agent_container, PhysiCell_globals.current_time);
            ((Cell_Container *) microenvironment.agent_container)->update_all_cells(PhysiCell_globals.current_time);

            /*
              Custom add-ons could potentially go here.
            */

            // svc, 2.12.2019, microenvironment expansion
            // update 18.03.2020 (dynamic_space_expansion config)
            if (default_microenvironment_options.dynamic_space_expansion) {
                if (check_conditions_to_expand()) {
                    std::cout << "Domain expansion" << std::endl;
                    expand_microenvironment(mechanics_voxel_size);
                    num_expansions++;

                    if (expand_vascularity) {
                        double seed_multiplier = pow(num_expansions, expansion_scale_factor);
                        expand_vascular_network(seed_multiplier);
                    }

                    expand_healthy_cells(num_expansions);

                    cells_report();
                }
            }

            // check NP activation
            if (NP_simul) {
                int cc = count_cancer_cells();
                if (!NP_active && cc > cancer_cells_to_activate_NP) {
                    std::cout << std::endl << "# cancer cells = " << cc << " --- ACTIVATING NP" << std::endl << std::endl;
                    NP_active = true;
                }
            }

            PhysiCell_globals.current_time += diffusion_dt;
        }

        if( PhysiCell_settings.enable_legacy_saves == true )
        {
            log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);
            report_file.close();
        }
    }
    catch( const std::exception& e )
    { // reference to the base of a polymorphic object
        std::cout << e.what(); // information from length_error printed
    }

    // save a final simulation snapshot

    sprintf( filename , "%s/final" , PhysiCell_settings.folder.c_str() );
    save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time );

    sprintf( filename , "%s/final.svg" , PhysiCell_settings.folder.c_str() );
    SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function );


    // timer

    std::cout << std::endl << "Total simulation runtime: " << std::endl;
    BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() );
    std::cout << std::endl;

    microenvironment.release_opencl();

    return EXIT_SUCCESS;
}
