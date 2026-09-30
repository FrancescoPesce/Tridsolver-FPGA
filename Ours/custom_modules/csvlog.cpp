/**
 * CSV cell log
 * 
 * Author: Milos Savic (svc@dmi.uns.ac.rs)
 */


#include "csvlog.h"

CSV_log::CSV_log(std::string file_name) {
    this->file_name = file_name;

    std::ofstream output_file;
    output_file.open(file_name, std::fstream::out);

    output_file << " \n";
    output_file << "Simulation parameters:\n";
    output_file << "Space: ," << (default_microenvironment_options.simulate_2D ? "2D" : "3D") << "\n";
    output_file << "x: ," << default_microenvironment_options.X_range[1] << "\n";
    output_file << "y: ," << default_microenvironment_options.Y_range[1] << "\n";
    output_file << "z: ," << default_microenvironment_options.Z_range[1] << "\n";
    output_file << "Simulation length: ," << PhysiCell_settings.max_time << "\n";
    output_file << "Starting number of cells: ," << all_cells->size() << "\n";

    output_file.close();

    for (int i = 0; i <= 5; i++) 
        hist_dead[i] = 0;
}


void CSV_log::log() {
    std::ofstream output_file;
    output_file.open(file_name, std::fstream::app);

    int hist_live[6] = {0, 0, 0, 0, 0, 0};

    for (int i = 0; i < all_cells->size(); i++) {
        Cell* c = (*all_cells)[i];
		int type = c->type;

        if (type > 5) {
            std::cout << "[csvlog.cpp] Warning!!! Unknown cell type " << type << std::endl;
            continue;
        }

        if (c->phenotype.death.dead)
            hist_dead[type]++;
        else
            hist_live[type]++;
	}
    
    output_file << " \n";
    output_file << "Simulation:\n";
    output_file << "Time: " << PhysiCell_globals.current_time << "\n";
    output_file << "Total amount of cells: " << all_cells->size() << "\n";
    output_file << "Number of CSC: "         << hist_live[CANCER_STEM_CELL_TYPE] << "\n";
    output_file << "Number of CC: "          << hist_live[CANCER_CELL_TYPE] << "\n";
    
    //output_file << "Number of VS: "        << hist_live[VASCULAR_CELL_TYPE] << "\n"
    //output_file << "Number of HC: "        << hist_live[HEALTHY_CELL_TYPE] << "\n";
    
    output_file << "Number of dead CSC: "    <<  hist_dead[CANCER_STEM_CELL_TYPE] << "\n";
    output_file << "Number of dead CC: "     <<  hist_dead[CANCER_CELL_TYPE] << "\n";
    output_file << "Cell, type, volume, cycle model, cycle phase, x, y, z, OXYGEN\n";
    
    std::string row = "";
    int row_counter = 0;
    for (int i = 0; i < all_cells->size(); i++) {
        //int cell_type = (*all_cells)[i]->type;
        row_counter++;
        row = "";
        row += std::to_string(row_counter);
        row += ",";
        row += (*all_cells)[i]->type_name;
        row += ",";
        row += std::to_string((*all_cells)[i]->phenotype.volume.total);
        row += ",";
        row += (*all_cells)[i]->phenotype.cycle.pCycle_Model->name;
        row += ",";
        row += (*all_cells)[i]->phenotype.cycle.current_phase().name;
        row += ",";
        for (int j = 0; j < (*all_cells)[i]->position.size(); j++) {
            row += std::to_string((*all_cells)[i]->position[j]);
            if (j <= (*all_cells)[i]->position.size() - 2) {
                row += ",";
            }
        }

        // svc, oxygen logging, 12.10.2019
        row += ",";
        int o2_index = microenvironment.find_density_index("oxygen");
        if (o2_index == -1) {
            std::cout << "WARNING [microenvironment]: invalid oxygen density index (-1)";
            break;
        } 
        double o2 = (*all_cells)[i]->nearest_density_vector()[o2_index];
        row += std::to_string(o2);

        row += "\n";
        output_file << row << "\n";
    }

    output_file.close();
}