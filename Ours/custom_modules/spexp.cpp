/**
 * Simulation space expansion
 * 
 * Author: Milos Savic (svc@dmi.uns.ac.rs)
 */

#include "spexp.h"

bool critical_proximity(double x, double y) {
    double d = x - y;
    if (d < 0)
        d *= -1.0;

    return d <= default_microenvironment_options.critical_proximity;
}

double dist(double x, double y) {
    double d = x - y;
    if (d < 0)
        d *= -1.0;
    
    return d;
}

bool check_conditions_to_expand() {
    double xmin = default_microenvironment_options.X_range[0];
    double xmax = default_microenvironment_options.X_range[1];
    double ymin = default_microenvironment_options.Y_range[0];
    double ymax = default_microenvironment_options.Y_range[1];
    double zmin = default_microenvironment_options.Z_range[0];
    double zmax = default_microenvironment_options.Z_range[1];

    double x_size = xmax - xmin;
    double y_size = ymax - ymin;
    double z_size = zmax - zmin;
    double max_size = default_microenvironment_options.max_size;
    if (x_size >= max_size && y_size >= max_size && z_size >= max_size) {
        return false;
    }
    
    for (int i = 0; i < all_cells->size(); i++) {
        double x = (*all_cells)[i]->position[0];
        double y = (*all_cells)[i]->position[1];
        double z = (*all_cells)[i]->position[2];
        std::string type = (*all_cells)[i]->type_name;

        bool critical = 
            critical_proximity(x, xmin) || critical_proximity(x, xmax) ||
            critical_proximity(y, ymin) || critical_proximity(y, ymax) ||
            critical_proximity(z, zmin) || critical_proximity(z, zmax);

        if (critical) {
            std::cout << "A cell with critical proximity detected [" << x << ", " << y << ", " << z << "], type = " << type << std::endl;
            return true;
        }
    }

    return false;
}

typedef struct {
    std::vector<double> cords;
    std::vector<double> density;
} space_density_state;

std::vector<space_density_state> density_state;

void save_current_densities() {
    density_state.clear();

    for (unsigned int i = 0; i < microenvironment.mesh.x_coordinates.size(); i++) {
        for (unsigned int j = 0; j < microenvironment.mesh.y_coordinates.size(); j++) {
            for (unsigned int k = 0; k < microenvironment.mesh.z_coordinates.size(); k++) {
                double x = microenvironment.mesh.x_coordinates[i];
                double y = microenvironment.mesh.y_coordinates[j];
                double z = microenvironment.mesh.z_coordinates[k];
                
                space_density_state ss;
                ss.cords.clear();
                ss.cords.push_back(x); 
                ss.cords.push_back(y); 
                ss.cords.push_back(z); 
                std::vector<double> den = microenvironment(i, j, k);
                ss.density.clear();
                for (int d = 0; d < den.size(); d++)
                    ss.density.push_back(den[d]);
                density_state.push_back(ss); 
            }
        } 
    }
}

void restore_previous_densities() {
    for (int i = 0; i < density_state.size(); i++) {
        std::vector<double> cords = density_state[i].cords;
        std::vector<double> dens = density_state[i].density;
        microenvironment.update_density_vector(cords, dens);
    }
}

void expand_microenvironment(double mechanics_voxel_size) {
    int expansion_delta = default_microenvironment_options.expansion_delta;
    std::cout << "Expanding microenvironment, current time " << PhysiCell_globals.current_time << std::endl;
    save_current_densities();

    /*
    std::cout << "Before expandning " << std::endl;
    for (unsigned int i = 0; i < microenvironment.mesh.x_coordinates.size(); i++) {
        for (unsigned int j = 0; j < microenvironment.mesh.y_coordinates.size(); j++) {
            for (unsigned int k = 0; k < microenvironment.mesh.z_coordinates.size(); k++) {
                std::vector<double> den = microenvironment(i, j, k);
                std::cout << den[0] << " ";
            }
        }
    }
    std::cout << std::endl;
    */

    default_microenvironment_options.X_range[0] -= expansion_delta;
    default_microenvironment_options.X_range[1] += expansion_delta;
    default_microenvironment_options.Y_range[0] -= expansion_delta;
    default_microenvironment_options.Y_range[1] += expansion_delta;  
    default_microenvironment_options.Z_range[0] -= expansion_delta;
    default_microenvironment_options.Z_range[1] += expansion_delta;  
            
    setup_microenvironment();
    restore_previous_densities();

    /*
    std::cout << "After restoring " << std::endl;
    for (unsigned int i = 0; i < microenvironment.mesh.x_coordinates.size(); i++) {
        for (unsigned int j = 0; j < microenvironment.mesh.y_coordinates.size(); j++) {
            for (unsigned int k = 0; k < microenvironment.mesh.z_coordinates.size(); k++) {
                std::vector<double> den = microenvironment(i, j, k);
                std::cout << den[0] << " ";
            }
        }
    }
    std::cout << std::endl;
    */

    Cell_Container* new_cell_container = create_cell_container_for_microenvironment(microenvironment, mechanics_voxel_size);
    
    for (int i = 0; i < all_cells->size(); i++) {
        (*all_cells)[i]->reset_container();
    }

    for (int i = 0; i < all_cells->size(); i++) {
        (*all_cells)[i]->register_microenvironment(&microenvironment);
        (*all_cells)[i]->assign_position((*all_cells)[i]->position[0], (*all_cells)[i]->position[1], (*all_cells)[i]->position[2]);
        (*all_cells)[i]->phenotype.sync_to_microenvironment(&microenvironment);
    }

    microenvironment.compute_all_gradient_vectors();
    microenvironment.diffusion_solver_setup_done = false;
    
    std::cout << "Microenviroment expansion finished... " << std::endl;
}
