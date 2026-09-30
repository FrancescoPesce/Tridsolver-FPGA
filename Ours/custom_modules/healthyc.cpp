/**
 * EvoNano healthy cells
 * 
 * Author: Milos Savic (svc@dmi.uns.ac.rs)
 */

#include "healthyc.h"

Cell_Definition healthy_cell_def;

double HC_NP_resistance = 0;             // the resistance of HCs to NP
int NP_HC_deaths = 0;                    // the number of HC deaths caused by NP 


void init_healthy_cell_params() {
    if (NP_simul) {
        HC_NP_resistance = parameters.doubles("healthy_cell_NP_resistance");
    }
}

void healthy_cell_handle_NP(Cell* c) {
    // is NP simulated?
    if (!NP_simul)
        return;

    // is NP active?
    if (!NP_active)
        return;

    int NP_index = microenvironment.find_density_index("NP");
    double NP_near = microenvironment.nearest_density_vector(c->position)[NP_index];

    if (NP_near > HC_NP_resistance) {
        int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index("Apoptosis");
        c->start_death(apoptosis_model_index);
        ++NP_HC_deaths;
    }
}

void healthy_cell_update_phenotype(Cell* c, Phenotype& phenotype, double dt) {
    // don't bother with dying cells (cells for which a death cycle started)
    if (c->phenotype.death.dead)
        return;

    update_cell_and_death_parameters_O2_based(c, phenotype, dt);
    healthy_cell_handle_NP(c);
}


void create_healthy_cell_definition() {
    healthy_cell_def = cell_defaults; 
	healthy_cell_def.type = HEALTHY_CELL_TYPE; 
	healthy_cell_def.name = HEALTHY_CELL_NAME;
    healthy_cell_def.parameters.pReference_live_phenotype = &(healthy_cell_def.phenotype);
    
    // turn off motility
    healthy_cell_def.phenotype.motility.is_motile = false;

    // turn off cell division
    int cycle_start_index = live.find_phase_index(PhysiCell_constants::live); 
	int cycle_end_index = live.find_phase_index(PhysiCell_constants::live);
    healthy_cell_def.phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = 0.0;  

    // turn off volume updates
    healthy_cell_def.functions.volume_update_function = NULL;

    // turn off cell death
    int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index("Apoptosis");
    healthy_cell_def.phenotype.death.rates[apoptosis_model_index] = 0.0;
    int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index("Necrosis");
    healthy_cell_def.phenotype.death.rates[necrosis_model_index] = 0.0;

    healthy_cell_def.functions.update_phenotype = healthy_cell_update_phenotype;
}


void init_healthy_cells() {
    init_healthy_cell_params();
    create_healthy_cell_definition();
}


Cell* create_healthy_cell(double x, double y, double z) {
    Cell* c = create_cell(healthy_cell_def);
	c->assign_position(x, y, z);
    c->is_movable = false;
    return c;
}


void expand_healthy_cells(int num_expansions) {
    int rhc_index = parameters.ints.find_index("random_healthy_cells");
    int ef_index = parameters.ints.find_index("random_healthy_cells_expansion_factor");
    if (rhc_index == -1)
        return;

    int numhc = parameters.ints(rhc_index);

    double expansion_factor = 0;
    if (ef_index != -1)
        expansion_factor = parameters.doubles(ef_index);

    double multiplier = pow(num_expansions, expansion_factor);

	std::vector<double> borders = space_extension_for_random_sampling();

    int total = numhc * multiplier;
    std::cout << "Expanding healthy cells, creating " << total << " healthy cells" << std::endl;

    for (int i = 0; i < total; i++) {
        double x = uniformABCD(borders[0], borders[1], borders[2], borders[3]);
		double y = uniformABCD(borders[4], borders[5], borders[6], borders[7]);
		double z = uniformABCD(borders[8], borders[9], borders[10], borders[11]);
		//std::cout << "EXP healthy cell " << x << ", " << y << ", " << z << std::endl;
		create_healthy_cell(x, y, z); 
    }
}