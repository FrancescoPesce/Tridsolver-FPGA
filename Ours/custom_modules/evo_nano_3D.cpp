/**
 * EvoNano main module
 * 
 * Author: Milos Savic (svc@dmi.uns.ac.rs) & Vladimir Kurbalija (Phase 2)
 */

#include "evo_nano_3D.h"

// CAF dynamics
int vascular_expansions = 0;
int caf_cells = 0;
double CAF_vascular_ratio;

// CSC deatachment
bool CSC_detachment_active = false;
double CSC_detachment_prob = 0;
double CSC_detachment_dist = std::numeric_limits<double>::max();
int detached_CSC = 0;

void create_cell_types(void)
{
	SeedRandom(parameters.ints("random_seed")); // or specify a seed here

	// housekeeping

	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );

	// Name the default cell type
	cell_defaults.type = 0;
	cell_defaults.name = "evo-nano cell";

	// set default cell cycle model
    live.transition_rate(0,0) = 0.00432;            // default rate is 0.0432 / 60.0
	cell_defaults.functions.cycle_model = live;

	// set default_cell_functions;
	cell_defaults.functions.update_phenotype = update_cell_and_death_parameters_O2_based;

	// make sure the defaults are self-consistent.

	cell_defaults.phenotype.secretion.sync_to_microenvironment(&microenvironment);
	cell_defaults.phenotype.sync_to_functions(cell_defaults.functions);

	// first find index for a few key variables.
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index("Apoptosis");
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index("Necrosis");
	int oxygen_substrate_index = microenvironment.find_density_index("oxygen");

	// initially no necrosis
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0;

	// set oxygen uptake / secretion parameters for the default cell type
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 10;
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 38;
	
    init_params_and_cell_definitions();
    build_cell_definitions_maps(); 
	//display_cell_definitions( std::cout );
}


void init_params_and_cell_definitions() {
	// phase 1
	NP_simul = parameters.bools("NP_active");
    if (NP_simul) {
        int NP_index = microenvironment.find_density_index("NP");
        NP_simul = NP_index != -1;
        if (!NP_simul) {
            std::cout << "[WARNING], NP active but the substrate description is missing, deactivating..." << std::endl;
		} else {
			std::cout << "NP simulated..." << std::endl;
		}
	}
	
	// phase 2 -- prostaglandin
	prostaglandin_simul = parameters.bools("prostaglandin_active");
    if (prostaglandin_simul) {
        int prostaglandin_index = microenvironment.find_density_index("prostaglandin");
        prostaglandin_simul = prostaglandin_index != -1;
        if (!prostaglandin_simul) {
            std::cout << "[WARNING], prostaglandin active but the substrate description is missing, deactivating..." << std::endl;
		} else {
			std::cout << "prostaglandin simulated..." << std::endl;
		}
	}

	// phase 3	
	cytokine_caf_simul = parameters.bools("CAF_active");
    if (cytokine_caf_simul){
        int cytokine_index = microenvironment.find_density_index("cytokine");
        cytokine_caf_simul = cytokine_index != -1;
        if (!cytokine_caf_simul) {
            std::cout << "[WARNING], CYTOKINE & CAF active but the substrate description is missing, deactivating..." << std::endl;
        } else {
			std::cout << "CYTOKINE & CAF simulated..." << std::endl;
			CAF_vascular_ratio = parameters.doubles("CAF_vascular_ratio");
		}
    }

	init_cancer_cells();
	init_vascular_cells();
	init_healthy_cells();


    CSC_detachment_prob = parameters.doubles("CSC_detachment_probability");
    CSC_detachment_active = CSC_detachment_prob > 0;
    if (CSC_detachment_active) {
        CSC_detachment_dist = parameters.doubles("CSC_detachment_critical_proximity");
		std::cout << "CSC detachment simulated..." << std::endl;
	}
}


void setup_microenvironment(void)
{
	// make sure to override and go back to 2D
	if( default_microenvironment_options.simulate_2D == true )
	{
		std::cout << "Warning: overriding XML config option and setting to 3D!" << std::endl;
		default_microenvironment_options.simulate_2D = false;
	}

	// initialize BioFVM
	initialize_microenvironment();

	return;
}


void create_randomly_placed_cells() {
	std::vector<double> borders = space_borders_for_random_sampling();

	int rhc_index = parameters.ints.find_index("random_healthy_cells");
	if (rhc_index != -1) {
		int numhc = parameters.ints(rhc_index);
		std::cout << "Creating " << numhc << " randomly placed healthy cells" << std::endl;
		
		for (int i = 0; i < numhc; i++) {
			double x = uniformAB(borders[0], borders[1]);
			double y = uniformAB(borders[2], borders[3]);
			double z = uniformAB(borders[4], borders[5]);
			//std::cout << "Healthy cell " << x << ", " << y << ", " << z << std::endl;

			create_healthy_cell(x, y, z);
		}
	}

	
	int rcc_index = parameters.ints.find_index("random_cancer_cells");
	if (rcc_index != -1) {
		int numcc = parameters.ints(rcc_index);
		std::cout << "Creating " << numcc << " randomly placed cancer cells" << std::endl;
		
		for (int i = 0; i < numcc; i++) {
			double x = uniformAB(borders[0], borders[1]);
			double y = uniformAB(borders[2], borders[3]);
			double z = uniformAB(borders[4], borders[5]);
			//std::cout << "Cancer cell " << x << ", " << y << ", " << z << std::endl;

			create_cancer_cell(x, y, z, false);
		}
	}

	// and our favorite cancer cell in the center
	create_cancer_cell(0, 0, 0, false);
}


void setup_tissue(void)
{
	init_vascular_network();
	create_randomly_placed_cells();
}


void cells_report() {
	int num_active = 0, num_movable = 0, num_outdomain = 0, num_dead = 0;
	int hist[6] = {0, 0, 0, 0, 0, 0};

	for (int i = 0; i < all_cells->size(); i++) {
        Cell* c = (*all_cells)[i];
		if (c->is_active)            num_active++;
		if (c->is_movable)           num_movable++;
		if (c->is_out_of_domain)     num_outdomain++;
		if (c->phenotype.death.dead) num_dead++;

		int type = c->type;
		if (type > 5) {
			std::cout << "[Warning!!!] Unknown cell type " << type << std::endl;
		} else {
			hist[type]++;
		}
	}

	std::cout << "Active = " << num_active << ", Movable = " << num_movable 
	        << ", Out of domain = " << num_outdomain << ", Dead = " << num_dead << std::endl;

	if (hist[0] > 0)
		std::cout << "Warning!!! There are cells with the default cell type, #total = " << hist[0] << std::endl;
	
	// histogram
	std::cout << CANCER_CELL_NAME << " = " << hist[CANCER_CELL_TYPE] << ", "
		      << CANCER_STEM_CELL_NAME << " = " << hist[CANCER_STEM_CELL_TYPE] << ", "
			  << CAF_CELL_NAME << " = " << hist[CAF_CELL_TYPE] << ", " 
			  << VASCULAR_CELL_NAME << " = " << hist[VASCULAR_CELL_TYPE] << ", "
			  << HEALTHY_CELL_NAME << " = " << hist[HEALTHY_CELL_TYPE] << std::endl;

	if (CSC_detachment_active) {
		std::cout << "Detached CSC = " << detached_CSC << std::endl; 
	}
}


std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with flow cytometry coloring
	std::vector<std::string> output = false_cell_coloring_live_dead(pCell);

	// if the cell is motile and not dead, paint it black

	if (!pCell->phenotype.death.dead) {
		switch (pCell->type) {
			case CANCER_CELL_TYPE:
				output[0] = "grey";
        		output[2] = "grey";
				break;
			case CANCER_STEM_CELL_TYPE:
				output[0] = "blue";
        		output[2] = "blue";
				break;
			case VASCULAR_CELL_TYPE:
				output[0] = "maroon";
        		output[2] = "maroon";
				break;
			case HEALTHY_CELL_TYPE:
				output[0] = "green";
        		output[2] = "green";
				break;
			case CAF_CELL_TYPE:
				output[0] = "yellow";
				output[2] = "yellow";
				break;
			default:
				output[0] = "black";
        		output[2] = "black";
		}
	}
	
	return output;
}


void create_CAF_near_vascular(Cell* vascular_cell) {
	double ratio = caf_cells / vascular_expansions;
	if (ratio < CAF_vascular_ratio) {
		double vx = vascular_cell->position[0];
		double vy = vascular_cell->position[1];
		double vz = vascular_cell->position[2];

		std::vector<double> rndvec = UniformOnUnitSphere();
		rndvec *= 1.5 * cell_defaults.phenotype.geometry.radius;

		double x = vx + rndvec[0];
		double y = vy + rndvec[1];
		double z = vz + rndvec[2];

		create_CAF_cell(x, y, z);

		++caf_cells;
	}
}

// custom evo-nano cell division rules
void evonano_cell_division(Cell_Container* cell_container, double t) {
	double phenotype_dt_tolerance = 0.001 * phenotype_dt;
	double time_since_last_cycle = t - cell_container->last_cell_cycle_time;
	if (fabs(time_since_last_cycle - phenotype_dt) >= phenotype_dt_tolerance)
		return;

	for (int i = 0; i < cell_container->cells_ready_to_divide.size(); i++){
		Cell* c = cell_container->cells_ready_to_divide[i];
		int type = c->type;

		if (type == CANCER_CELL_TYPE) {
			divide_cancer_cell(c);
		} 
		else if (type == CANCER_STEM_CELL_TYPE) {
			divide_cancer_stem_cell(c);
		} 
		else if (type == VASCULAR_CELL_TYPE) {
			divide_vascular_cell(c);
			++vascular_expansions;

			if (cytokine_caf_simul) {
				create_CAF_near_vascular(c);
			}
		}
		else {
			c->divide();
		}
	}

	cell_container->cells_ready_to_divide.clear();
	//cells_report();
}


bool critical_distance_to_vascular(Cell* csc_cell) {
	std::vector<double> csc_pos = csc_cell->position;
	
	for (int i = 0; i < vascular_cells->size(); i++) {
		Cell* vc = (*vascular_cells)[i];
		std::vector<double> vc_pos = vc->position;
		double d = dist(csc_pos, vc_pos);
		if (d < CSC_detachment_dist)
			return true;
	}

	return false;
}

// CSC deatachment rules
void evonano_CSC_deatachment(Cell_Container* cell_container, double t) {
	if (!CSC_detachment_active)
		return;

	double phenotype_dt_tolerance = 0.001 * phenotype_dt;
	double time_since_last_cycle = t - cell_container->last_cell_cycle_time;
	if (fabs(time_since_last_cycle - phenotype_dt) >= phenotype_dt_tolerance)
		return;

	std::vector<Cell*> detached_cells;

	for (int i = 0; i < all_cells->size(); i++) {
        Cell* c = (*all_cells)[i];
		int type = c->type;

		if (type == CANCER_STEM_CELL_TYPE) {
			bool detach = UniformRandom() < CSC_detachment_prob;
			if (detach && critical_distance_to_vascular(c)) {
				detached_CSC++;
				detached_cells.push_back(c);
			}
		}
	}

	for (int i = 0; i < detached_cells.size(); i++) {
		Cell* c = detached_cells[i];
		c->die();
	}
}




