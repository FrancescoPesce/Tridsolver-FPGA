/**
 * Evo-Nano cancer cells
 * 
 * Author: Milos Savic (svc@dmi.uns.ac.rs) & Vladimir Kurbalija (Phase 2)
 */

#include "cancc.h"

Cell_Definition cancer_cell_def;
Cell_Definition cancer_stem_cell_def;
Cell_Definition caf_cell_def;

bool cancer_stem_cells_enabled = true;      

double CSCP = 0;             // the probability that CC will be divided into 1 CC + 1 CSC 
double SDP = 0;              // the probability that CSC division will be symmetric                  
double DOUBLE_CC_PROB = 0;   // the probability that CSC will be divided into two CC

// angiogen properties
double ANG_IN_CANCER = 0;    // density of angiogen in cancer cells

// phase 1 -- NP properties
double CC_NP_LOWER = 0;      // the lower bound of NP threshold for cancer cells to activate apoptosis 
double CC_NP_UPPER = 0;      // the upper bound of NP threshold for cancer cells to activate apoptosis 
int NP_CC_deaths = 0;        // the number of CC deaths caused by NP 

// phase 2 -- NP for stem cells (Make DCC more susceptible to apoptosis than CSC)
double CSC_NP_LOWER = 0;      // the lower bound of NP threshold for cancer stem cells to activate apoptosis 
double CSC_NP_UPPER = 0;      // the upper bound of NP threshold for cancer stem cells to activate apoptosis 

// phase 2 -- prostaglandin properties
bool prostaglandin_simul;          			// is prostaglandin simulated 
double PROSTAGLANDIN_IN_CANCER = 0;       	// density of prostaglandin in CC & CSC cells
double cc_prostaglandin_thr = 0;       		// prostaglandin threshold for CC -> CSC
int prostaglandin_conversions = 0;			// counter

// phase 3 -- cytokine properties
bool cytokine_caf_simul;          // are cytokine & CAF simulated 
double CYTOKINE_IN_CAF = 0;       // density of cytokine in CAF cells
double cc_cytokine_thr = 0;       // cytokine threshold for CC -> CSC
int cytokine_conversions = 0;


double EPS = 0.0001;         // a small number close to zero (< EPS used instead of == 0)

void init_cancer_cell_params() {
    cancer_stem_cells_enabled = parameters.bools("cancer_stem_cells_enabled");
	CSCP = parameters.doubles("CSC_probability_CC_division");
	SDP = parameters.doubles("symmetric_division_probability");
	DOUBLE_CC_PROB = parameters.doubles("symmetric_division_double_CC_probability");

    int angiogen_index = microenvironment.find_density_index("angiogen");
	bool angiogen_simul = angiogen_index != -1;
    if (angiogen_simul) {
        ANG_IN_CANCER = parameters.doubles("angiogen_in_cancer_cell");
        std::cout << "ANGIOGEN simulated" << std::endl;
    }

    // phase 1 & 2
    if (NP_simul) {
        CC_NP_LOWER = parameters.doubles("cancer_cell_NP_resistance_lower_bound");
        CC_NP_UPPER = parameters.doubles("cancer_cell_NP_resistance_upper_bound");
        CSC_NP_LOWER = parameters.doubles("cancer_stem_cell_NP_resistance_lower_bound");
        CSC_NP_UPPER = parameters.doubles("cancer_stem_cell_NP_resistance_upper_bound");
    }
    
    // phase 2
    if (prostaglandin_simul) {
        PROSTAGLANDIN_IN_CANCER = parameters.doubles("prostaglandin_in_cancer_cell");
        cc_prostaglandin_thr = parameters.doubles("cancer_cell_prostaglandin_thr");
    }
    
    // phase 3
    if (cytokine_caf_simul) {
        CYTOKINE_IN_CAF = parameters.doubles("cytokine_in_CAF");
        cc_cytokine_thr = parameters.doubles("cancer_cell_cytokine_resistance");
    }
}


/*
void setup_angiogen_in_cancer_cell(Cell* c) {
	int angiogen_index = microenvironment.find_density_index("angiogen");
	if (angiogen_index != -1) {
		microenvironment.nearest_density_vector(c->position)[angiogen_index] = ANG_IN_CANCER;
	}
}
*/


void cancer_cell_handle_NP(Cell* c) {
    // is NP simulated?
    if (!NP_simul)
        return;

    // is NP active?
    if (!NP_active)
        return;

    int ind = c->custom_data.find_variable_index("NP_resistance");
    double NP_resistance = c->custom_data[ind];
    int NP_index = microenvironment.find_density_index("NP");
    double NP_near = microenvironment.nearest_density_vector(c->position)[NP_index];

    int cycle_start_index = live.find_phase_index(PhysiCell_constants::live); 
	int cycle_end_index = live.find_phase_index(PhysiCell_constants::live);
    double original_division_rate = parameters.doubles("cancer_cell_division_rate");
    if (c->type == CANCER_STEM_CELL_TYPE)
        original_division_rate = parameters.doubles("cancer_stem_cell_division_rate");

    int NP_prev_index = c->custom_data.find_variable_index("NP_prev");

    if (NP_near > NP_resistance) {
        int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index("Apoptosis");
        c->start_death(apoptosis_model_index);
        NP_CC_deaths++;
    } else if (NP_near < EPS) {
        c->phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = original_division_rate; 
    } else {
        double current_transition_rate = c->phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index);
        double NP_prev = c->custom_data[NP_prev_index];
        if (NP_prev > EPS) {
            double NP_increase = NP_near / NP_prev; 
            double adjusted_division_rate = current_transition_rate / NP_increase;
            c->phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = adjusted_division_rate;
        }
    }

    c->custom_data[NP_prev_index] = NP_near;
}


void cancer_cell_handle_cytokine(Cell* c) {
    int cytokine_index = microenvironment.find_density_index("cytokine");
    double cytokine_near = microenvironment.nearest_density_vector(c->position)[cytokine_index];
    if (cytokine_near > cc_cytokine_thr) {
        c->convert_to_cell_definition(cancer_stem_cell_def);
        cytokine_conversions++;
    }
}

void cancer_cell_handle_prostaglandin(Cell* c) {
    int prostaglandin_index = microenvironment.find_density_index("prostaglandin");
    double prostaglandin_near = microenvironment.nearest_density_vector(c->position)[prostaglandin_index];
    if (prostaglandin_near > cc_prostaglandin_thr) {
        c->convert_to_cell_definition(cancer_stem_cell_def);
        prostaglandin_conversions++;
    }
}

void cancer_cell_phenotype_update(Cell* c, Phenotype& phenotype, double dt) {
    // diffuse prostaglandin in dead cell
    if (c->phenotype.death.dead) {
    	if (prostaglandin_simul) {
			int prostaglandin_index = microenvironment.find_density_index("prostaglandin");
			microenvironment.nearest_density_vector(c->position)[prostaglandin_index] = PROSTAGLANDIN_IN_CANCER;
		}
	
        return;
    }

    update_cell_and_death_parameters_O2_based(c, phenotype, dt);
    //cancer_cell_handle_angiogen(c);
    
	cancer_cell_handle_NP(c);   
    
    if (prostaglandin_simul) {
    	cancer_cell_handle_prostaglandin(c);
    }
    
    if (cytokine_caf_simul) {
    	cancer_cell_handle_cytokine(c);
    }
}


void create_cancer_cell_definition() {
    cancer_cell_def = cell_defaults; 
	cancer_cell_def.type = CANCER_CELL_TYPE; 
	cancer_cell_def.name = CANCER_CELL_NAME;
    cancer_cell_def.parameters.pReference_live_phenotype = &(cancer_cell_def.phenotype); 

    // cancer cell motility
    cancer_cell_def.phenotype.motility.persistence_time = parameters.doubles("cancer_cell_persistence_time"); 
    cancer_cell_def.phenotype.motility.migration_speed = parameters.doubles("cancer_cell_migration_speed");  
    
    // cancer cell adhesion
    cancer_cell_def.phenotype.mechanics.cell_cell_adhesion_strength *=
        parameters.doubles("cancer_cell_relative_adhesion"); 

    // cancer cell apoptosis rate
    int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index("Apoptosis");
    cancer_cell_def.phenotype.death.rates[apoptosis_model_index] =
            parameters.doubles("cancer_cell_apoptosis_rate"); 

    // set cancer cell division rate
    int cycle_start_index = live.find_phase_index(PhysiCell_constants::live); 
	int cycle_end_index = live.find_phase_index(PhysiCell_constants::live);
    cancer_cell_def.phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = 
        parameters.doubles("cancer_cell_division_rate");

    
    cancer_cell_def.functions.update_phenotype = cancer_cell_phenotype_update;

    cancer_cell_def.custom_data.add_variable("NP_resistance", "dimensionless", 0);
    cancer_cell_def.custom_data.add_variable("NP_prev", "dimensionless", 0);

    // set angiogen uptake and secretion parameters
    int ang_index = microenvironment.find_density_index("angiogen");
    if (ang_index != -1) {
        cancer_cell_def.phenotype.secretion.uptake_rates[ang_index] = 0;
	    cancer_cell_def.phenotype.secretion.secretion_rates[ang_index] = ANG_IN_CANCER;
	    cancer_cell_def.phenotype.secretion.saturation_densities[ang_index] = ANG_IN_CANCER;
    }
}


void create_cancer_stem_cell_definition() {
    cancer_stem_cell_def = cell_defaults; 
	cancer_stem_cell_def.type = CANCER_STEM_CELL_TYPE; 
	cancer_stem_cell_def.name = CANCER_STEM_CELL_NAME;
    cancer_stem_cell_def.parameters.pReference_live_phenotype = &(cancer_stem_cell_def.phenotype);  

    // cancer stem cell necrosis params
    cancer_stem_cell_def.parameters.o2_necrosis_threshold = parameters.doubles("cancer_stem_cell_o2_necrosis_threshold");
    cancer_stem_cell_def.parameters.o2_necrosis_max = parameters.doubles("cancer_stem_cell_o2_necrosis_max");
    
    // cancer stem cell motility
    cancer_stem_cell_def.phenotype.motility.persistence_time = parameters.doubles("cancer_stem_cell_persistence_time"); 
    cancer_stem_cell_def.phenotype.motility.migration_speed = parameters.doubles("cancer_stem_cell_migration_speed"); 

    // cancer stem cell adhesion
    cancer_stem_cell_def.phenotype.mechanics.cell_cell_adhesion_strength *=
            parameters.doubles("cancer_stem_cell_relative_adhesion"); 

    // cancer stem cell apoptosis rate
    int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index("Apoptosis");
    cancer_stem_cell_def.phenotype.death.rates[apoptosis_model_index] =
            parameters.doubles("cancer_stem_cell_apoptosis_rate"); 

    // cancer stem cell division rate
    int cycle_start_index = live.find_phase_index(PhysiCell_constants::live); 
	int cycle_end_index = live.find_phase_index(PhysiCell_constants::live);
    cancer_stem_cell_def.phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = 
        parameters.doubles("cancer_stem_cell_division_rate");

    cancer_stem_cell_def.functions.update_phenotype = cancer_cell_phenotype_update;

    cancer_stem_cell_def.custom_data.add_variable("NP_resistance", "dimensionless", 0);
    cancer_stem_cell_def.custom_data.add_variable("NP_prev", "dimensionless", 0);
}


void create_caf_cell_definition() {
    caf_cell_def = cell_defaults; 
	caf_cell_def.type = CAF_CELL_TYPE; 
	caf_cell_def.name = CAF_CELL_NAME;
    caf_cell_def.parameters.pReference_live_phenotype = &(caf_cell_def.phenotype);

    // turn off cell division
    int cycle_start_index = live.find_phase_index(PhysiCell_constants::live); 
	int cycle_end_index = live.find_phase_index(PhysiCell_constants::live);
    caf_cell_def.phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = 0.0; 

    // turn off cell death
    int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index("Apoptosis");
    caf_cell_def.phenotype.death.rates[apoptosis_model_index] = 0.0;
    int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index("Necrosis");
    caf_cell_def.phenotype.death.rates[necrosis_model_index] = 0.0;

    // turn off volume updates
    caf_cell_def.functions.volume_update_function = NULL;

    // turn off phenotype updates
    caf_cell_def.functions.update_phenotype = NULL;

    // set cytokine uptake and secretion parameters
    int cytokine_index = microenvironment.find_density_index("cytokine");
    if (cytokine_index != -1) {
        caf_cell_def.phenotype.secretion.uptake_rates[cytokine_index] = 0;
        caf_cell_def.phenotype.secretion.secretion_rates[cytokine_index] = CYTOKINE_IN_CAF;
        caf_cell_def.phenotype.secretion.saturation_densities[cytokine_index] = CYTOKINE_IN_CAF;
    }
}

void init_cancer_cells() {
    init_cancer_cell_params();
    create_cancer_cell_definition();
    create_cancer_stem_cell_definition();
    create_caf_cell_definition();
}


Cell* create_cancer_cell(double x, double y, double z, bool stem) {
	Cell* c;
	
	if (stem)
		c = create_cell(cancer_stem_cell_def);
	else
		c = create_cell(cancer_cell_def);

	c->assign_position(x, y, z);
	//setup_angiogen_in_cancer_cell(c);

    // phase 1 & 2: creating NP resistance levels for CC and CSC
	if (NP_simul) {
        int ind = c->custom_data.find_variable_index("NP_resistance");
        if (stem)
        	c->custom_data[ind] = CSC_NP_LOWER + UniformRandom() * (CSC_NP_UPPER - CSC_NP_LOWER);
        else
        	c->custom_data[ind] = CC_NP_LOWER + UniformRandom() * (CC_NP_UPPER - CC_NP_LOWER);
        //std::cout << "Cancer cell created, NP resistance = " << c->custom_data[ind] << std::endl;
    }

    return c;
}

Cell* create_CAF_cell(double x, double y, double z) {
    Cell* c = create_cell(caf_cell_def);
    c->assign_position(x, y, z);
    return c;
}


void divide_cancer_cell(Cell* c) {
	if (!cancer_stem_cells_enabled) {
		Cell* new_cell = c->divide();
		//setup_angiogen_in_cancer_cell(new_cell);
		return;
	}

	double pCSC = UniformRandom();    // cancer stem cell probability
	Cell* new_cell = c->divide();
	if (pCSC <= CSCP) {
		new_cell->convert_to_cell_definition(cancer_stem_cell_def);	
	} 

	//setup_angiogen_in_cancer_cell(new_cell);
}


void divide_cancer_stem_cell(Cell* c) {
	if (!cancer_stem_cells_enabled) {
		std::cout << "Warning!!! Cancer stem cell division, but cancer stem cells disabled!" << std::endl;
		return;
	}

	double symetricDivision = UniformRandom();
	Cell* new_cell = c->divide();
	if (symetricDivision <= SDP) {
		// symetric division
		double pCC = UniformRandom();
		if (pCC <= DOUBLE_CC_PROB) {
			c->convert_to_cell_definition(cancer_cell_def);
			new_cell->convert_to_cell_definition(cancer_cell_def);
		}
	} else {
		// asymetric division
		new_cell->convert_to_cell_definition(cancer_cell_def);	
	}

	//setup_angiogen_in_cancer_cell(new_cell);
}

