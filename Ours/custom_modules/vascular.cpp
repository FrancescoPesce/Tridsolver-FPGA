/**
 * Evo-Nano vascularity
 * 
 * Author: Milos Savic (svc@dmi.uns.ac.rs)
 */


#include "vascular.h"

Cell_Definition vascular_cell_def;

bool NP_simul = false;                    // true if NP is simulated
bool NP_active = false;                   // true if NP is active
double NP_in_vascular_cell;               // NP density in vascular cells

double ANG_THRESHOLD_VASCULAR = 0;        // angiogen critical threshold for vascular cell division
double VCEP = 0;                          // vascular cell expansion probability

int num_seed_points = 0;                  // the number of seed vascular cells
int vascular_chain_max_length = 0;        // max length of initial vascular chains
double vascularity_branching_factor = 0;  // the probability of branching in initial vascular chains
int max_branching_depth = 0;              // maximal branching depth
double expp_init_vasc;                    // the expansion probability for initial vascular cells 

std::vector<Cell*>* vascular_cells;

void init_vascular_params() {
    ANG_THRESHOLD_VASCULAR = parameters.doubles("vascular_cell_angiogen_threshold");
    VCEP = parameters.doubles("vascular_cell_expansion_probability");

    num_seed_points = parameters.ints("vascular_seed_points");
    vascular_chain_max_length = parameters.ints("vascular_chain_max_length");
    vascularity_branching_factor = parameters.doubles("vascularity_branching_factor");
    max_branching_depth = parameters.ints("vascularity_max_branching_depth");
    expp_init_vasc = parameters.doubles("initial_vascularity_expansion_probability");

    NP_simul = parameters.bools("NP_active");
    if (NP_simul)
        NP_in_vascular_cell = parameters.doubles("NP_in_vascular_cell");

    vascular_cells = new std::vector<Cell*>;
}


void vascular_cell_handle_angiogen(Cell* c) {
    int angiogen_index = microenvironment.find_density_index("angiogen");
    if (angiogen_index == -1)
        return;

    double ang_near = microenvironment.nearest_density_vector(c->position)[angiogen_index];
    if (ang_near >= ANG_THRESHOLD_VASCULAR) {
        int ffd_index = c->custom_data.find_variable_index("flagged_for_division");
        bool ffd = c->custom_data[ffd_index] == 1;
        if (ffd)
            return;

        bool expand = UniformRandom() <= VCEP;
        if (!expand)
            return;

        int numexp_index = c->custom_data.find_variable_index("num_expansions");
        int num_expansions = c->custom_data[numexp_index];
        
        expand = num_expansions == 0;
        if (expand) {
            c->custom_data[ffd_index] = 1;
            c->flag_for_division();
        }
    }
}


void vascular_cell_handle_NP(Cell* c) {
    if (!NP_simul)
        return;

    int NP_index = microenvironment.find_density_index("NP");
    if (NP_index == -1)
        return;

    if (!NP_active)
        return;

    int voxel_index = c->get_current_voxel_index();
    microenvironment.update_dirichlet_node(voxel_index, NP_index, NP_in_vascular_cell);   
}


/* 
void vascular_cell_handle_NP(Cell* c) {
    int NP_index = microenvironment.find_density_index("NP");
    if (NP_index == -1)
        return;

    if (!NP_active)
        return;

    c->phenotype.secretion.uptake_rates[NP_index] = 0;
	c->phenotype.secretion.secretion_rates[NP_index] = NP_in_vascular_cell;
	c->phenotype.secretion.saturation_densities[NP_index] = NP_in_vascular_cell;
}
*/

void vascular_cell_update_phenotype(Cell* c, Phenotype& phenotype, double dt) {
    //int voxel_index = c->get_current_voxel_index();
	//int oxygen_substrate_index = microenvironment.find_density_index("oxygen");
	//microenvironment.update_dirichlet_node(voxel_index, oxygen_substrate_index, 38.0);

    vascular_cell_handle_angiogen(c);    
    vascular_cell_handle_NP(c);
}


void create_vascular_cell_definition() {
    vascular_cell_def = cell_defaults; 
	vascular_cell_def.type = VASCULAR_CELL_TYPE; 
	vascular_cell_def.name = VASCULAR_CELL_NAME;
    vascular_cell_def.parameters.pReference_live_phenotype = &(vascular_cell_def.phenotype);
    
    // turn off motility
    vascular_cell_def.phenotype.motility.is_motile = false;

    // turn off cell division
    int cycle_start_index = live.find_phase_index(PhysiCell_constants::live); 
	int cycle_end_index = live.find_phase_index(PhysiCell_constants::live);
    vascular_cell_def.phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = 0.0; 

    // turn off cell death
    int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index("Apoptosis");
    vascular_cell_def.phenotype.death.rates[apoptosis_model_index] = 0.0;
    int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index("Necrosis");
    vascular_cell_def.phenotype.death.rates[necrosis_model_index] = 0.0;

    // set phenotype updates
    vascular_cell_def.functions.update_phenotype = vascular_cell_update_phenotype;

    // turn off volume updates
    vascular_cell_def.functions.volume_update_function = NULL;

    // set oxygen uptake and secretion parameters for vascular cells
    int oxygen_substrate_index = microenvironment.find_density_index("oxygen");
	vascular_cell_def.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 0;
	vascular_cell_def.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 38;
	vascular_cell_def.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 38;

    // the number of expansions
    vascular_cell_def.custom_data.add_variable("num_expansions", "dimensionless", 0);
    vascular_cell_def.custom_data.add_variable("flagged_for_division", "dimensionless", 0);

    // set NP uptake and secretion parameters
    //int NP_index = microenvironment.find_density_index("NP");
    //if (NP_index != -1) {
    //    vascular_cell_def.phenotype.secretion.uptake_rates[NP_index] = 0;
	//    vascular_cell_def.phenotype.secretion.secretion_rates[NP_index] = NP_in_vascular_cell;
	//    vascular_cell_def.phenotype.secretion.saturation_densities[NP_index] = NP_in_vascular_cell;
    //}
}


void init_vascular_cells() {
    init_vascular_params();
    create_vascular_cell_definition();
}


void setup_vascular_cell(Cell* c) {
	c->is_movable = false;
	
	//int voxel_index = c->get_current_voxel_index();
	//int oxygen_substrate_index = microenvironment.find_density_index("oxygen");
	//microenvironment.update_dirichlet_node(voxel_index, oxygen_substrate_index, 38.0);

	int ffd_index = c->custom_data.find_variable_index("flagged_for_division");
    c->custom_data[ffd_index] = 0;
    
	int numexp_index = c->custom_data.find_variable_index("num_expansions");
    c->custom_data[numexp_index] = 0;
}


Cell* create_vascular_cell(double x, double y, double z) {
	Cell* c = create_cell(vascular_cell_def);
	c->assign_position(x, y, z);
	setup_vascular_cell(c);
    vascular_cells->push_back(c);
    return c;
}


void divide_vascular_cell(Cell* c) {
	int ffd_index = c->custom_data.find_variable_index("flagged_for_division");
	c->custom_data[ffd_index] = 0;
    
	int numexp_index = c->custom_data.find_variable_index("num_expansions");
    c->custom_data[numexp_index]++;

	int angiogen_index = microenvironment.find_density_index("angiogen");
	std::vector<double> grad = c->nearest_gradient(angiogen_index);
	
	std::vector<double> rand_vec = {grad[0], grad[1], grad[2]}; 
	double l = sqrt(grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]);
	rand_vec /= l;
	rand_vec *= c->phenotype.geometry.radius;
	
	double x = c->position[0] + rand_vec[0];
	double y = c->position[1] + rand_vec[1];
	double z = c->position[2] + rand_vec[2];	

	create_vascular_cell(x, y, z);
}



//
// vascular network
//

double r, xmin, xmax, ymin, ymax, zmin, zmax;

int vascular_network_size = 0;
int num_branches = 0;

bool cell_in_cube(double x, double y, double z) {
    if (x - r < xmin || x + r > xmax)
        return false;

    if (y - r < ymin || y + r > ymax) 
        return false;

    if (z - r < zmin || z + r > zmax)
        return false;

    return true;
}


bool cell_in_cube(std::vector<double> coords) {
    return cell_in_cube(coords[0], coords[1], coords[2]);
}


void add_random_noise(std::vector<double>& direction) {
    double dx = UniformRandom() / 10;
    if (UniformRandom() < 0.5) 
        dx *= -1;

    double dy = UniformRandom() / 10;
    if (UniformRandom() < 0.5)
        dy *= -1;

    double dz = UniformRandom() / 10;
    if (UniformRandom() < 0.5)
        dz *= -1;

    direction[0] += dx;
    direction[1] += dy;
    direction[2] += dz;

    double l = sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);
	direction /= l;
}


void start_seeding(std::vector<double> seed, std::vector<double> direction, int branching_depth) {
    bool finished = false;
    int chain_len = 0;

    do {
        Cell* vascular = create_vascular_cell(seed[0], seed[1], seed[2]);
        ++vascular_network_size;
        ++chain_len;

        // initialize the expandability of the vascular cell
        if (chain_len > 1) {
            bool expandable_cell = UniformRandom() < expp_init_vasc;
            if (!expandable_cell) {
                int numexp_index = vascular->custom_data.find_variable_index("num_expansions");
                vascular->custom_data[numexp_index] = 1;
            }
        }

        // update the seed node        
        seed[0] += direction[0] * r;
        seed[1] += direction[1] * r;
        seed[2] += direction[2] * r;
        
        // do branching
        if (branching_depth < max_branching_depth) {
            bool branch = UniformRandom() < vascularity_branching_factor;
            if (branch) {
                ++num_branches;
                start_seeding(seed, UniformOnUnitSphere(), branching_depth + 1);
            }
        }

        // add a small amount of random noise to the direction vector
        add_random_noise(direction);

        finished = !cell_in_cube(seed) || chain_len >= vascular_chain_max_length;       
    } while (!finished); 
}

void refresh_globals() {
    r = vascular_cell_def.phenotype.geometry.radius;
    
    std::vector<double> borders = space_borders_for_random_sampling();
    xmin = borders[0]; xmax = borders[1];
    ymin = borders[2]; ymax = borders[3];
    zmin = borders[4]; zmax = borders[5]; 

    /*
    xmin = default_microenvironment_options.X_range[0] + default_microenvironment_options.expansion_delta + 10;
    xmax = default_microenvironment_options.X_range[1] - default_microenvironment_options.expansion_delta - 10;
    ymin = default_microenvironment_options.Y_range[0] + default_microenvironment_options.expansion_delta + 10;
    ymax = default_microenvironment_options.Y_range[1] - default_microenvironment_options.expansion_delta - 10;
    zmin = default_microenvironment_options.Z_range[0] + default_microenvironment_options.expansion_delta + 10;
    zmax = default_microenvironment_options.Z_range[1] - default_microenvironment_options.expansion_delta - 10;
    */
}

void init_vascular_network() {
    std::cout << "\nCreating vascularity network [";
    std::cout << "seed vascular cells = " << num_seed_points;
    std::cout << ", vascular chain max length = " << vascular_chain_max_length;
    std::cout << ", max branching depth = " << max_branching_depth;
    std::cout << ", branching factor = " << vascularity_branching_factor;
    std::cout << ", expansion probability = " << expp_init_vasc << "]" << std::endl;

    refresh_globals();
    
    SeedRandom(1UL);

    for (int i = 0; i < num_seed_points; i++) {
        double xa = (3 * xmin + xmax) / 4, xb = (3 * xmax + xmin) / 4;
        double ya = (3 * ymin + ymax) / 4, yb = (3 * ymax + ymin) / 4; 
        double za = (3 * zmin + zmax) / 4, zb = (3 * zmax + zmin) / 4;

        std::vector<double> seed = {uniformAB(xa, xb), uniformAB(ya, yb), uniformAB(za, zb)};
        std::vector<double> direction = UniformOnUnitSphere();

        start_seeding(seed, direction, 0); 
    }

    std::cout << "Finished [vascular network size = " << vascular_network_size;
    std::cout << ", #branches = " << num_branches << "]" << std::endl << std::endl;
}


void expand_vascular_network(double seed_multiplier) {
    std::cout << "Expanding vascularity network" << std::endl;

    refresh_globals();

    /*
    double xleft_start = default_microenvironment_options.X_range[0] + default_microenvironment_options.expansion_delta + 10;
    double xleft_end = default_microenvironment_options.X_range[0] + 2 * default_microenvironment_options.expansion_delta;
    double xright_start = default_microenvironment_options.X_range[1] - 2 * default_microenvironment_options.expansion_delta; 
    double xright_end = default_microenvironment_options.X_range[1] - default_microenvironment_options.expansion_delta - 10;

    double yleft_start = default_microenvironment_options.Y_range[0] + default_microenvironment_options.expansion_delta + 10;
    double yleft_end = default_microenvironment_options.Y_range[0] + 2 * default_microenvironment_options.expansion_delta;
    double yright_start = default_microenvironment_options.Y_range[1] - 2 * default_microenvironment_options.expansion_delta; 
    double yright_end = default_microenvironment_options.Y_range[1] - default_microenvironment_options.expansion_delta - 10;

    double zleft_start = default_microenvironment_options.Z_range[0] + default_microenvironment_options.expansion_delta + 10;
    double zleft_end = default_microenvironment_options.Z_range[0] + 2 * default_microenvironment_options.expansion_delta;
    double zright_start = default_microenvironment_options.Z_range[1] - 2 * default_microenvironment_options.expansion_delta; 
    double zright_end = default_microenvironment_options.Z_range[1] - default_microenvironment_options.expansion_delta - 10;
    */

    std::vector<double> borders = space_extension_for_random_sampling();

    SeedRandom(1UL);

    vascular_network_size = 0;
    num_branches = 0;

    for (int i = 0; i < num_seed_points * seed_multiplier; i++) {
        std::vector<double> seed = {
            /*
            uniformABCD(xleft_start, xleft_end, xright_start, xright_end), 
            uniformABCD(yleft_start, yleft_end, yright_start, yright_end), 
            uniformABCD(zleft_start, zleft_end, zright_start, zright_end)
            */

            uniformABCD(borders[0], borders[1], borders[2], borders[3]), 
            uniformABCD(borders[4], borders[5], borders[6], borders[7]), 
            uniformABCD(borders[8], borders[9], borders[10], borders[11])
        };

        std::vector<double> direction = UniformOnUnitSphere();
        
        // orient the direction vector towards the inner 
        if (seed[0] * direction[0] > 0)
            direction[0] *= -1;

        if (seed[1] * direction[1] > 0)
            direction[1] *= -1;

        if (seed[2] * direction[2] > 0)
            direction[2] *= -1;
        
        start_seeding(seed, direction, 0); 
    }

    std::cout << "Vasular expansion finished [#new vascular cells = " << vascular_network_size;
    std::cout << ", #branches = " << num_branches << "]" << std::endl << std::endl;
}