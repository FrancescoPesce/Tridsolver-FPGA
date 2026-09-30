/**
 * EvoNano main module
 * 
 * Author: Milos Savic (svc@dmi.uns.ac.rs)
 */

#include <fstream>
#include <sstream>
#include <string>

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"

#include "cancc.h"
#include "vascular.h"
#include "healthyc.h"
#include "utils.h"

using namespace BioFVM; 
using namespace PhysiCell;

// setup functions 
void create_cell_types( void );
void setup_tissue( void );
void init_params_and_cell_definitions();
void setup_microenvironment( void ); 

// custom pathology coloring function 
std::vector<std::string> my_coloring_function(Cell*);

// periodic reporting about cells
void cells_report(void);

// EvoNano cell division rules
void evonano_cell_division(Cell_Container* cell_container, double t);

// EvoNano CSC detachment rules
void evonano_CSC_deatachment(Cell_Container* cell_container, double t);



