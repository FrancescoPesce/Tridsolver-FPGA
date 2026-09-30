/**
 * Evo-Nano vascularity
 * 
 * Author: Milos Savic (svc@dmi.uns.ac.rs)
 */

#ifndef __vascular_h__
#define __vascular_h__

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"

#include "utils.h"

using namespace BioFVM; 
using namespace PhysiCell;

const std::string VASCULAR_CELL_NAME = "vascular cell";
const int VASCULAR_CELL_TYPE = 3;

// creates the initial vasularity network
void init_vascular_network();

// expands the vascularity network when the simulation space is expanded
void expand_vascular_network(double seed_multiplier);    

// setup vascular cell definition
void init_vascular_cells();

// create a new vascular cell at the given position
Cell* create_vascular_cell(double x, double y, double z);

// "divide" a vascular cell (create a new vascular cell extending the given cell)
void divide_vascular_cell(Cell* c);

// is NP active when simulated?
extern bool NP_active;

// is NP simulated at all?
extern bool NP_simul;

// vector containing all vascular cells
extern std::vector<Cell*>* vascular_cells;

#endif