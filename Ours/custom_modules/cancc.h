/**
 * Evo-Nano cancer cells
 * 
 * Author: Milos Savic (svc@dmi.uns.ac.rs) & Vladimir Kurbalija (Phase 2)
 */

#ifndef __cancc_h__
#define __cancc_h__

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"

#include "vascular.h"

using namespace BioFVM; 
using namespace PhysiCell;

const std::string CANCER_CELL_NAME = "cancer cell";
const int CANCER_CELL_TYPE = 1;

const std::string CANCER_STEM_CELL_NAME = "cancer stem cell";
const int CANCER_STEM_CELL_TYPE = 2;

const std::string CAF_CELL_NAME = "CAF";
const int CAF_CELL_TYPE = 5;


void init_cancer_cells();

Cell* create_cancer_cell(double x, double y, double z, bool stem);

Cell* create_CAF_cell(double x, double y, double z);

void divide_cancer_cell(Cell* c);
void divide_cancer_stem_cell(Cell* c);

extern int NP_CC_deaths;

// the number of CC converted to CSC based on cytokine 
extern int cytokine_conversions;

// are cytokine and CAF cells simulated
extern bool cytokine_caf_simul;

// is prostaglandin simulated at all?
extern bool prostaglandin_simul;

// the number of CC converted to CSC based on prostaglandin 
extern int prostaglandin_conversions;

#endif
