/**
 * EvoNano healthy cells
 * 
 * Author: Milos Savic (svc@dmi.uns.ac.rs)
 */

#ifndef __healthyc_h__
#define __healthyc_h__

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"

#include "utils.h"
#include "vascular.h"

using namespace BioFVM; 
using namespace PhysiCell;

const std::string HEALTHY_CELL_NAME = "healthy cell";
const int HEALTHY_CELL_TYPE = 4;


void init_healthy_cells();
Cell* create_healthy_cell(double x, double y, double z);

extern int NP_HC_deaths;

void expand_healthy_cells(int num_expansions);

#endif