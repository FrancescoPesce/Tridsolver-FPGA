/**
 * Simulation space expansion
 * 
 * Author: Milos Savic (svc@dmi.uns.ac.rs)
 */

#ifndef __spexp_h__
#define __spexp_h__

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"

using namespace BioFVM; 
using namespace PhysiCell;

#include "evo_nano_3D.h"

void expand_microenvironment(double mechanics_voxel_size);
bool check_conditions_to_expand();

#endif