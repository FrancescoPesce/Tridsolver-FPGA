/**
 * Evo-Nano shared utility functions
 * 
 * Author: Milos Savic (svc@dmi.uns.ac.rs)
 */

#ifndef __utils_h__
#define __utils_h__

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"

using namespace BioFVM; 
using namespace PhysiCell;

double uniformAB(double a, double b);
double uniformABCD(double a, double b, double c, double d);

std::vector<double> space_borders_for_random_sampling();
std::vector<double> space_extension_for_random_sampling();

#endif