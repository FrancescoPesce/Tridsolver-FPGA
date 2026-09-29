/**
 * Evo-Nano shared utility functions
 * 
 * Author: Milos Savic (svc@dmi.uns.ac.rs)
 */

#include "utils.h"


double uniformAB(double a, double b) {
    return a + UniformRandom() * (b - a);
}


double uniformABCD(double a, double b, double c, double d) {
    bool left = UniformRandom() < 0.5;
    if (left)
        return uniformAB(a, b);
    else
        return uniformAB(c, d);
}


std::vector<double> space_borders_for_random_sampling() {
    std::vector<double> ret;

    // xmin, xmax
    ret.push_back(default_microenvironment_options.X_range[0] + default_microenvironment_options.expansion_delta + 10);
    ret.push_back(default_microenvironment_options.X_range[1] - default_microenvironment_options.expansion_delta - 10);
    
    // ymin, ymax
    ret.push_back(default_microenvironment_options.Y_range[0] + default_microenvironment_options.expansion_delta + 10);
    ret.push_back(default_microenvironment_options.Y_range[1] - default_microenvironment_options.expansion_delta - 10);
    
    // zmin, zmax
    ret.push_back(default_microenvironment_options.Z_range[0] + default_microenvironment_options.expansion_delta + 10);
    ret.push_back(default_microenvironment_options.Z_range[1] - default_microenvironment_options.expansion_delta - 10);

    return ret;
}


std::vector<double> space_extension_for_random_sampling() {
    std::vector<double> ret;
    
    // xleft_start, xleft_end, xright_start, xright_end
    ret.push_back(default_microenvironment_options.X_range[0] + default_microenvironment_options.expansion_delta + 10);
    ret.push_back(default_microenvironment_options.X_range[0] + 2 * default_microenvironment_options.expansion_delta);
    ret.push_back(default_microenvironment_options.X_range[1] - 2 * default_microenvironment_options.expansion_delta); 
    ret.push_back(default_microenvironment_options.X_range[1] - default_microenvironment_options.expansion_delta - 10);

    // yleft_start, yleft_end, yright_start, yright_end
    ret.push_back(default_microenvironment_options.Y_range[0] + default_microenvironment_options.expansion_delta + 10);
    ret.push_back(default_microenvironment_options.Y_range[0] + 2 * default_microenvironment_options.expansion_delta);
    ret.push_back(default_microenvironment_options.Y_range[1] - 2 * default_microenvironment_options.expansion_delta); 
    ret.push_back(default_microenvironment_options.Y_range[1] - default_microenvironment_options.expansion_delta - 10);

    // zleft_start, zleft_end, zright_start, zright_end
    ret.push_back(default_microenvironment_options.Z_range[0] + default_microenvironment_options.expansion_delta + 10);
    ret.push_back(default_microenvironment_options.Z_range[0] + 2 * default_microenvironment_options.expansion_delta);
    ret.push_back(default_microenvironment_options.Z_range[1] - 2 * default_microenvironment_options.expansion_delta); 
    ret.push_back(default_microenvironment_options.Z_range[1] - default_microenvironment_options.expansion_delta - 10);

    return ret;
}