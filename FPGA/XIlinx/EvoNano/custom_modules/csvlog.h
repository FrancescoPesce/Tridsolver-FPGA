/**
 * CSV cell log
 * 
 * Author: Milos Savic (svc@dmi.uns.ac.rs)
 */

#ifndef __csvlog_h__
#define __csvlog_h__

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"

using namespace BioFVM;
using namespace PhysiCell;

#include "evo_nano_3D.h"

class CSV_log {
public:
    CSV_log(std::string file_name);
    void log();
private:
    std::string file_name;
    int hist_dead[5];
};

#endif