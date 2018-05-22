#ifndef VOR_DETECTION_HPP
#define VOR_DETECTION_HPP

#include "stdinclude.hpp"
#include "rw_restart.hpp"

MatrixXd Vortex_detection ( unsigned int Nr, int Nc, vector<int> col_grads, string filename ) ;

#endif