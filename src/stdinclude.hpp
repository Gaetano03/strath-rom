
#ifndef STDINCLUDE_HPP
#define STDINCLUDE_HPP
// #define PI  3.14159265358979323846 


#include "ctime"
#include <iostream>
#include <stdio.h> 
#include <vector>
#include <fstream>
#include <sstream>
#include <complex> 
#include <cmath>
#include <iomanip>
#include <string>
#include "smartmath.h"
#include "smartuq.h"
#include "LinearAlgebra/Eigen/Dense"
#include "LinearAlgebra/Eigen/Eigenvalues"
#include "LinearAlgebra/Eigen/Core"
#include "LinearAlgebra/Eigen/Householder"

using namespace std;
using namespace Eigen;

int fibonacci( int i){

  if (i < 0) return -1; /* F(i) non e' definito per interi i negativi! */

  if (i == 0) return 0;
  else if (i == 1) return 1;
  else return fibonacci(i-1) + fibonacci(i-2);

}



void eig_sort(VectorXd &eig_val, MatrixXd &eig_vec)
{
    unsigned int swap_count = 1;
    double temp;
    VectorXd temp_vec(eig_vec.rows());

    while (swap_count > 0)
    {
        swap_count = 0;

        for(unsigned int index = 1; index < eig_val.size(); index++)
        {
            if (eig_val(index) > eig_val(index-1))
            {
                temp = eig_val(index-1);
                eig_val(index-1) = eig_val(index);
                eig_val(index) = temp;

                temp_vec = eig_vec.col(index-1);
                eig_vec.col(index-1) = eig_vec.col(index);
                eig_vec.col(index) = temp_vec;

                swap_count++;
            }
        }
    }
}



#endif