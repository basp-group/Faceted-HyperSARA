/*
 *  This mex-file wraps the projection onto the simplex proposed by 
 *  L. Condat in "Fast projection onto the simplex and the l1 ball". 
 * 
 *
 *  COMPILATION
 * =============
 *   1. Download the file "condat_simplexproj.c" from Condat's web page:
 *
 *        https://www.gipsa-lab.grenoble-inp.fr/~laurent.condat/download/condat_simplexproj.c
 *
 *   2. Open the file and COMMENT OUT these lines:
 * 
 *         87 #include <gsl/gsl_rng.h>
 *         88 #include <gsl/gsl_randist.h>
 *
 *      as well as the fonction "main()":
 *
 *        572 int main() {
 *           ...
 *        620 }
 *
 *   3. Compile the mex-file as follows:
 *
 *        mex project_simplex_condat.cpp
 *
 *  A C++ compiler is required.
 *
 *
 *  USAGE
 * =======
 *
 *  >> eta = 1;
 *  >> x = rand(10);
 *  >> y = project_simplex_condat(x, eta);
 *
 **************************************************************************
 * Version : 1.0 (27-04-2017)
 * Author  : Giovanni Chierchia
 **************************************************************************
 * Copyright (C) 2017
 *
 * This file is part of the codes provided at http://proximity-operator.net
 *
 * By downloading and/or using any of these files, you implicitly agree to 
 * all the terms of the license CeCill-B (available online).
 **************************************************************************
 */


#include "Wrapper/mx_image.hpp"
#include "condat_simplexproj.c"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // check the params
    if (nrhs != 2)
        mexErrMsgTxt("This function takes 2 inputs");
    if (nlhs != 1)
        mexErrMsgTxt("This function gives 1 output");
       
    // read the inputs
    MxImage<double> x( prhs[0] );
    const double a      = (double) mxGetScalar( prhs[1] );
    
    // init. the output
    int N1 = x.size(0);
    int N2 = x.size(1);
    int N3 = x.size(2);
    int N4 = x.size(3);
    MxImage<double> y( plhs[0], N1, N2, N3, N4 );
       
    // write the output
    simplexproj_Condat( x.begin(), y.begin(), x.total(), a );
}