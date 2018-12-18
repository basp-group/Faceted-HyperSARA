/*
 * This mex-file implements the Pool Adjacent Violators Algorithm, which 
 * finds in O(N) the solution to the following optimization problem:
 *
 *  minimize_y  sum_i w_i (y_i - x_i)^2   s.t.  y_1 <= ... <= y_N
 *
 * where 'x' and 'w' are input arrays of the same size. When the inputs 
 * are matrices, the above problem is solved column-by-colum.
 *
 *  USAGE
 * =======
 *
 * x = 50 * randn(8,2);
 * w = 10 * rand(size(x));
 * y = pava(x,w);
 * 
 *************************************************************************
 * Version : 1.0 (07-06-2017)
 * Author  : Giovanni Chierchia
 **************************************************************************
 * Copyright (C) 2017
 *
 * This file is part of the codes provided at http://proximity-operator.net
 *
 * By downloading and/or using any of these files, you implicitly agree to 
 * all the terms of the license CeCill-B (available online).
 **************************************************************************/

#include "Wrapper/mx_image.hpp"


template <typename T>
void pav(Image<T> &y, const Image<T> &x, const Image<T> &w);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // check the params
    if (nrhs != 2)
        mexErrMsgTxt("This function takes 2 input");
    if (nlhs != 1)
        mexErrMsgTxt("This function gives 1 output");
    
    // pixel type
    typedef double T;
    
    // read the inputs
    MxImage<T> x( prhs[0] );
    MxImage<T> w( prhs[1] );
    
    // check the inputs
    if ( x.rows() != w.rows() || x.cols() != w.cols() )
        mexErrMsgTxt("Inputs must have the same size");
    
    // init the outputs
    int N1 = x.rows();
    int N2 = x.cols();
    MxImage<T> y( plhs[0], N1, N2 );

    // write the output
    pav(y, x, w);
}


template <typename T>
void pav(Image<T> &y, const Image<T> &x, const Image<T> &w)
{
    // sizes
    int N = y.rows(); // length of each block
    int L = y.cols(); // number of blocks
    
    // trivial case
    if (N <= 1) {
        x.copy_to(y);
        return;
    }
    
    // scrolling columns
    for(int l=0; l<L; l++)
    {
        // pointers
              T* p_y = y.ptr(0,l);
        const T* p_x = x.ptr(0,l);
        const T* p_w = w.ptr(0,l);

        // copy
        std::copy(p_x, p_x+N, p_y);

        // loop
        bool pooled;
        do {
            int i = 0;
            pooled = false;
            while (i < N-1) {
                int k = i;    
                while (k < N-1 && p_y[k] >= p_y[k+1]) k += 1;
                if (p_y[i] != p_y[k]) 
                {
                    double num = 0.0;
                    double den = 0.0;
                    for (int j = i; j <= k; j++) {
                        num += p_w[j] * p_y[j];
                        den += p_w[j];
                    }
                    for (int j = i; j <= k; j++) {
                        p_y[j] = num / den;
                    }
                    pooled = true;
                }
                i = k+1;
            }
        } while (pooled); 
    }
}