/**************************************************************************
 * Version : 1.0 (27-04-2017)
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
#include "Eigen/Eigen"
using namespace Eigen;


template <typename T>
void sv_rec(Image<T> &y, const Image<T> &u, const Image<T> &s, const Image<T> &v);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // check the params
    if (nrhs != 3)
        mexErrMsgTxt("This function takes 3 input");
    if (nlhs != 1)
        mexErrMsgTxt("This function gives 1 output");
    
    // pixel type
    typedef double T;
    
    // read the inputs
    MxImage<T> u( prhs[0] );
    MxImage<T> s( prhs[1] );
    MxImage<T> v( prhs[2] );
    
    // init. the outputs
    int N1 = u.size(0);
    int N2 = v.size(0);
    int N3 = v.size(2);
    int N4 = v.size(3);
    MxImage<T> y( plhs[0], N1, N2, N3, N4 );

    // write the output
    sv_rec(y, u, s, v);
}


template <typename T>
void sv_rec(Image<T> &y, const Image<T> &u, const Image<T> &s, const Image<T> &v)
{
    // get the size
    int N1 = y.size(0);
    int N2 = y.size(1);
    int B = std::min(N1,N2);
    int N = y.size(2)*y.size(3);
    
    // scroll the field
    for(int b=0; b < N; b++) 
    {
        // pointers
        T* p_y = y.ptr(0,0,b);
        const T* p_u = u.ptr(0,0,b);
        const T* p_s = s.ptr(0,b);
        const T* p_v = v.ptr(0,0,b);
        
        // wrapping for Eigen
        Map<MatrixXd> yy(p_y, N1, N2);
        Map<const MatrixXd> uu(p_u, N1, B);
        Map<const VectorXd> ss(p_s, B);
        Map<const MatrixXd> vv(p_v, N2, B);
    
        // s.v. reconstruction
        yy = uu * ss.asDiagonal() * vv.transpose();
    }
}