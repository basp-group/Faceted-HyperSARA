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
void sv_dec(Image<T> &u, Image<T> &s, Image<T> &v, const Image<T> &x);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // check the params
    if (nrhs != 1)
        mexErrMsgTxt("This function takes 1 input");
    if (nlhs != 3)
        mexErrMsgTxt("This function gives 3 output");
    
    // pixel type
    typedef double T;
    
    // read the inputs
    MxImage<T> x( prhs[0] );
    
    // init. the outputs
    int N1 = x.size(0);
    int N2 = x.size(1);
    int N3 = x.size(2);
    int N4 = x.size(3);
    int M = std::min(N1,N2);
    MxImage<T> u( plhs[0], N1, M, N3, N4 );
    MxImage<T> s( plhs[1],  M, N3, N4 );
    MxImage<T> v( plhs[2], N2, M, N3, N4 );

    // write the output
    sv_dec(u, s, v, x);
}


template <typename T>
void sv_dec(Image<T> &u, Image<T> &s, Image<T> &v, const Image<T> &x)
{
    // get the size
    int N1 = x.size(0);
    int N2 = x.size(1);
    int B = std::min(N1,N2);
    int N = x.size(2)*x.size(3);
    
    // scroll the field
    for(int b=0; b < N; b++) 
    {
        // pointers
        T* p_u = u.ptr(0,0,b);
        T* p_s = s.ptr(0,b);
        T* p_v = v.ptr(0,0,b);
        const T* p_x = x.ptr(0,0,b);
        
            // wrapping for Eigen
        Map<MatrixXd> uu(p_u, N1, B);
        Map<VectorXd> ss(p_s, B);
        Map<MatrixXd> vv(p_v, N2, B);
        Map<const MatrixXd> xx(p_x, N1, N2);
        
        // s.v. decomposition
        JacobiSVD<MatrixXd,HouseholderQRPreconditioner> svd(xx, ComputeThinU | ComputeThinV);
        
        // write the output
        uu = svd.matrixU();
        ss = svd.singularValues();
        vv = svd.matrixV();
    }
}