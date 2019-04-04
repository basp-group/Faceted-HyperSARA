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


#ifndef EIGEN_PLUGIN_HPP
#define EIGEN_PLUGIN_HPP


#include <Eigen/Core>
using Eigen::MatrixXd;
using Eigen::Map;
using Eigen;    // to be commented in case of conflicts

//
// NOTE: Eigen uses a column-major order, like matlab !!!
//

template <typename T>
void convert_to_eigen(MatrixXd &y, const Image<T> &x)
{
    // map the input
    Map<const MatrixXd> xx( x.ptr(0), x.rows(), x.cols() );
    
    // copy (the resizing is automatic)
    y = xx.cast<double>();
}

template <typename T>
void convert_from_eigen(Image<T> &y, const MatrixXd &x)
{
    if( x.rows() != y.rows() || x.cols() != y.cols() )
        mexErrMsgTxt("conversion error");
    
    // map the output
    Map<MatrixXd> yy( y.ptr(0), x.rows(), x.cols() );
    
    // copy (the resizing is NOT automatic)
    yy = x.cast<T>();
}

#endif