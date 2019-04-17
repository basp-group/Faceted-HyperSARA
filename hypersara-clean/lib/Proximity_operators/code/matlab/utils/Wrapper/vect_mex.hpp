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

#ifndef VECT_MEX_HPP_
#define VECT_MEX_HPP_

#include "vect_image.hpp"
#include "mex_traits.hpp"
#include <stdint.h>


template <typename T>
VectImage<T> mx2im(const mxArray* mx)
{
	switch( mxGetClassID(mx) )
	{
        case mxLOGICAL_CLASS: return mexim::read_logical<T>(mx);
        case mxUINT8_CLASS:	  return mexim::read_numeric<uint8_t, T>(mx);
		case mxINT8_CLASS: 	  return mexim::read_numeric< int8_t, T>(mx);
		case mxUINT16_CLASS:  return mexim::read_numeric<uint16_t,T>(mx);
		case mxINT16_CLASS:   return mexim::read_numeric< int16_t,T>(mx);
		case mxUINT32_CLASS:  return mexim::read_numeric<uint32_t,T>(mx);
		case mxINT32_CLASS:   return mexim::read_numeric< int32_t,T>(mx);
		case mxUINT64_CLASS:  return mexim::read_numeric<uint64_t,T>(mx);
		case mxINT64_CLASS:	  return mexim::read_numeric< int64_t,T>(mx);
		case mxSINGLE_CLASS:  return mexim::read_numeric<   float,T>(mx);
		case mxDOUBLE_CLASS:  return mexim::read_numeric<  double,T>(mx);

		// mxUNKNOWN_CLASS, mxCELL_CLASS, mxSTRUCT_CLASS, mxCHAR_CLASS, mxVOID_CLASS, mxFUNCTION_CLASS, mxOPAQUE_CLASS, mxOBJECT_CLASS
		default: mexErrMsgTxt("[mx2im] mx_type not expected.");
	}
}


template<typename T>
mxArray* im2mx(const Image<T> &im, mxClassID mx_type)
{
    switch(mx_type)
    {
        case mxLOGICAL_CLASS: return mexim::write_logical<T>(im);
        case mxUINT8_CLASS:	  return mexim::write_numeric<T, uint8_t>(im);
		case mxINT8_CLASS: 	  return mexim::write_numeric<T,  int8_t>(im);
		case mxUINT16_CLASS:  return mexim::write_numeric<T,uint16_t>(im);
		case mxINT16_CLASS:   return mexim::write_numeric<T, int16_t>(im);
		case mxUINT32_CLASS:  return mexim::write_numeric<T,uint32_t>(im);
		case mxINT32_CLASS:   return mexim::write_numeric<T, int32_t>(im);
		case mxUINT64_CLASS:  return mexim::write_numeric<T,uint64_t>(im);
		case mxINT64_CLASS:	  return mexim::write_numeric<T, int64_t>(im);
		case mxSINGLE_CLASS:  return mexim::write_numeric<T,   float>(im);
		case mxDOUBLE_CLASS:  return mexim::write_numeric<T,  double>(im);

		// mxUNKNOWN_CLASS, mxCELL_CLASS, mxSTRUCT_CLASS, mxCHAR_CLASS, mxVOID_CLASS, mxFUNCTION_CLASS, mxOPAQUE_CLASS, mxOBJECT_CLASS
		default: mexErrMsgTxt("[mx2im] mx_type not expected.");
	}
}


template <typename T> 
mxArray* im2mx(const Image<T> &im) {
    return im2mx(im, mexim::mxClassIDTrait<T>::type);
}



namespace mexim {
    
template <typename mx_type, typename T>
VectImage<T> read_numeric(const mxArray* mx)
{
    // matrix size
    mwSize Ndim = mxGetNumberOfDimensions(mx);
    const mwSize* sz  = mxGetDimensions(mx);
    int size[4];
    size[0] = (Ndim>0) ? sz[0] : 1;
    size[1] = (Ndim>1) ? sz[1] : 1;
    size[2] = (Ndim>2) ? sz[2] : 1;
    size[3] = (Ndim>3) ? sz[3] : 1;

    // check the size
    if(Ndim > 4)
    	mexErrMsgTxt("The array cannot have more than 4 dimensions");

    // check the type
    if( mxIsComplex(mx) )
        mexErrMsgTxt("Complex values are not supported (yet)");

    // create the matrix
    VectImage<T> im(size);

    // fill the matrix
    T* im_ptr = im.ptr(0);
    const mx_type* mx_ptr = (mx_type*) mxGetData(mx);
    for(int n=0; n < im.total(); n++) {
    	im_ptr[n] = (T) mx_ptr[n];
    }
    
    return im;
}


template <typename typename T>
VectImage<T> read_logical(const mxArray* mx)
{
    // matrix size
    mwSize Ndim = mxGetNumberOfDimensions(mx);
    const mwSize* sz  = mxGetDimensions(mx);
    int size[4];
    size[0] = (Ndim>0) ? sz[0] : 1;
    size[1] = (Ndim>1) ? sz[1] : 1;
    size[2] = (Ndim>2) ? sz[2] : 1;
    size[3] = (Ndim>3) ? sz[3] : 1;

    // check the size
    if(Ndim > 4)
    	mexErrMsgTxt("The array cannot have more than 4 dimensions");

    // create the matrix
    VectImage<T> im(size);

    // fill the matrix
    T* im_ptr = im.ptr(0);
    const mxLogical* mx_ptr = mxGetLogicals(mx);
    for(int n=0; n < im.total(); n++) {
    	im_ptr[n] = (T) mx_ptr[n];
    }
    
    return im;
}


template<typename T, typename mx_type>
mxArray* write_numeric(const Image<T> &im)
{
    // get the size
    const int* sz = im.size();

    // create the matrix
    mxArray* mx = mxCreateNumericArray(4, sz, mxClassIDTrait<mx_type>::type, mxREAL);

    // fill the matrix
    const T* im_ptr = im.ptr(0);
    mx_type* mx_ptr = (mx_type*) mxGetData(mx);
    for(int n=0; n < im.total(); n++) {
        mx_ptr[n] = (mx_type) im_ptr[n];
    }

    return mx;
}


template<typename T>
mxArray* write_logical(const Image<T> &im)
{
    // get the size
    const int* sz = im.size();

    // create the matrix
    mxArray* mx = mxCreateLogicalArray(4, sz);

    // fill the matrix
    const T* im_ptr = im.ptr(0);
    mxLogical* mx_ptr = mxGetLogicals(mx);
    for(int n=0; n < im.total(); n++) {
        mx_ptr[n] = (mxLogical) im_ptr[n];
    }

    return mx;
}

}


#endif /* MY_MEX_HPP_ */
