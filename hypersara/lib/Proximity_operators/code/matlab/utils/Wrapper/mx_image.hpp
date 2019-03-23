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

#ifndef MX_MEX_HPP
#define MX_MEX_HPP

#include "image.hpp"
#include "mex_traits.hpp"

template <typename T>
struct MxImage : public Image<T>
{
    mxArray* owned_mx;
            
    /*
     * This constructor wraps a "previously allocated" mxArray
     */
    MxImage(const mxArray* mx) : owned_mx(nullptr) {
        super(mx);
    }

    /*
     * This constructor allocates and wraps a new mxArray. 
     * The pointer is returned in 'mx'. 
     * The caller takes care of the deallocation.
     */
    MxImage(mxArray*& mx, int N_rows, int N_cols, int N_bands = 1, int N_planes = 1) : owned_mx(nullptr) {
        const mwSize sz[4] = {N_rows, N_cols, N_bands, N_planes};
        mx = mxCreateNumericArray(4, sz, mxClassIDTrait<T>::type, mxREAL);
        super(mx);
    }
       
    /*
     * This constructor allocates and wraps a new mxArray. 
     * The pointer is not returned, as the destructor handles the deallocation.
     */
    MxImage(int N_rows, int N_cols, int N_bands = 1, int N_planes = 1) : owned_mx(nullptr) { 
         const mwSize sz[4] = {N_rows, N_cols, N_bands, N_planes};
         owned_mx = mxCreateNumericArray(4, sz, mxClassIDTrait<T>::type, mxREAL);
         super(owned_mx);
    }
    
    /*
     * This constructor allocates and wraps a new mxArray. 
     * The pointer is not returned, as the destructor handles the deallocation.
     */
    MxImage(const int* size) : owned_mx(nullptr) { 
         const mwSize sz[4] = {size[0], size[1], size[2], size[3]};
         owned_mx = mxCreateNumericArray(4, sz, mxClassIDTrait<T>::type, mxREAL);
         super(owned_mx);
    }
     
     ~MxImage() {
         if( owned_mx != nullptr ) mxDestroyArray(owned_mx);
     }

private:

    void super(const mxArray *mx)
    {        
        // complex checking
        if( mxIsComplex(mx) )
            mexErrMsgTxt("Complex values not supported");
        
        // type checking
        if( mxClassIDTrait<T>::type != mxGetClassID(mx) )
            mexErrMsgTxt("The template type does not match the matrix typed");
            
        // size checking
        mwSize Ndim = mxGetNumberOfDimensions(mx);
        if(Ndim > 4)
            mexErrMsgTxt("The array cannot have more than 4 dimensions");

        // data pointer
        T* data = (T*) mxGetData(mx);

        // matrix size
        const mwSize* sz  = mxGetDimensions(mx);
        int size[4];
        size[0] = (Ndim>0) ? sz[0] : 1;
        size[1] = (Ndim>1) ? sz[1] : 1;
        size[2] = (Ndim>2) ? sz[2] : 1;
        size[3] = (Ndim>3) ? sz[3] : 1;

        // roi handling
        size_t step[3];
        step[0] = size[0];
        step[1] = size[0]*size[1];
        step[2] = size[0]*size[1]*size[2];
        
        // configure the base class
        Image<T>::init(data, size, step);
    }
};


template <>
struct MxImage<mxLogical> : public Image<mxLogical>
{   
    mxArray* owned_mx;
            
    /*
     * This constructor wraps a "previously allocated" mxArray
     */
    MxImage(const mxArray* mx) : owned_mx(nullptr) {
        super(mx);
    }

    /*
     * This constructor allocates and wraps a new mxArray. 
     * The pointer is returned in 'mx'. 
     * The caller takes care of the deallocation.
     */
    MxImage(mxArray*& mx, int N_rows, int N_cols, int N_bands = 1, int N_planes = 1) : owned_mx(nullptr) {
        const mwSize sz[4] = {N_rows, N_cols, N_bands, N_planes};
        mx = mxCreateLogicalArray(4, sz);
        super(mx);
    }
       
    /*
     * This constructor allocates and wraps a new mxArray. 
     * The pointer is not returned, as the destructor handles the deallocation.
     */
    MxImage(int N_rows, int N_cols, int N_bands = 1, int N_planes = 1) : owned_mx(nullptr) { 
         const mwSize sz[4] = {N_rows, N_cols, N_bands, N_planes};
         owned_mx = mxCreateLogicalArray(4, sz);
         super(owned_mx);
    }
    
    /*
     * This constructor allocates and wraps a new mxArray. 
     * The pointer is not returned, as the destructor handles the deallocation.
     */
    MxImage(const int* size) : owned_mx(nullptr) { 
         const mwSize sz[4] = {size[0], size[1], size[2], size[3]};
         owned_mx = mxCreateLogicalArray(4, sz);
         super(owned_mx);
    }
     
     ~MxImage() {
         if( owned_mx != nullptr ) mxDestroyArray(owned_mx);
     }

private:

    void super(const mxArray *mx)
    {                
        // type checking
        if( mxLOGICAL_CLASS != mxGetClassID(mx) )
            mexErrMsgTxt("The matrix is not logical");
            
        // size checking
        mwSize Ndim = mxGetNumberOfDimensions(mx);
        if(Ndim > 4)
            mexErrMsgTxt("The array cannot have more than 4 dimensions");

        // data pointer
        mxLogical* data = mxGetLogicals(mx);

        // matrix size
        const mwSize* sz  = mxGetDimensions(mx);
        int size[4];
        size[0] = (Ndim>0) ? sz[0] : 1;
        size[1] = (Ndim>1) ? sz[1] : 1;
        size[2] = (Ndim>2) ? sz[2] : 1;
        size[3] = (Ndim>3) ? sz[3] : 1;

        // roi handling
        size_t step[3];
        step[0] = size[0];
        step[1] = size[0]*size[1];
        step[2] = size[0]*size[1]*size[2];
        
        // configure the base class
        Image<mxLogical>::init(data, size, step);
    }
};


#endif