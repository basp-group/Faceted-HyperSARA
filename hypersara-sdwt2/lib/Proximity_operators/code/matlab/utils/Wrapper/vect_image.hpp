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

#ifndef VECT_IMAGE_HPP
#define VECT_IMAGE_HPP

#include "image.hpp"

template<typename T>
struct VectImage : public Image<T>
{
	typedef typename std::vector<T>::iterator iterator;
	typedef typename std::vector<T>::const_iterator const_iterator;

    std::vector<T> vect;

    /* constructors */
    VectImage() {
    	create(0,0,0,0);
    }
    VectImage(int N_rows, int N_cols, int N_bands = 1, int N_planes = 1) {
    	create(N_rows, N_cols, N_bands, N_planes);
    }
    VectImage(const int* sz) {
    	create(sz);
    }
    VectImage(const VectImage &that) {
    	create( that.size() );
    	std::copy( that.begin(), that.end(), this->begin() );
    }
    
    void create(const int* sz) {
    	create(sz[0], sz[1], sz[2], sz[3]);
    }
    void create(int N_rows, int N_cols, int N_bands = 1, int N_planes = 1)
    {
        vect.resize(N_rows*N_cols*N_bands*N_planes, T());

        // data pointer
        Image<T>::data = vect.data();

        // matrix size
        Image<T>::size_[0] = N_rows;
        Image<T>::size_[1] = N_cols;
        Image<T>::size_[2] = N_bands;
        Image<T>::size_[3] = N_planes;

        // roi handling
        Image<T>::step[0] = Image<T>::size(0);
        Image<T>::step[1] = Image<T>::size(0)*Image<T>::size(1);
        Image<T>::step[2] = Image<T>::size(0)*Image<T>::size(1)*Image<T>::size(2);
    }

    /* iterators*/
    typename std::vector<T>::const_iterator begin() const {
    	return vect.begin();
    }
    typename std::vector<T>::iterator begin() {
    	return vect.begin();
    }
    typename std::vector<T>::const_iterator end() const {
    	return vect.end();
    }
    typename std::vector<T>::iterator end() {
    	return vect.end();
    }
};


#endif