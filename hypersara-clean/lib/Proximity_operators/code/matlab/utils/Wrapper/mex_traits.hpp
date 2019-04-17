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

#ifndef MEX_TRAITS_HPP
#define MEX_TRAITS_HPP

#include <mex.h>
#include <matrix.h>


template <typename T>
struct mxClassIDTrait {
	static const mxClassID type = mxUNKNOWN_CLASS;
};
template <>
struct mxClassIDTrait<mxLogical> {
	static const mxClassID type = mxLOGICAL_CLASS;
};
template <>
struct mxClassIDTrait<uint8_t> {
	static const mxClassID type = mxUINT8_CLASS;
};
template <>
struct mxClassIDTrait<int8_t> {
	static const mxClassID type = mxINT8_CLASS;
};
template <>
struct mxClassIDTrait<uint16_t> {
	static const mxClassID type = mxUINT16_CLASS;
};
template <>
struct mxClassIDTrait<int16_t> {
	static const mxClassID type = mxINT16_CLASS;
};
template <>
struct mxClassIDTrait<uint32_t> {
	static const mxClassID type = mxUINT32_CLASS;
};
template <>
struct mxClassIDTrait<int32_t> {
	static const mxClassID type = mxINT32_CLASS;
};
template <>
struct mxClassIDTrait<uint64_t> {
	static const mxClassID type = mxUINT64_CLASS;
};
template <>
struct mxClassIDTrait<int64_t> {
	static const mxClassID type = mxINT64_CLASS;
};
template <>
struct mxClassIDTrait<float> {
	static const mxClassID type = mxSINGLE_CLASS;
};
template <>
struct mxClassIDTrait<double> {
	static const mxClassID type = mxDOUBLE_CLASS;
};

#endif