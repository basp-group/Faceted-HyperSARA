"""
 Version : 1.0 (27-04-2017)
 
 Author  : Giovanni Chierchia

 Copyright (C) 2017
 
 This file is part of the codes provided at http://proximity-operator.net
 
 By downloading and/or using any of these files, you implicitly agree to 
 all the terms of the license CeCill-B (available online).
 """
   
# -*- coding: utf-8 -*-

import numpy as np;

class Bathtub :
    """
    This class computes the proximity operator of the function

            f(x) = gamma * max{|x|-w,0}     with w > 0

    When the input 'x' is an array, the output is computed element-wise.

     INPUTS
    ========
     x     - ND array
     gamma - positive, scalar or ND array with the same size as 'x'
     w     - positive, scalar or ND array with the same size as 'x'
    """
    
    
    gamma = None
    w     = None
    
    
    def __init__(self, gamma, w):

        if np.any( gamma <= 0 ):
            raise Exception("'gamma' must be positive")
            
        if np.any( w <= 0 ):
            raise Exception("'w' must be positive")
            
        if (not np.isscalar(gamma)) and (not np.isscalar(w)) and (np.shape(gamma) != np.shape(w)):
            raise Exception("'gamma' and 'w' must be either scalar-scalar, scalar-vector, or vector-vector with the same size")
            
        self.gamma = gamma;
        self.w     = w;


    
    def prox(self, x, tau):
        
        self.__check(x)
        
        scale = tau * self.gamma;
        
        # preliminaries
        abs_x  = np.abs(x);
        sign_x = np.sign(x);

        # 2nd branch
        p = self.w * sign_x;

        # 1st branch
        mask = abs_x < self.w;
        p[mask] = x[mask];

        # 3rd branch
        mask = abs_x > self.w + scale;
        p3 = x - scale * sign_x;
        p[mask] = p3[mask];
            
        return p;
        
        
        
    def fun(self, x):
        
        self.__check(x)
        
        return np.sum( self.gamma * np.maximum(0, np.abs(x) - self.w) );
        
        
        
    def __check(self, x):
            
        if (not np.isscalar(self.gamma)) and (np.shape(self.gamma) != np.shape(x)):
            raise Exception("'gamma' must be either scalar or the same size as 'x'")

        if (not np.isscalar(self.w)) and (np.shape(self.w) != np.shape(x)):
            raise Exception("'w' must be either scalar or the same size as 'x'")