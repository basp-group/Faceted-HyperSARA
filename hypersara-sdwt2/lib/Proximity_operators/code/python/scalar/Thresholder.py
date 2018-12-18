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

class Thresholder :
    """
    This class computes the proximity operator of the function

                   / gamma * a * x   if x < 0
            f(x) = | 0               if x = 0         with a <= b
                   \ gamma * b * x   otherwise

    When the input 'x' is an array, the output is computed element-wise.

     INPUTS
    ========
     x     - ND array
     gamma - positive, scalar or ND array with the same size as 'x'
     a     - scalar or ND array with the same size as 'x'
     b     - scalar or ND array with the same size as 'x'
    """
    
    a     = None
    b     = None
    
    def __init__(self, gamma, a, b):

        if np.any( gamma <= 0 ):
            raise Exception("'gamma' must be positive")
            
        if (not np.isscalar(a)) and (not np.isscalar(b)) and (np.shape(a) != np.shape(b)):
            raise Exception("'a' and 'b' must be either scalar-scalar, scalar-vector, or vector-vector with the same size")
            
        if np.any(a > b):
            raise Exception("'a' must be less than or equal to 'b'")
            
        self.a = a * gamma;
        self.b = b * gamma;


    
    def prox(self, x, tau):
        
        self.__check(x)
            
        return np.minimum(0, x - self.a * tau) + np.maximum(0, x - self.b * tau);
        
        
        
    def fun(self, x):
        
        self.__check(x)
        
        return self.a * np.minimum(0,x) + self.b * np.maximum(0,x);
        
        
        
    def __check(self, x):
            
        if (not np.isscalar(self.a)) and (np.shape(self.a) != np.shape(x)):
            raise Exception("'a' must be either scalar or the same size as 'x'")

        if (not np.isscalar(self.b)) and (np.shape(self.b) != np.shape(x)):
            raise Exception("'b' must be either scalar or the same size as 'x'")