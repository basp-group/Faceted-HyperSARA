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

class BoxConstraint :
    """
    This class computes the projection onto the constraint set:

                           low <= x <= high

    When the input 'x' is an array, the output is computed element-wise.

     INPUTS
    ========
     x    - ND array
     low  - scalar or ND array with the same size as 'x'
     high - scalar or ND array with the same size as 'x'
     """ 
    
    low  = None
    high = None
    
    def __init__(self, low, high):

        if np.any( low >= high ):
            raise Exception("'low' must be lower than 'high'")
            
        self.low  = low;
        self.high = high;


    
    def prox(self, x, tau=None):
        
        self.__check(x)
            
        return np.maximum( self.low, np.minimum(x, self.high) )
        
        
        
    def fun(self, x):
        
        self.__check(x)
        
        if np.all(self.low <= x) and np.all(x <= self.high):
            return 0
        else:
            return np.inf
        
        
        
    def __check(self, x):
            
        if (not np.isscalar(self.low)) and (np.shape(self.low) != np.shape(x)):
            raise Exception("'low' must be either scalar or the same size as 'x'")

        if (not np.isscalar(self.high)) and (np.shape(self.high) != np.shape(x)):
            raise Exception("'high' must be either scalar or the same size as 'x'")