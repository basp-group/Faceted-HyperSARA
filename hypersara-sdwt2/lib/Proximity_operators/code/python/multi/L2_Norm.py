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

class L2_Norm :
    """
    This class computes the proximity operator of the function

                        f(x) = gamma * ||x||_2

    When the input 'x' is an array, the computation can vary as follows:
     - direction <  0 --> 'x' is processed as a single vector [DEFAULT]
     - direction => 0 --> 'x' is processed block-wise along the specified direction
                          (0 -> columns, 1 -> rows, etc.)

     INPUTS
    ========
     x         - ND array
     gamma     - positive, scalar or ND array compatible with the blocks of 'x'
     direction - integer, direction of block-wise processing
    """
    
    gamma     = None
    direction = None
    
    def __init__(self, gamma, direction=None):

        if np.any( gamma <= 0 ):
            raise Exception("'gamma' must be positive")
        
        if (direction is not None) and (direction < 0):
            direction = None
            
        self.gamma     = gamma;
        self.direction = direction;


    
    def prox(self, x, tau):
        
        self.__check(x)
        
        xx = np.sqrt(np.sum( np.square(x), self.direction, keepdims=True ))
        xx = np.maximum(0, 1 - tau*self.gamma / xx)
        p  = x * xx
        
        return p
        
        
        
    def fun(self, x):
        
        self.__check(x)
        
        xx = np.sqrt(np.sum( np.square(x), self.direction, keepdims=True ))
        p  = np.sum( self.gamma * xx )        
        
        return p
        
        
        
    def __check(self, x):
            
        sz = np.shape(x)
        if self.direction is not None:
            sz = list(sz)
            sz[self.direction] = 1
            sz = tuple(sz)
        
        if (not np.isscalar(self.gamma)) and ( (self.direction is None) or np.shape(self.gamma) != sz ):
            raise Exception("'gamma' must be either scalar or compatible with the blocks of 'x'")