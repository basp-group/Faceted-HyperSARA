"""
 Version : 1.0 (27-04-2017)
 
 Author  : Giovanni Chierchia

 Copyright (C) 2017
 
 This file is part of the codes provided at http://proximity-operator.net
 
 By downloading and/or using any of these files, you implicitly agree to 
 all the terms of the license CeCill-B (available online).
 """
   
# -*- coding: utf-8 -*-

from Thresholder import Thresholder

class HingeLoss(Thresholder) :
    """
    This class computes the proximity operator of the function

                      f(x) = gamma * max{x,0}

    When the input 'x' is an array, the output is computed element-wise.

     INPUTS
    ========
     x     - ND array
     gamma - positive, scalar or ND array with the same size as 'x'
    """
        
    def __init__(self, gamma):

        Thresholder.__init__(self, gamma, 0, 1)