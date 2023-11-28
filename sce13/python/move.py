# -*- coding: utf-8 -*-
"""
Created on Sat Jul  9 14:44:19 2022

@author: Jungle
"""
import numpy as np
def smoothmoveavg(init_list,N,string):
    '''coculate the smoothmoveavg of list
       init_list: the initial list
       N: move window
       M: length of list
       string: mode of avg
    '''
    n = np.ones(N)
    weights = n/N
    new_list = np.convolve(weights, init_list,mode = string)
    return new_list



