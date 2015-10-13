'''
@author: Moritz Augustin
'''

import numba
import numpy
import mwmc

@numba.jit
def alpha_count_center(data_win):
    """ counts positive data values over alle species of (only) the
        center cell
    """
    ws_h = (data_win.shape[0]-1)/2 # ws_h: winsize half: (winsize-1)/2
    offset_i = ws_h
    offset_j = ws_h
    
    sum_ = 0.0
    for k in range(data_win.shape[2]): # loop over all species
        if data_win[offset_i, offset_j, k] > 0: # center cell, species k
            sum_ += 1
            
    return sum_

@numba.jit
def alpha_shannon_center_percellnormalized(data_win):
    """ computes the shannon measure -sum_k p_k*ln(p_k) (where k is 
        the species index) of (only) the center cell
        the p_k are not the probabilities from the max. entropy method 
        but are normalized such that sum_k p_k = 1 in each cell
    """
    ws_h = (data_win.shape[0]-1)/2 # ws_h: winsize half: (winsize-1)/2
    offset_i = ws_h
    offset_j = ws_h
    
    center_species = data_win[offset_i, offset_j, :]
    # neglect nodata for normalization
    normalization = numpy.sum(center_species[center_species>0])
    
    sum_ = 0.0
    if normalization > 0: # if scale is 0 there is nothing to do...
        for k in range(data_win.shape[2]): # loop over all species
            center_k = center_species[k]/normalization
            # the lim p->0 of p*log(p) = 0 (log(0) does not exist)
            # and negative entries corresponding to NODATA values
            # are neglected
            if center_k > 0: # center cell, species k
                sum_ += -center_k*numpy.log(center_k)
    
    if sum_ > 0:
#         alpha = numpy.exp(sum_)
        return sum_
    else:  # not a single species occurred in that cell
        return mwmc.nodata_global

@numba.jit
def alpha_shannon_center_unnormalized(data_win):
    """ computes the shannon measure -sum_k p_k*ln(p_k) (where k is 
        the species index) of (only) the center cell
        the p_k are the probabilities from the max. entropy method 
        i.e. they are NOT normalized per cell 
    """
    ws_h = (data_win.shape[0]-1)/2 # ws_h: winsize half: (winsize-1)/2
    offset_i = ws_h
    offset_j = ws_h
    
    center_species = data_win[offset_i, offset_j, :]
    
    sum_ = 0.0
    for k in range(data_win.shape[2]): # loop over all species
        center_k = center_species[k]
        # the lim p->0 of p*log(p) = 0 (log(0) does not exist)
        # and negative entries corresponding to NODATA values
        # are neglected
        if center_k > 0: # center cell, species k
            sum_ += -center_k*numpy.log(center_k)
    
    if sum_ > 0:
#         alpha = numpy.exp(sum_)
        return sum_
    else:  # not a single species occurred in that cell
        return mwmc.nodata_global

@numba.jit
def alpha_simpson_center_percellnormalized(data_win):
    """ computes the simspon concentration measure 1/(sum_k p_k*p_k) (where k is 
        the species index) of (only) the center cell
        the p_k are not the probabilities from the max. entropy method 
        but are normalized such that sum_k p_k = 1 in each cell
    """
    ws_h = (data_win.shape[0]-1)/2 # ws_h: winsize half: (winsize-1)/2
    offset_i = ws_h
    offset_j = ws_h
    
    center_species = data_win[offset_i, offset_j, :]
    # neglect nodata for normalization
    normalization = numpy.sum(center_species[center_species>0])
    
    sum_ = 0.0
    if normalization > 0: # if scale is 0 there is nothing to do...
        for k in range(data_win.shape[2]): # loop over all species
            center_k = center_species[k]/normalization
            # negative entries corresponding to NODATA values
            # are neglected
            if center_k > 0: # center cell, species k
                sum_ += center_k*center_k
    if sum_ > 0:
        return 1/sum_
    else:  # not a single species occurred in that cell
        return mwmc.nodata_global

# dict for function selection by cmd line argument
measures_dict = {'alpha_count_center': 
                    alpha_count_center,
                 'alpha_shannon_center_percellnormalized': 
                    alpha_shannon_center_percellnormalized,
                 'alpha_shannon_center_unnormalized':
                    alpha_shannon_center_unnormalized,
                 'alpha_simpson_center_percellnormalized': 
                    alpha_simpson_center_percellnormalized}