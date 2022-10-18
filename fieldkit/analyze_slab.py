from .manipulate import compress_dimension

import numpy as np
import scipy.interpolate

def binodal_from_slab(fields, field_idx_for_interfaces=1, interface_scale=1.0):
    ''' Tool to estimate phase coexistance from a slab-like simulation
        given fields estimate binodal compositions and volumes of each phase 

        WARNING: the fields used in this function should come from DensityMolecule not DensitySpecies!!!!

    '''

    # --------------------------------- 
    # first compress fields into 1d
    # --------------------------------- 
     
    # Confirm not square. Assumes all fields have same cell size
    field = fields[0]
    assert(field.is_orthorhombic), "phase_coexistence_slab expects field to be orthorhombic"
    assert(not field.is_cubic()), "cell cannot be square for slab, need a long axis"
   
    # get longest dimension 
    longest_dim = np.diag(field.h).argmax()
    dims_to_compress = [0,1,2] #start with all
    dims_to_compress.remove(longest_dim) # then remove longest
   
    # compress to 1d 
    fields_1d = compress_dimension(fields, dims_to_compress) # from manipulate

    # ----------------------------------------------
    # now find different phases and their volumes
    # ----------------------------------------------
    field = fields_1d[field_idx_for_interfaces]
    idx_of_phases, vol_of_phases = _slab_estimate_grid_index_and_volume_of_phases(field,interface_scale)


    # ---------------------------------------------------
    # calculate average density of each phase and return
    # ---------------------------------------------------

    state = {} 
    nmolecules = len(fields_1d)
    nphases = len(idx_of_phases)

    # set nu (the relative volume of each phase)    
    state['nu'] = []
    for vol in vol_of_phases:
      state['nu'].append(vol)

    # initialize phi_i's
    for imol in range(nmolecules):
      state[f'phi_{imol}'] = [0] * nphases

    # set phi of each molecule in each phase
    for iphase in range(len(idx_of_phases)):
      idx = idx_of_phases[iphase]
      sum_phi = 0
      for imol in range(nmolecules):
        phi = np.average(fields_1d[imol].data[idx].real)
        sum_phi += phi
        state[f'phi_{imol}'][iphase] = phi # average over idx of each phase

      # normalize to phi sum to unity (due to rounding errors in averaging)
      for imol in range(nmolecules):
        state[f'phi_{imol}'][iphase] /= sum_phi
      
    return state 

def _slab_estimate_grid_index_and_volume_of_phases(field, interface_scale=1.0):
    # load useful values
    assert(field.dim == 1), "fields must be 1 dimensional"
    Lx = field.h[0,0] 
    M = field.npw[0]
    x = field.coords
    y = field.data
    
    # compute derivative using splines 
    tck = scipy.interpolate.splrep(x,y)
    dydx = scipy.interpolate.splev(x, tck, der=1)

    # find interface using the max and min of 1st deriv
    idx_interface_left = dydx.argmax()
    idx_interface_right = dydx.argmin()
    
    # calculate center of phaseI (used to normalize interface widths)
    # shift right idx by PBC until it is actually on the right
    tmp_idx_interface_right = idx_interface_right
    while (tmp_idx_interface_right < idx_interface_left): tmp_idx_interface_right += M
    idx_center_phaseI = int(0.5*(idx_interface_left + tmp_idx_interface_right))
    while (idx_center_phaseI <  0): idx_center_phaseI += M
    while (idx_center_phaseI >= M): idx_center_phaseI -= M
    volI = (tmp_idx_interface_right - idx_interface_left) / M

    # calculate center of phaseII (used to normalize interface widths)
    # shift right idx by PBC until it is actually on the LEFT
    tmp_idx_interface_right = idx_interface_right
    while (tmp_idx_interface_right > idx_interface_left): tmp_idx_interface_right -= M
    idx_center_phaseII = int(0.5*(idx_interface_left + tmp_idx_interface_right))
    while (idx_center_phaseI <  0): idx_center_phaseI += M
    while (idx_center_phaseI >= M): idx_center_phaseI -= M
    volII = (idx_interface_left - tmp_idx_interface_right) / M

    # find width of interface
    # note: using center of phaseI/II for normalization
    deltay = y[idx_center_phaseI] - y[idx_center_phaseII]
    width_interface_left = deltay / abs(dydx[idx_interface_left]) * interface_scale
    width_interface_right = deltay / abs(dydx[idx_interface_right]) * interface_scale

    # storage to determine whether 
    idx_phaseI = np.zeros(M,dtype=bool)
    idx_phaseII = np.zeros(M,dtype=bool)
    if idx_interface_left < idx_interface_right:
      is_phaseI = False
    else:
      is_phaseI = True

    # now assign grid points to phaseI/II    
    for i in range(M):
      
      dx_left = x[i] - x[idx_interface_left] # how far past left interface
      while (dx_left <  0.5*Lx): dx_left += Lx # apply min image convention
      while (dx_left >= 0.5*Lx): dx_left -= Lx
      dx_left = abs(dx_left) # then take abs value (this makes it invariant to phaseI/II)

      dx_right = x[idx_interface_right] - x[i] # how far to right interface
      while (dx_right <  0.5*Lx): dx_right += Lx # apply min image convention
      while (dx_right >= 0.5*Lx): dx_right -= Lx
      dx_right = abs(dx_right) # then take abs value (this makes it invariant to phaseI/II)

      if dx_right > width_interface_right and dx_left > width_interface_left :
        if is_phaseI:
          idx_phaseI[i] = True
        else:
          idx_phaseII[i] = True
      if i == idx_interface_left:
        is_phaseI = True
      elif i == idx_interface_right:
        is_phaseI = False
    
    # package data for output
    idx_of_phases = [idx_phaseI, idx_phaseII]
    vol_of_phases = [volI, volII] # normalized volumes (they should sum to one)
   
    # plotting  -- for debugging
    #import matplotlib.pyplot as plt
    #plt.plot(x,y) 
    #plt.plot(x[idx_center_phaseI],y[idx_center_phaseI],label='centerI',marker='*',ms=10) 
    #plt.plot(x[idx_center_phaseII],y[idx_center_phaseII],label='centerII',marker='*',ms=10) 
    #plt.plot(x[idx_phaseI],y[idx_phaseI],marker='.',ls='',label='I') 
    #plt.plot(x[idx_phaseII],y[idx_phaseII],marker='.',ls='',label='II') 
    #plt.plot(x,dydx,marker='.', label='derivative') 
    #plt.legend()
    #plt.show()
    return idx_of_phases, vol_of_phases
   
     
