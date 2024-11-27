from .manipulate import compress_dimension

import numpy as np
import scipy.interpolate
from copy import deepcopy

def binodal_from_slab(fields, field_idx_for_interfaces=0, interface_scale=1.0, plot=False):
    ''' Tool to estimate phase coexistance from a slab-like simulation
        given fields estimate binodal compositions and volumes of each phase 

        WARNING: the fields used in this function should come from DensityMolecule not DensitySpecies!!!!
    '''

    # --------------------------------------------
    # first compress fields into 1d if applicable
    # --------------------------------------------
     
    # Confirm not square. Assumes all fields have same cell size
    field = fields[0]
    assert(field.is_orthorhombic), "phase_coexistence_slab expects field to be orthorhombic"
   
    # get longest dimension and compress if needed
    if field.dim != 1:
        assert(not field.is_cubic()), "cell cannot be square for slab, need a long axis"
        longest_dim = np.diag(field.h).argmax()
        dims_to_compress = [0,1,2] # start with all
        dims_to_compress.remove(longest_dim) # then remove longest
        fields_1d = compress_dimension(fields, dims_to_compress) # from manipulate
    else:
        fields_1d = fields

    # --------------------------------------------
    # now find different phases and their volumes
    # --------------------------------------------

    field = fields_1d[field_idx_for_interfaces]
    idx_of_phases, vol_of_phases = _slab_estimate_grid_index_and_volume_of_phases(field,interface_scale, plot=plot)

    # ---------------------------------------------------
    # calculate average density of each phase and return
    # ---------------------------------------------------
    
    state = {} 
    nmolecules = len(fields_1d)
    nphases = len(idx_of_phases)

    # check if slab is valid
    valid_slab = True
    min_phase_width = 3 # width in data points, this helps exclude frames that have hotspots
    for idx_of_phase in idx_of_phases:
        if sum(idx_of_phase) < min_phase_width:
            valid_slab = False
    if not valid_slab:
        raise RuntimeError("Error when calculating binodal from slab, one of the phases is too narrow")

    # set nu (the relative volume of each phase)    
    state['nu'] = []
    for vol in vol_of_phases:
      state['nu'].append(vol)

    # TODO: instead of hardcoding phi and C I should have the user pass an argument so that the raw densities can be output (if appropriate)
  
    # set C of each phase (use the total sum of all fields)
    state['C'] = [0] * nphases
    for iphase in range(nphases):
      idx = idx_of_phases[iphase]
      total_density = 0.0
      for imol in range(nmolecules):
        total_density += np.average(fields_1d[imol].data[idx].real) # FIXME: causing runtime warnings and Nan?
      state['C'][iphase] = total_density

    # initialize phi_i's
    for imol in range(nmolecules):
      state[f'phi_{imol}'] = [0] * nphases

    # set phi of each molecule in each phase
    for iphase in range(nphases):
      idx = idx_of_phases[iphase]
      sum_phi = 0
      for imol in range(nmolecules):
        # FIXME: this is where runtime warnings and NANs come from
        phi = np.average(fields_1d[imol].data[idx].real)
        sum_phi += phi
        state[f'phi_{imol}'][iphase] = phi # average over idx of each phase

      # normalize to phi sum to unity (due to rounding errors in averaging)
      for imol in range(nmolecules):
        state[f'phi_{imol}'][iphase] /= sum_phi
      
    return state 

def center_1D_density_fields(density_fields, step_size=1, plot=False):
    '''
    Center a list of 1D density fields such that each dense phase's center of 
    mass is located at the center of each field's coordinates.
    
    Args:
        density_fields: A list of 1D density fields to be centered.
        step_size: The step size used when estimating the first derivative.
                   For noiser data, a higher step size can help.
        plot: A bool which specifies whether to plot the centering.
            
    Returns:
        The list of centered 1D density fields.
    '''
    # Find the indices of the interface for each density field.
    
    centered_density_fields = []
    for density_field in density_fields:
        density_field_initial = deepcopy(density_field)

        # Try finding the interface for various rolls of the density field.

        interface_found = False
        npw = density_field.npw[0]
        shift = npw // 10
        while not interface_found:
            try:
                phase_data = find_phase_densities(density_field, 
                                                  step_size=step_size)
                interface_found = True
            except:  
                density_field.data.real = np.roll(density_field.data.real, 
                                                  shift=shift)
        left_ind = phase_data[4]
        right_ind = phase_data[5]

        # Center the density field such that the dense phase's center of mass
        # is at the center of the field coordinates.

        density_data = density_field.data.real
        dense_phase_indices = np.arange(len(density_data))[left_ind + 1 :\
                                                           right_ind]
        dense_phase_points = density_data[dense_phase_indices]
        com_ind = round(np.sum(dense_phase_points * dense_phase_indices) /\
                        np.sum(dense_phase_points))
        shift = density_data.size // 2 - com_ind
        centered_density_data = np.roll(density_data, shift)

        # Plot the found centering if applicable.

        if plot:
            density_coords = density_field.coords
            import matplotlib.pyplot as plt
            fig = plt.figure(dpi=150) 
            ax1, ax2 = fig.subplots(1, 2, sharey=True)
            ax1.plot(density_coords, density_field_initial.data.real, 
                     marker='.', color='tab:blue')
            ax2.plot(density_coords, centered_density_data, 
                     marker='.', color='tab:green')
            ax2.axvline(density_coords[com_ind + shift],
                        color='k', linestyle='--')
            ax1.set_xlabel('Position')
            ax1.set_ylabel('Density')
            ax1.set_title('Initial density field', fontsize='medium')
            ax2.set_xlabel('Position')
            ax2.set_title('Centered density field', fontsize='medium')
            plt.tight_layout()
            plt.show()

        # Save the centered density field.

        density_field.data.real = centered_density_data
        centered_density_fields.append(density_field)

    # Return the centered density fields.

    return centered_density_fields

def find_phase_densities(density_field, step_size=1, dense_trim=0,
                         dilute_trim=0, plot=False, 
                         output_filename=None):
    '''
    Find the dense and dilute phases of a 1D density field. This is 
    accomplished by estimating the first derivative and finding the maximum and
    minimum values, which correspond to the left and right interface 
    respectively.

    Args:
        density_field: The 1D density field to be centered.
        step_size: The step size used when estimating the first derivative.
                   For noiser data, a higher step size can help.
        dense_trim: The amount (in grid points) to reduce each side of the 
                    dense phase.
        dilute_trim: The amount (in grid points) to reduce each side of the 
                     dilute phase.
        plot: A bool which specifies whether to plot the found phases.
        output_filename: The filename to write phase data to.
            
    Returns:
        (As a tuple) the dense phase density, its standard deviation, the 
        dilute phase density, its standard deviation, the left interface index,
        and the right interface index.
    '''
    # Pull data and preallocate first derivative array.

    density_coords = density_field.coords
    density_data = density_field.data.real
    first_derivative = np.zeros(density_data.size - step_size)

    # Calculate the first derivative at each viable point.

    for i in range(density_data.size - step_size):
        delta_density = density_data[i + step_size] - density_data[i]
        delta_position = density_coords[i + step_size] - density_coords[i]
        first_derivative[i] = delta_density / delta_position

    # Get the dense phase points and dilute phase points.

    left_ind = first_derivative.argmax()  # Left interface.
    right_ind = first_derivative.argmin()  # Right interface.
    if left_ind > right_ind:
        raise RuntimeError("Dense phase is straddling the periodic boundary\nTry centering the density profile with fieldkit's center_1D_density_fields function")
    dense_inds = np.zeros(len(density_data), dtype=bool)
    dense_inds[left_ind + dense_trim + 1 : right_ind - dense_trim] = True
    dense_phase_points = density_data[dense_inds]
    dilute_inds = np.ones(len(density_data), dtype=bool)
    dilute_inds[left_ind - dilute_trim  + 1 : right_ind + dilute_trim] = False
    dilute_phase_points = density_data[dilute_inds]

    # Calculate the dense and dilute phase densities and standard deviations.

    dense_phase_density = np.average(dense_phase_points)
    dense_phase_stdev = np.std(dense_phase_points, ddof=1)
    dilute_phase_density = np.average(dilute_phase_points)
    dilute_phase_stdev = np.std(dilute_phase_points, ddof=1)

    # Plot the found phases if applicable.

    if plot:
        import matplotlib.pyplot as plt
        fig = plt.figure(dpi=150) 
        ax1, ax2 = fig.subplots(1, 2, sharey=True)
        ax1.plot(density_coords, density_data, marker='.', color='tab:blue')
        ax1.axvline(density_coords[left_ind], color='tab:cyan', linestyle='--')
        ax1.axvline(density_coords[right_ind],color='tab:cyan', linestyle='--')
        ax2.plot(density_coords[dense_inds], density_data[dense_inds], 
                 marker='.', color='tab:green')
        ax2.plot(density_coords[:left_ind - dilute_trim], 
                 density_data[:left_ind - dilute_trim], 
                 marker='.', color='tab:green')  # Left dilute phase.
        ax2.plot(density_coords[right_ind + dilute_trim + 1:], 
                 density_data[right_ind + dilute_trim + 1:], 
                 marker='.', color='tab:green')  # Right dilute phase.
        ax1.set_xlabel('Position')
        ax1.set_ylabel('Density')
        ax1.set_title('Density field with found interface', fontsize='medium')
        ax2.set_xlabel('Position')
        ax2.set_title('Dense and dilute phases', fontsize='medium')
        plt.tight_layout()
        plt.show()

    # Write the phase data (if applicable) and return it.

    phase_data = (dense_phase_density, dense_phase_stdev, dilute_phase_density,
                  dilute_phase_stdev, left_ind, right_ind)
    if output_filename:
        with open(output_filename, 'w') as f:
            f.write('# dense density | dense stdev | dilute density | dilute stdev | left interface index | right interface index\n')
            for data in phase_data:
                f.write(f'{data} ')
            f.write('\n')
    return phase_data

def _slab_estimate_grid_index_and_volume_of_phases(field, interface_scale=1.0, plot=False):
    # load useful values
    assert(field.dim == 1), "fields must be 1 dimensional"
    Lx = field.h[0,0] 
    M = field.npw[0]
    x = field.coords.real # only use real 
    y = field.data.real # only use real 
    
    # compute derivative using splines 
    tck = scipy.interpolate.splrep(x,y)
    dydx = scipy.interpolate.splev(x, tck, der=1)

    # Below I will follow the following convention when I say "left" and "right"
    # "left" is MAX of deriv (i.e. where density is increasing INTO slab)
    # "right" is MIN of deriv (i.e. where density is decreasing OUT OF slab)
    #
    #    ^ 
    #    |              _______                
    # y  |             /       \
    #    | phaseII    / phaseI  \  phaseII
    #    |           /           \             
    #    | ----------             ---------  
    #    |           ^           ^           
    #    |          LEFT        RIGHT        
    #     ------------------------------------>
    #                                         x
    #
    # Due to PBCs, sometimes "left" will not appear on the left. e.g
    #    ^ 
    #    |______                         ______
    # y  |      \                       /   
    #    |phaseI \      phaseII        / phaseI
    #    |        \                   /     
    #    |         -------------------      
    #    |       ^                      ^   
    #    |      RIGHT                  LEFT 
    #     ------------------------------------>
    #                                         x
    # 
    # this situation is accounted for below. These sketches are mostly in case I 
    # need to debug this function in the future...

    
    # find interface using the max and min of 1st deriv
    idx_interface_left = dydx.argmax()  # note: "left" convention above
    idx_interface_right = dydx.argmin() # note: "right" convention above
    
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
      
      # technically I could use HALF the interface width here, but using the entire interface width
      # on both sides of the interface is more conservative
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
    if plot:
      import matplotlib.pyplot as plt
      import matplotlib
      #matplotlib.use('Agg') # use different backend to avoid picotte issues
      fig = plt.figure()
      ax = fig.subplots(1,1)
        
      ax.plot(x,y,c='grey') 
      ax.plot(x[idx_center_phaseI],y[idx_center_phaseI],label='center of I',marker='*',ls='',ms=10,c='tab:blue') 
      ax.plot(x[idx_phaseI],y[idx_phaseI],marker='.',ls='',label='phase I',c='tab:blue') 

      ax.plot(x[idx_center_phaseII],y[idx_center_phaseII],label='center of II',marker='*',ls='',ms=10,c='tab:red') 
      ax.plot(x[idx_phaseII],y[idx_phaseII],marker='.',ls='',label='II',c='tab:red') 
      ax.plot(x,dydx,marker='.', label='derivative', color='tab:purple') 
      ylim = ax.get_ylim()
      ax.plot([x[idx_interface_left]]*2, ylim, ls='--',c='k')
      ax.plot([x[idx_interface_right]]*2, ylim, ls='--',c='k')
      ax.legend()
      plt.show()
    return idx_of_phases, vol_of_phases
   
