#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

from .field import *
from .field_io import *
    
def update_species_field_gaussians(field, centers, npw, magnitude, x, y, z):
    #parameters for all spheres phases
    gaussian_width = 0.1
    inv_width2 = 1.0 / gaussian_width / gaussian_width 
        
    for c in centers:
        for ix,myx in enumerate(x):
            dx = (myx - c[0])
            #apply pbc
            if dx > 0.5:   dx -= 1.0
            if dx < -0.5:  dx += 1.0
            
            for iy,myy in enumerate(y):
                dy = (myy - c[1])
                if dy > 0.5:   dy -= 1.0
                if dy < -0.5:  dy += 1.0
            
                for iz,myz in enumerate(z):
                    dz = (myz - c[2])
                    if dz > 0.5:   dz -= 1.0
                    if dz < -0.5:  dz += 1.0

                    dr2 = dx*dx + dy*dy + dz*dz

                    field[ix,iy,iz] -= np.exp(-0.5*dr2*inv_width2)
                    
    field = (field - np.min(field)) / (np.max(field) - np.min(field)) * magnitude
    field -= np.average(field)  
    
    return field


def update_species_field_levelset(phase, field, magnitude, originshift, x, y, z):
    pi = np.pi

    xx, yy, zz = np.meshgrid(x, y, z)
    
    x = xx + originshift[0] #0.125
    y = yy + originshift[1] # 0.125
    z = zz + originshift[2] # 0.125

    if phase == 'double-gyroid':
        field += 0.8*(np.sin(4*pi*x)*np.sin(2*pi*z)*np.cos(2*pi*y) \
                 + np.sin(4*pi*y)*np.sin(2*pi*x)*np.cos(2*pi*z) \
                 + np.sin(4*pi*z)*np.sin(2*pi*y)*np.cos(2*pi*x))\
                 - 0.2*(np.cos(4*pi*x)*np.cos(4*pi*y) \
                   + np.cos(4*pi*y)*np.cos(4*pi*z) \
                   + np.cos(4*pi*z)*np.cos(4*pi*x))
    elif phase == 'single-diamond' or phase == 'alt-diamond':
    # Fig 3a, diamond (Eq 13)
      field -= np.cos(2*pi*x)*np.cos(2*pi*y)*np.cos(2*pi*z) \
          + np.sin(2*pi*x)*np.sin(2*pi*y)*np.cos(2*pi*z)\
          + np.sin(2*pi*x)*np.cos(2*pi*y)*np.sin(2*pi*z)\
          + np.cos(2*pi*x)*np.sin(2*pi*y)*np.sin(2*pi*z)  
    else:
     raise RuntimeError(f"Invalid phase {phase}")    
     
    field = (field - np.min(field)) / (np.max(field) - np.min(field)) * magnitude
    field -= np.average(field)  
    return field

def initialize_sphere_positions(phase, field, npw, magnitude): 
    ndim = 3
    #create centers
    centersA = []
    centersC = []
    if phase == "A15":
        # A15, space group 223 (Pm-3m)  
        # wyckoff positions 2a and 6c
        centersA = [[0,0,0],[0.5,0.5,0.5],  # position 2a: 0,0,0
        [0.25,0,0.5], [0.75,0,0.5], [0.5,0.25,0],[0.5,0.75,0],[0,0.5,0.25],[0,0.5,0.75] # position 6c: 0.25, 0, 0.5
        ]
        
    elif phase == "C15" or "alt-C15": #creates centersA 
        # set centersA, centerC using C15 positions
        general_positions = [[0,0,0],[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]] # C15, space group 227 (Fd-3m)  
        special_positions_8a = [[0,0,0],[0.75,0.25,0.75]] # Wyckoff position 8a
        special_positions_16d = [[5./8,5./8,5./8],[3./8,7./8,1./8],[7./8,1./8,3./8],[1./8,3./8,7./8]]
            
        for gp in general_positions:
          for sp in special_positions_8a:
            pos = np.add(gp,sp)
            for i in range(ndim):
              while (pos[i] >= 1): pos[i] -= 1
              while (pos[i] <  0): pos[i] += 1
            centersA.append(list(pos))
    
          for sp in special_positions_16d:
            pos = np.add(gp,sp)
            for i in range(ndim):
              while (pos[i] >= 1): pos[i] -= 1
              while (pos[i] <  0): pos[i] += 1
            
            if phase == "alt-C15":
                centersC.append(list(pos))
            else:
                centersA.append(list(pos))
        
    npw=tuple(npw)
    x = np.linspace(0,1,npw[0])
    y = np.linspace(0,1,npw[1])
    z = np.linspace(0,1,npw[2])
        
    if phase == "alt-C15":
        field = [update_species_field_gaussians(field, centersA, npw, magnitude, x, y, z), update_species_field_gaussians(field, centersC, npw, magnitude, x, y, z)]
    else:
        field = update_species_field_gaussians(field, centersA, npw, magnitude, x, y, z)
    return field
    
     
def initialize_phase(phase, npw, h):    
    """ Creates a list of Field objects for 3D phases with 2 or 3 species. 
    
    Adapted from FTS-tools/field/init_fields_3d_2spec.py and from FTS-tools/fields/ini_fields_3d_3spec.py.
    2 species 3D phases implemented in this function: A15. C15, double-gyroid, single-diamond
    3 species 3D phases implemented in this function: alt- diamond, alt-C15.
    
    Args:
        phase: name of the phase (see description above for phases implemented in this function)
        npw: number of grid points for each dimension
        h: a cell tensor, in which the columns are box vectors that describe the coordinate system
    
    Returns:
        field_list: list of Field objects, each object for each species.
        
    Raises:
        NotImplementedError: if phase is not a phases implemented in this function (see description above), the function will exit.
    """
    
    sphere_phases = ["A15", "C15", "alt-C15"]
    levelset_phases = ["double-gyroid", "single-diamond", "alt-diamond"]
    
    #parameters for all phases:
    npw = npw
    magnitude = 10
    
    x = np.linspace(0,1,npw[0])
    y = np.linspace(0,1,npw[1])
    z = np.linspace(0,1,npw[2])
    xx, yy, zz = np.meshgrid(x, y, z)
    
    #fields
    w_A = np.zeros(npw)
    w_B = np.zeros(npw)
    w_C = np.zeros(npw)

    if phase in sphere_phases:
        if phase == "alt-C15":
            nspecies = 3
            w_A =  initialize_sphere_positions(phase, w_A, npw, magnitude)[0]
            w_C =  initialize_sphere_positions(phase, w_C, npw, magnitude)[1]

        else:    
            nspecies = 2
            w_A = initialize_sphere_positions(phase, w_A, npw, magnitude)
        
    elif phase in levelset_phases:
        nspecies = 2
        originshift = (0, 0, 0)
        if phase == "double-gyroid":
            originshift = (0.125, 0.125, 0.125)
        w_A = update_species_field_levelset(phase, w_A, magnitude, originshift, x, y, z)

        if phase == "alt-diamond":
            nspecies = 3
            w_A = update_species_field_levelset(phase, w_A, magnitude, (0.125,0.125,0.125), x, y, z)
            w_C = update_species_field_levelset(phase, w_C, magnitude, (0.625,0.625,0.625), x, y, z)
    else:
        raise NotImplementedError("Invalid phase: {}".format(phase))
    
    if nspecies == 2:
        field_list = [Field(npw_Nd=npw, data=w_A, h = h), Field(npw_Nd=npw, data=w_B, h = h)]
    elif nspecies == 3:
        w_B = np.random.random(npw)
        field_list = [Field(npw_Nd=npw, data=w_A, h=h), Field(npw_Nd=npw, data=w_B, h=h), Field(npw_Nd=npw, data=w_C, h=h)]
    
    return field_list
        
    
    
    
def add_ellipse(field, center, axis_lengths, smoothing_width = 0.05,height=1):
    """ adds an elipse
    
    Args:
        field: field to add ellipse to (edited in place)
        center: center of ellipse (in scaled [0-1])
        axis_lengths: lengths of ellipse axis in each dimension (in scaled [0-1])
        smoothing_width: sigma to use for smoothing gaussian filter
        height : height of ellipse
     
    Returns:
        field: field edited in place
        
    Raises:
        None currently.

    """
    
  
    npw = field.npw_Nd
    dim = len(npw)
    assert(len(axis_lengths) == dim)
    dx = np.zeros(dim)
    data = np.zeros(npw)

    for index in np.ndindex(npw):
        scaled_index = index / np.array(npw)
        # Equation of ellipse: x**2/a**2 + y**2/b**2 + z**2/c**2  = 1
        dist2 = 0
        for idim in range(dim):
          dx[idim] = scaled_index[idim] - center[idim]
          if dx[idim] > 0.5:   dx[idim] -= 1.0
          if dx[idim] < -0.5:  dx[idim] += 1.0
          dist2 += dx[idim]*dx[idim] / axis_lengths[idim]**2
  
        if dist2 < 1: # is inside ellipse
          data[index] += height
        else: # outside ellipse 
          pass
    
    # now smear using gaussian filter
    import scipy.ndimage 
    data_smooth=np.zeros(npw)
    sigma = [int(i*smoothing_width) for i in npw] 
    scipy.ndimage.gaussian_filter(data,sigma=sigma,output=data_smooth, mode='wrap')

    # return in place 
    field.data += data_smooth
     
        
        
        
