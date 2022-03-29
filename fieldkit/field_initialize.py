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
        
    
    
    
def add_ellipse(field, center, axis_lengths, smoothing_width = 0.05,height=1,units="scaled"):
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
    h = field.h
    assert(len(axis_lengths) == dim)
    dx = np.zeros(dim)
    data = np.zeros(npw)
    
    if units == "unscaled":
      assert(field.is_orthorhombic()), f"Error: h must be orthorhombic in add_gaussian function. {h}"
      boxl = np.diag(h)
      boxh = 0.5*boxl
    elif units == "scaled":
      boxl = [1.0]*dim
      boxh = [0.5]*dim


    for index in np.ndindex(npw):
        scaled_index = index / np.array(npw)
        if units == "unscaled":
            rgrid = np.dot(h, scaled_index)
        elif units == "scaled":
            rgrid = scaled_index
        # Equation of ellipse: x**2/a**2 + y**2/b**2 + z**2/c**2  = 1
        dist2 = 0
        for idim in range(dim):
          dx[idim] = rgrid[idim] - center[idim]
          while dx[idim] >  boxh[idim]: dx[idim] -= boxl[idim]
          while dx[idim] < -boxh[idim]: dx[idim] += boxl[idim]
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
     
def add_gaussian(field, center, sigma, height=1):
    ''' add a gaussian to a field. 

    Args: 
        field: field to add ellipse to (edited in place)
        units: units to use for args `center` and `sigma`. 
               `units='real'` is default currently and `center` must fall within `h` cell tensor of field
        center: center of Gaussian (in real units)
        height: height of Gaussian
        
    '''
    npw = field.npw_Nd
    invM = 1.0/np.prod(npw)
    h = field.h # cell tensor
    dim = field.dim
    dx = np.zeros(dim)
    a = sigma # warning: the sigma used for smearing in FTS contians an extra factor of sqrt(2)
    rcut = 5.0*sigma
    rcut2 = rcut*rcut
    center = np.array(center) # recast center to array
   
    # box must be orthorhombic
    #assert(field.is_orthorhombic()), f"Error: h must be orthorhombic in add_gaussian function. {h}"
    
    boxl = np.diag(h)
    boxh = 0.5*boxl

    # check that center is within box
    # option 1: assert that center myst be within box
    #assert (np.all(center >= 0.0)), f"center must be within box {center}, {h}"
    #assert (np.all(center < boxl)), f"center must be within box {center}, {h}"

    # option 2: apply PBC to center
    for idim in range(dim):
        while center[idim] <  0.0:        center[idim] += boxl[idim]
        while center[idim] >= boxl[idim]: center[idim] -= boxl[idim]

    data = np.zeros(npw)
    # TODO: its very inefficient to manually loop over all grid points
    # should create a neighborlist of the relevant grid points nearby "center"
    
    # original (likely slow) version
    for index in np.ndindex(npw):
        scaled_index = index / np.array(npw)
        rgrid = np.dot(h, scaled_index)

        # Equation of gaussian: Gamma(r) = 1 / (2*pi)**1.5 / a**3 * exp(-r2 / 2 / a**2)
        dr2 = 0
        for idim in range(dim):
          dx[idim] = rgrid[idim] - center[idim]
          while dx[idim] >  boxh[idim]: dx[idim] -= boxl[idim]
          while dx[idim] < -boxh[idim]: dx[idim] += boxl[idim]
          dr2 += dx[idim] * dx[idim]

        # use cutoff to marginally speed up calculation
        if dr2 < rcut2:
            gamma = invM * 1.0 / (np.sqrt(2.0*np.pi)*a)**dim * np.exp(-0.5*dr2 / a**2) # note normalization depends on dim
            data[index] += gamma
    
    # scale so that sum of data == height. This helps fix discritization errors if the full Gaussian isn't fully resolved on grid
    data *= height / np.sum(data)
    #print(f"{np.sum(data) = }") # should be ~1
    
    # now add to field
    field.data += data


def particle_to_field(trjfile, topfile, frames_to_average, npw, P):
    ''' Initialize a field using a the particle coordinates and Hockney/Eastwood function

    Args: 
        trjfile: trajectory file. Any file format compatable with mdtraj should work.
        topfile: topology file. Any file format compatable with mdtraj should work.
        frames_to_average: frame indicies to average over when converting to field. list, or int (if single frame)
        npw: dimension of grid for output field (note that the voxels must be approx. cubic to use current Hockney-Eastwood mapping function
        P: order of Hockney-Eastwood assignment function (see Deserno and Holm 1998 for more details)
    '''

    # open trajectory using MD traj
    import mdtraj as md
    #t = md.load_lammpstrj(trjfile,top=topfile)
    t = md.load(trjfile,top=topfile)

    # if only a single frame is specified, turn it into a list
    if type(frames_to_average) == int:
      frames_to_average = [frames_to_average]
    nframes_to_average = len(frames_to_average)
  
    # extract useful parameters from first frame (e.g. atomtypes, cell size, etc)
    # WARNING: assumes that all atom types are present in 1st frame
    frame = t[frames_to_average[0]]
    table, bonds = frame.topology.to_dataframe()
    #key = 'resSeq'
    key = 'name' # TODO would be nice is key wasn't hardcoded. Ideally the code would try several options
    atomtypes_list = []
    for atomtype in table[key]:
        if not atomtype in atomtypes_list:
            atomtypes_list.append(atomtype)
    natomtypes = len(atomtypes_list)
    boxl = frame.unitcell_lengths[0]
    h = np.diag(boxl) # cell size of first frame

    # perform some checks over the frames
    for i in range(1,nframes_to_average):
      frame_index = frames_to_average[i]
      frame = t[frame_index] 
      boxl_tmp = frame.unitcell_lengths[0]
      h_tmp = np.diag(boxl_tmp)
      assert(np.all(h_tmp == h)),f"cell size must be constant across all frames used in averaging {h} {h_tmp}"
      assert(np.all(frame.unitcell_angles == 90.0)), "box must be orthorhombic"

    # using box size, initialize new fields (one for each type of atom)
    print(f"Creating {natomtypes} fields for {natomtypes} found atom types")
    fields = []
    for itype in range(natomtypes):
      fields.append(Field(npw_Nd = npw, h=h))

    for frame_index in frames_to_average:
      print(f"Processing frame {frame_index} containing {frame.n_atoms} atoms")
      # find specified frame
      frame = t[frame_index] 

      # loop through all particles in frame and call add_gaussian function
      for iatom in range(frame.n_atoms):
          myP = P # all atom types currently use same sigma. Should be able to generalize...
          pos = frame.xyz[0][iatom]
          atomtype = table[key][iatom] # TODO: consider using something other than resSeq
          atomtype_index = atomtypes_list.index(atomtype)
          
          # TODO: calling this function within a double for-loop over particles and frames is killing performance
          #       it would be worth thinking about a more efficient implementation
          # initial idea: try pybinding this function or numba?
          add_hockney_eastwood_function(fields[atomtype_index], center=pos, P=myP, height=1)
          
          # TODO: in principle it should be possible to also initialize using Gaussians
          # however add_gaussians function currently searches over all grid points which makes it much to slow
          #add_gaussian(field, center=pos, sigma=mysigma, height=1)
    
    # fields contain total over ALL frames, need to compute average
    for field in fields:
      field.data /= nframes_to_average

    # output (for debugging)
    #write_to_file("fields.dat",fields)
    #write_to_VTK("fields.vtk",fields)

    # return field
    return fields

def particle_to_field_gsd(gsdfile, frames_to_average, npw, P, normalize=False):
    ''' Initialize a field using a the particle coordinates and Hockney/Eastwood function

    Args: 
        gsdfile: trajectory file in gsd format 
        frames_to_average: frame indicies to average over when converting to field. list, or int (if single frame)
        npw: dimension of grid for output field (note that the voxels must be approx. cubic to use current Hockney-Eastwood mapping function
        P: order of Hockney-Eastwood assignment function (see Deserno and Holm 1998 for more details)
        normalize: whether or not to normalize densities by rho0
    '''

    # open trajectory using gsd
    import gsd.hoomd
    t = gsd.hoomd.open(gsdfile,'rb')
    
    # if only a single frame is specified, turn it into a list
    if type(frames_to_average) == int:
      frames_to_average = [frames_to_average]
    nframes_to_average = len(frames_to_average)
  
    # extract useful parameters from first frame (e.g. atomtypes, cell size, etc)
    # WARNING: assumes that all atom types are present in 1st frame
    frame = t[frames_to_average[0]]

    atomtypes_list = frame.particles.types
    natomtypes = len(atomtypes_list)
    boxl = frame.configuration.box[0:3]
    h = np.diag(boxl) # cell size of first frame
    assert(np.all(frame.configuration.box[3:] == 0)), "Cell must be orthorhombic"
    dim = frame.configuration.dimensions 
    natoms = frame.particles.N
    V = np.linalg.det(h)
    rho0 = natoms / V
    
    # perform some checks over the frames
    for i in range(1,nframes_to_average):
      frame_index = frames_to_average[i]
      frame = t[frame_index] 
      natoms_tmp = frame.particles.N
      boxl_tmp = frame.configuration.box[0:3]
      h_tmp = np.diag(boxl_tmp)
      assert(natoms_tmp == natoms),f"natoms must be constant across all frames used in averaging {natoms} {natoms_tmp}"
      assert(np.all(h_tmp == h)),f"cell size must be constant across all frames used in averaging {h} {h_tmp}"
      assert(np.all(frame.configuration.box[3:] == 0.0)), "box must be orthorhombic"

    # using box size, initialize new fields (one for each type of atom)
    print(f"Creating {natomtypes} fields for {natomtypes} found atom types")
    fields = []
    for itype in range(natomtypes):
      fields.append(Field(npw_Nd = npw, h=h))

    for frame_index in frames_to_average:
      print(f"Processing frame {frame_index} containing {frame.particles.N} atoms")
      # find specified frame
      frame = t[frame_index] 

      # loop through all particles in frame and call add_gaussian function
      for iatom in range(frame.particles.N):
          myP = P # all atom types currently use same sigma. Should be able to generalize...
          pos = frame.particles.position[iatom]
          atomtype_index = frame.particles.typeid[iatom]

          # TODO: calling this function within a double for-loop over particles and frames is killing performance
          #       it would be worth thinking about a more efficient implementation
          # initial idea: try pybinding this function or numba?
          add_hockney_eastwood_function(fields[atomtype_index], center=pos, P=myP, height=1)
          
          # TODO: in principle it should be possible to also initialize using Gaussians
          # however add_gaussians function currently searches over all grid points which makes it much to slow
          #add_gaussian(field, center=pos, sigma=mysigma, height=1)
    
    # fields contain total over ALL frames, need to compute average
    for field in fields:
      field.data /= nframes_to_average
      if normalize:
        field.data /= rho0 # I'm not sure if this normalization is quite right...seems to work reasonably though...

    # output (for debugging)
    #write_to_file("fields.dat",fields)
    #write_to_VTK("fields.vtk",fields)

    # return field
    return fields


def add_hockney_eastwood_function(field,center, P, height = 1):
    ''' add a Hockney Eastwood function of order P to field

    Args: 
        field: field to add Hockney Eastwood function to (edited in place)
               `units='real'` is default currently and `center` must fall within `h` cell tensor of field
        center: center of Hockney-Eastwood function (in real units). This is mapped into "index units" within function
        height: total integral of function TODO: probably better to rename this from "height" to something else
        
    '''
    npw = field.npw_Nd
    h = field.h # cell tensor
    dim = field.dim

    center = np.array(center) # recast center to array
    
    # box must be orthorhombic
    assert(field.is_orthorhombic()), f"Error: h must be orthorhombic in add_gaussian function. {h}"
    boxl = np.diag(h)
    
    # voxels must be cubic 
    grid_spacings = np.diag(h) / npw # since field is orthorhombic, gridspacings will be of length dim
    assert (np.allclose(grid_spacings, grid_spacings[0],atol=1e-2)), f"voxels must be cubic {grid_spacings = }"

    
    xparticle = center / boxl * npw # xparticle should be in range [0,npw)

    # check that 'center' and 'xparticle' are in appropriate range
    #for d in range(dim):
    #  assert(center[d] >= 0 and center[d] < boxl[d]), f"Invalid 'center' position specified {center = }"
    #  assert(xparticle[d] >= 0 and xparticle[d] < npw[d]), f"Invalid {xparticle = }"
     
    # apply PBC to xparticle
    for idim in range(dim):
        while xparticle[idim] <  0:         xparticle[idim] += npw[idim]
        while xparticle[idim] >= npw[idim]: xparticle[idim] -= npw[idim]

    grid_indicies, weights = hockney_eastwood_function(xparticle, npw, P)
   
    data = np.zeros(npw)
    for i,index in enumerate(grid_indicies):
        data[index] += weights[i]

    # now add to field
    field.data += data


    
def hockney_eastwood_function(xparticle, ngridpts, P):
    ''' Hockney-Eastwood Function as described in Deserno and Holm 1998

    Args:
        xparticle: position of particle in each dimension. Note that xparticle is given in "index units" so that spacking between grid points is = 1
        ngridpts: the number of grid point that make up the total grid   
        P: order of the assignment function

    Returns:
        grid_indicies: a list of tuples (each tuple is of lengh dim) that describe the grid indicies updated by the function
        weights: the weight of the particle corresponding to each of those grid indicies
    '''

    dim = len(xparticle)
    assert(dim >= 1 and dim <= 3) 
    
    grid_indicies_per_dim = [[]]*dim # blank 2d list, len(1) == dim, len(2) == unspecified
    weights_per_dim = np.zeros((dim,P))

    for idim in range(dim): 

      # set xbar
      if ((P % 2) == 0): # P is even
        # set xbar as midpoint of nearest two mesh points
        if np.floor(xparticle[idim]) != np.ceil(xparticle[idim]): 
          xbar = 0.5*(np.floor(xparticle[idim]) + np.ceil(xparticle[idim]))
          grid_indicies_per_dim[idim] = [int(np.floor(xbar)), int(np.ceil(xbar))]
          nelem = 2
        else: # check for edge case if xparticle is exactly an integer (so floor and ceil are equal)
          xbar = 0.5*(np.floor(xparticle[idim]) + np.ceil(xparticle[idim]) + 1) # use higher grid index, it wont matter since its weight will be ~0
          grid_indicies_per_dim[idim] = [int(np.floor(xbar)), int(np.ceil(xbar))+1] # should the 2nd entry be np.floor instead?
          nelem = 2
      else: # P is odd
        # set xbar as location of nearest mesh point
        xbar = np.round(xparticle[idim])
        grid_indicies_per_dim[idim] = [int(xbar)]
        nelem = 1

      # now find the indicies of the P closest mesh points
      while (len(grid_indicies_per_dim[idim]) < P):
        idx_down = grid_indicies_per_dim[idim][0] - 1
        idx_up = grid_indicies_per_dim[idim][-1] + 1
        grid_indicies_per_dim[idim].insert(0,idx_down) # add to beginning
        grid_indicies_per_dim[idim].append(idx_up)     # add to end
        nelem += 2
     
      grid_indicies_per_dim[idim] = np.array(grid_indicies_per_dim[idim],dtype=np.int64)
      
      # grid_indicies should be an ascending sequence of integers
      assert(len(grid_indicies_per_dim[idim]) == P)
      
      # apply PBCs (but don't change order)
      for i in range(P):
        while grid_indicies_per_dim[idim][i] < 0:               grid_indicies_per_dim[idim][i] += ngridpts[idim]
        while grid_indicies_per_dim[idim][i] >= ngridpts[idim]: grid_indicies_per_dim[idim][i] -= ngridpts[idim]
     
      
      # These equations are from Appendix E of Deserno and Holm
      dx = xparticle[idim]-xbar
      if P == 1:
        weights_per_dim[idim][0] = 1 
      elif P == 2:
        weights_per_dim[idim][0] = 0.5*(1-2*dx)
        weights_per_dim[idim][1] = 0.5*(1+2*dx)
      elif P == 3:
        weights_per_dim[idim][0] = 0.125*(1 - 4*dx + 4*dx*dx)
        weights_per_dim[idim][1] = 0.25*(3 - 4*dx*dx)
        weights_per_dim[idim][2] = 0.125*(1 + 4*dx + 4*dx*dx)
      elif P == 4:
        weights_per_dim[idim][0] = 1.0/48.*(1 - 6*dx + 12*dx*dx - 8*dx*dx*dx)
        weights_per_dim[idim][1] = 1.0/48.*(23 - 30*dx - 12*dx*dx + 24*dx*dx*dx)
        weights_per_dim[idim][2] = 1.0/48.*(23 + 30*dx - 12*dx*dx - 24*dx*dx*dx)
        weights_per_dim[idim][3] = 1.0/48.*(1 + 6*dx + 12*dx*dx + 8*dx*dx*dx)
      elif P == 5:
        weights_per_dim[idim][0] = 1.0/384.*(1 - 8*dx + 24*dx*dx - 32*dx*dx*dx + 16*dx*dx*dx*dx)
        weights_per_dim[idim][1] = 1.0/96. *(19 - 44*dx + 24*dx*dx + 16*dx*dx*dx - 16*dx*dx*dx*dx)
        weights_per_dim[idim][2] = 1.0/192.*(115 - 120*dx*dx + 48*dx*dx*dx*dx)
        weights_per_dim[idim][3] = 1.0/96. *(19 + 44*dx + 24*dx*dx - 16*dx*dx*dx - 16*dx*dx*dx*dx)
        weights_per_dim[idim][4] = 1.0/384.*(1 + 8*dx + 24*dx*dx + 32*dx*dx*dx + 16*dx*dx*dx*dx)
      else:
        raise RuntimeError ("Invalid P")
    
    # noW need to find all combinations of grid_indicies and weights     
    grid_indicies = []
    weights = []
    
    # by using np.ndindex, this should work for any dimension
    shape_Pdim = tuple([P]*dim)
    for i_Pdim in np.ndindex(shape_Pdim):
        index_Nd = []
        weight = 1.0
        for idim in range(dim):
            index_Nd.append(grid_indicies_per_dim[idim][i_Pdim[idim]])
            weight *= weights_per_dim[idim][i_Pdim[idim]]
        grid_indicies.append(tuple(index_Nd))
        weights.append(weight)
        
    return grid_indicies,weights





