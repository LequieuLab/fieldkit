'''
Functions to manipulate fields
'''

import numpy as np
import copy

from .field import *

def change_resolution(fields,npw_new):
  """For a list of Field objects, change the resolution of each Field object.

    Args:
        field: a list of Field objects
        npw_new: a tuple that defines the new resolution of each Field object.
   
    Returns:
        field_new: a list of Field objects with each object set to a resolution of npw_new.
        
  """

  # if single entry, convert it to a list
  # TODO: I should probably adda utility function that does this
  try:
    nfields = len(fields)
  except TypeError:
    fields = [fields]
    nfields = len(fields)
 
  fields_new = []
  for field_old in fields: 
    
    npw_old = field_old.npw
    npw_new = np.array(npw_new)
    dim = len(npw_old)
    assert(dim == len(npw_new)), "dimension of fields must match resoltion"
    
    data_old = field_old.data
    data_old_k = np.fft.fftn(data_old)

    data_new_k = np.zeros(npw_new,dtype=complex) 

    # create map
    for kindex_old in np.ndindex(npw_old):
   
      # convert kindex (std FFT ordering) to kprime (unaliased)
      kprime_old = np.array(kindex_old)
      for i,k in enumerate(kindex_old):
        if k <= 0.5*npw_old[i]-1:
          kprime_old[i] = k
        else:
          kprime_old[i] = k - npw_old[i]
     
      if np.any(kprime_old >= 0.5*npw_new) or np.any(kprime_old < -0.5*npw_new):
        # dont copy into data_new, these are frequencies that are cut
        pass
      else: 
        # convert kprime (unaliased) to kindex (std FFT ordering)
        kprime_new = kprime_old
        kindex_new = np.zeros(dim,dtype=np.int64)
        for i,k in enumerate(kprime_new):
          if k >= 0:
            kindex_new[i] = k
          else:
            kindex_new[i] = k + npw_new[i]

        kindex_new = tuple(kindex_new)
        data_new_k[kindex_new] = data_old_k[kindex_old]
      
    # convert back to real space 
    data_new = np.fft.ifftn(data_new_k) 
    # need to renormalize data since npw changed from fft and ifft
    data_new *= np.prod(npw_new/npw_old)

    # create a new field
    field_new = Field(h=field_old.h, npw=npw_new, data=data_new)
    fields_new.append(field_new)

  return fields_new

def cubic_voxels(fields_old,tol=0.01):
    """For a list of Field objects, change the resolution of each Field object so that the voxels are approximately cubic

      Args:
          field_old: a list of Field objects
          tol: tolerable error for how cubic you want the voxels to be
     
      Returns:
          field_new: a list of Field objects with each object set to a new resolution that is approximately cubic
          
    """

    # perform checks on all fields
    for i in range(len(fields_old)):
        assert(fields_old[i].is_orthorhombic), "For voxels to be cubic, cell must be orthorhombic"
        assert(fields_old[i].dim == 3),"dimension must be 3"
        if i != 1:
          assert(np.all(fields_old[i].h == fields_old[i-1].h)),"h must be equal for all fields"
          assert(np.all(fields_old[i].npw == fields_old[i-1].npw)),"h must be equal for all fields"

    # since all fields have same h and npw, just use first field
    field_old = fields_old[0] 

    L = np.diag(field_old.h) # L will be constant     
    npw = np.array(field_old.npw) # npw will vary within while loop 
    #print(f"{L = }")

    # preliminaries 
    L_per_voxel = L / npw
    diff_per_dim = _calc_cubic_voxel_difference(L_per_voxel)
    error = np.max(np.abs(diff_per_dim))
    npw_attempts = [list(npw)] 

    while  error > tol:
        error_per_dim = np.abs(diff_per_dim)
        max_error_index = error_per_dim.argmax()
        if diff_per_dim[max_error_index] < 0:
            # max error voxel length is too small. need to DECREASE npw in this dim
            npw[max_error_index] -= 1
        else:
            # max error voxel length is too BIG. need to INCREASE npw in this dim
            npw[max_error_index] += 1

        # if this npw has already been tried, then multiply by two and continue
        if list(npw) in npw_attempts:
            npw *= 2 # double all npw and continue

        # now update error
        L_per_voxel = L / npw
        diff_per_dim = _calc_cubic_voxel_difference(L_per_voxel)
        error = np.max(np.abs(diff_per_dim))
        npw_attempts.append(list(npw)) # not casting as list, this makes comparison above easier
        #print(f"{npw = } {L_per_voxel=} {diff_per_dim=} {error=}")
      
    # using the updated npw create new fields
    fields_new = change_resolution(fields_old,npw)
    print(f"To make cubic voxels, new grid resolution is {npw}, voxel lengths are {L_per_voxel}, error = {error : .2e}")
    return fields_new
        
def _calc_cubic_voxel_difference(L_per_voxel):
    ''' Helper function to calculate how "uncubic" a voxel is
        Args:
            L_per_voxel: array storing the length of each side of voxel. should have 'dim' elements
        Returns: the error of each dimension (here the absolute value of the deviation of from the median)
    '''
    Lmedian = np.average(L_per_voxel)
    return L_per_voxel - Lmedian

def replicate_fields(fields, nreplicates):    
    """ For a list of Fields, replicate each Field object by nreplicates.
    Adapted from FTS-tools/replicate_fields.py and FTS-tools/lib/fieldtools.py.
    
    Args:
        fields: a list of Field objects
        nreplicates: number of replicates

    Returns: 
        fields_list: a list of Field objects, in which each Field object is replicated by nreplicates amount of times. 

    """

    # if single entry, convert it to a list
    # TODO: I should probably adda utility function that does this
    try:
      nfields = len(fields)
    except TypeError:
      fields = [fields]
      nfields = len(fields)
   
    field_list = []
    
    #replicate each Field object by nreplicates
    for f in range(len(fields)):
        dim = fields[0].dim
        npw = fields[0].npw
        h = fields[0].h
    
        if dim == 1:
            rep = nreplicates
            h = h * nreplicates
        else:
            assert(len(nreplicates) == dim)    
            rep = list(nreplicates)
            
        fields_new = np.tile(fields[f].data.real, rep)
        npw_new = np.array(npw) * nreplicates
        field_obj = Field(npw=npw_new, data=fields_new, h=h)
        

        coords_new = field_obj.CoordsFromH()            
        field_obj.set_coords(coords_new)
        field_list.append(field_obj)

    return field_list
        
def roll(fields, shift):
    """ roll the fields across the PBC using shift

    Args:
        fields: a list of Field objects
        shift: list of length dim. How much to roll the fields. In range [0-1].

    Returns: 
        fields_new: a list of Field objects, in which each Field object has been translated. 

    """
    # if single entry, convert it to a list
    # TODO: I should probably adda utility function that does this
    try:
      nfields = len(fields)
    except TypeError:
      fields = [fields]
      nfields = len(fields)
 

    fields_new = []
    for field in fields:    
        npw = field.npw
        dim = field.dim
        assert(len(shift) == dim)

        field_new = copy.deepcopy(field) # ensure deep copy
        
        myshift=np.zeros(dim,dtype=int) # must be an integer
        myaxis = tuple(np.arange(dim))
        for idim in myaxis:
            myshift[idim] = int(shift[idim]*npw[idim])
        
        field_new.data = np.roll(field.data, myshift, axis=myaxis)
        fields_new.append(field_new)
   
    return fields_new 


def expand_dimension(fields, dim_new, npw_new, boxl_new):
    """Convert a field to one of a higher dimension (e.g. 1d to 2d or 3d)

        Args:
          fields: list of field arguments
          dim_new (list): new dimension of fields (must be higher than input dimension)
          npw_new (list): npw for each new dimension created
          boxl_new (list): box length for each new dimension created
        
        Returns:
          fields_new: a new list of fields with expanded dimension      

    """

    assert(type(dim_new) == int), "dim_new must be int"

    # convert to lists if needed
    if type(npw_new) == int: npw_new = [npw_new]
    if type(boxl_new) == int: boxl_new = [boxl_new]
    if type(fields) == Field: field = [fields]

    fields_new = []
    for field in fields:
      dim = field.dim
      assert(dim_new <= 3), "new dimension must be <= 3"
      assert(dim <= dim_new), "old dimension must be < new_dim"
      assert(len(boxl_new) == dim_new - dim), "boxl_new needed for each new dimension"
      assert(len(npw_new) == dim_new - dim), "boxl_new needed for each new dimension"
      assert(field.is_orthorhombic()), "field must be orthorhombic"
   
      nexpand = dim_new - dim
      npw = list(field.npw) 
      boxl = list(np.diag(field.h))
      for i in range(nexpand) :
        npw.append(npw_new[i])
        boxl.append(boxl_new[i])
      h = np.diag(boxl)

      if nexpand == 1:
        if dim_new == 2:
          data = np.tile(field.data,(npw_new[0],1)).T
        if dim_new == 3:
          data = np.tile(field.data.T,(npw_new[0], 1 ,1)).T
      if nexpand == 2:
        assert(dim_new == 3 and dim == 1), "only way to expand by 2 is from 1d to 3d"
        data = np.tile(field.data.T,(npw_new[1], npw_new[0] ,1)).T

      #print(f"{npw = }, {h = } {data.shape = }")
      #breakpoint()
      fields_new.append(Field(npw=npw, h=h, data=data))

    return fields_new


def compress_dimension(fields, dims_to_compress):
    """Convert a field to one of a lower dimension (e.g. 1d to 2d or 3d)

        Args:
          fields: list of field arguments
          dims_to_compress (list): dimensions to compress
        
        Returns:
          fields_new: a new list of fields with compressed dimension      


    """

    # if single entry, convert to list
    if type(fields) == Field:
      fields = [fields]
    elif type(fields) == list:
      pass
    else:
      raise TypeError(f"Invalid type of fields {type(fields) = }")


    # now compress each field
    fields_new = []
    for field in fields:
      assert(field.is_orthorhombic), "compress_dimension requires fields to be orthorhombic"

      # extract npw and h from preserved dimensions
      npw = []
      h_diag = []
      for idim in range(field.dim):
        if not (idim in dims_to_compress):
          h_diag.append(field.h[idim,idim])
          npw.append(field.npw[idim])
      h = np.diag(h_diag) 

      # average data along dims to compress 
      data = np.average(field.data, axis=tuple(dims_to_compress))

      # create new field
      fields_new.append(Field(npw=npw, h=h, data=data))

    return fields_new


def approx_3dsphere_from_1d(fields_1d, npw_3d):
    """
       Transform a 1D field into a 3D field assuming a spherical geometry. 
       Further assumes that each field in fields_1d has the same cell tensor 
       and npw, that the field values are symmetrical across the cetner, and
       that any phase separated dense phase is at the center of the field's 
       coordinates.

       Args:
         fields_1d: list of 1D fields
         npw_3d: npw for the 3D fields

       Returns:
         fields_3d: list of 3D fields generated from the 1D fields
    """

    # Create a list for the 3D fields.

    fields_3d = []

    # Get the coordinates at the center of the 3D field and 1D field, and the
    # distance from the center for each 1D grid point.

    center_1d = [x // 2 for x in fields_1d[0].npw]
    r0_1d = fields_1d[0].coords[center_1d[0]]
    dists_1d = np.abs(fields_1d[0].coords - r0_1d)

    # Get the 1D cell tensor and make the 3D cell tensor.

    h_1d = fields_1d[0].h[0][0]
    h_3d = np.asarray([[h_1d, 0, 0],
                       [0, h_1d, 0],
                       [0, 0, h_1d]])

    # Loop through each species in the 1D field.

    nspecies = len(fields_1d)
    for l in range(nspecies):

        # Create a 3D field of the appropriate size and grid density.

        field_3d = Field(npw=npw_3d, h=h_3d)
        field_3d.data = field_3d.data.astype(type(fields_1d[l].data[0]))

        # Get the coordinates at the center of the 3D field and 1D field.

        center_3d = [x // 2 for x in npw_3d]
        r0_3d = field_3d.coords[center_3d[0]][center_3d[1]][center_3d[2]]

        # Loop through each 3D grid point.

        for i in range(npw_3d[0]):
            for j in range(npw_3d[1]):
                for k in range(npw_3d[2]):

                    # Get the distance between the center and the point.

                    r_3d = field_3d.coords[i][j][k]
                    dist_3d = np.linalg.norm(r_3d - r0_3d)

                    # Get the value from the 1D field at the same distance.

                    diff = np.abs(dist_3d - dists_1d)[:, 0]
                    inds = np.where(diff == diff.min())[0]
                    val = 0
                    for ind in inds:
                        val += fields_1d[l].data[ind]
                    val /= inds.size

                    # Set the data value in the 3D field.

                    field_3d.data[i][j][k] = val

        # Add the 3D field to the list of 3D fields.

        fields_3d.append(field_3d)

    # Return the list of 3D fields.

    return fields_3d

