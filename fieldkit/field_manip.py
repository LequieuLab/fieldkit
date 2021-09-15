import numpy as np
import copy

from .field import *

def change_resolution(fields_old,resolution_new):
  """For a list of Field objects, change the resolution of each Field object.

    Args:
        field_old: a list of Field objects
        resolution_new: a tuple that defines the new resolution of each Field object.
   
    Returns:
        field_new: a list of Field objects with each object set to a resolution of resolution_new.
        
  """
  fields_new = []
  for field_old in fields_old: 
    
    npw_old = field_old.npw_Nd
    npw_new = np.array(resolution_new)
    dim = len(npw_old)
    assert(dim == len(npw_new)), "dimension of fields must match resoltion"
    
    data_old = field_old.data
    data_old_k = np.fft.fftn(data_old)

    data_new_k = np.zeros(npw_new,dtype=np.complex) 

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
        kindex_new = np.zeros(dim,dtype=np.int)
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
    field_new = Field(h=field_old.h, npw_Nd=npw_new, data=data_new)
    fields_new.append(field_new)

  return fields_new

def replicate_fields(fields, nreplicates):    
    """ For a list of Fields, replicate each Field object by nreplicates.
    Adapted from FTS-tools/replicate_fields.py and FTS-tools/lib/fieldtools.py.
    
    Args:
        fields: a list of Field objects
        nreplicates: number of replicates

    Returns: 
        fields_list: a list of Field objects, in which each Field object is replicated by nreplicates amount of times. 

    """

    field_list = []
    
    #replicate each Field object by nreplicates
    for f in range(len(fields)):
        dim = fields[0].dim
        npw = fields[0].npw_Nd
        h = fields[0].h
    
        if dim == 1:
            rep = nreplicates
            h = h * nreplicates
        else:
            assert(len(nreplicates) == dim)    
            rep = list(nreplicates)
            
        fields_new = np.tile(fields[f].data.real, rep)
        npw_new = np.array(npw) * nreplicates
        field_obj = Field(npw_Nd=npw_new, data=fields_new, h=h)
        

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

    fields_new = []
    for field in fields:    
        npw = field.npw_Nd
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




