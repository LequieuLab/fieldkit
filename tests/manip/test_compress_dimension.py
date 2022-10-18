#!/usr/bin/env python3 

import fieldkit as fk
import numpy as np

#def test_3d_to_2d():
#    filename="slab_density.dat"
#
#    fields = fk.read_from_file(filename)
#
#    dim_new = 3
#    npw_new = 16
#    boxl_new = 5
#    fields_new = fk.expand_dimension(fields, dim_new, npw_new, boxl_new)
#   
#    h_target = np.diag([4.0, 6.928, boxl_new]) 
#    for field in fields_new:
#      assert(np.all(field.h == h_target))
#      assert(np.all(field.npw == np.array([32,32,16])))
#      assert(np.all(field.data.shape == np.array([32,32,16])))
#      for i in range(1,npw_new):
#        assert(np.all(field.data[:,:,i] == field.data[:,:,i-1])), "all of 3rd dimension shold be equal"

def test_3d_to_1d():
    filename="slab_density.dat"

    fields = fk.read_from_file(filename)

    dims_to_compress = [0,1] # compress x and y dimensions
    fields_new = fk.compress_dimension(fields, [0,1])
   
    for field, field_new in zip(fields,fields_new):
      assert(field_new.dim == 1)
      assert(field.h[2,2] == field_new.h[0,0])
      assert(field.npw[2] == field_new.npw[0])



if __name__ == '__main__':
  #test_3d_to_2d()
  test_3d_to_1d()
 

