#!/usr/bin/env python3 

import fieldkit as fk
import numpy as np

def test_2d_to_3d():
    filename="density.dat"

    fields = fk.read_from_file(filename)

    dim_new = 3
    npw_new = 16
    boxl_new = 5
    fields_new = fk.expand_dimension(fields, dim_new, npw_new, boxl_new)
   
    h_target = np.diag([4.0, 6.928, boxl_new]) 
    for field in fields_new:
      assert(np.all(field.h == h_target))
      assert(np.all(field.npw == np.array([32,32,16])))
      assert(np.all(field.data.shape == np.array([32,32,16])))
      for i in range(1,npw_new):
        assert(np.all(field.data[:,:,i] == field.data[:,:,i-1])), "all of 3rd dimension shold be equal"

def test_1d_to_3d():
    filename="density_1D.dat"

    fields = fk.read_from_file(filename)

    dim_new = 3
    npw_new = [16,32]
    boxl_new = [5,3]
    fields_new = fk.expand_dimension(fields, dim_new, npw_new, boxl_new)
   
    h_target = np.diag([4.26311, 5, 3]) 
    for field in fields_new:
      #breakpoint()
      assert(np.all(np.isclose(field.h, h_target,atol=1e-2)))
      assert(np.all(field.npw == np.array([64,16,32])))
      assert(np.all(field.data.shape == np.array([64,16,32])))
      for i in range(1,npw_new[1]):
        for j in range(1,npw_new[0]):
            assert(np.all(field.data[:,j,i] == field.data[:,j-1,i-1])), "all of 2nd and 3rd dimension shold be equal"



if __name__ == '__main__':
  test_2d_to_3d()
  test_1d_to_3d()
 

