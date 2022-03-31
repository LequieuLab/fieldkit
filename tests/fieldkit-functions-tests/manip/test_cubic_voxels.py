#!/usr/bin/env python3 

import fieldkit as fk 
from pytest import approx
import numpy as np

def test_cubic_voxels():

    filename = "density-hex-3d.dat"
    fields = fk.read_from_file(filename)
    new_fields = fk.cubic_voxels(fields, tol=1e-2)
    for new_field in new_fields:
      assert(np.all(new_field.npw == np.array([27,42,32])))

    new_fields = fk.cubic_voxels(fields, tol=1e-3)
    for new_field in new_fields:
      assert(np.all(new_field.npw == np.array([51,88,64])))

if __name__ == '__main__':
    test_cubic_voxels()
