import fieldkit as fk
from pytest import approx
import numpy as np

def test_2d():
    npw = (64,64)
    
    #h = np.array([[1,0,0],[0,1,0],[0,0,1]]) 
    h = np.array([[3,0],[0,1]]) 
    field = fk.Field(npw = npw, h = h)
    fk.add_gaussian(field, center=(0.0,0.0), sigma=0.10, height=1)
    fk.write_to_file("fields.dat",[field])

    assert(field.data[0,0] ==  approx(0.011656, abs=1e-4))
    assert(field.data[0,5] ==  approx(0.008591, abs=1e-4))
    assert(field.data[63,63] ==  approx(0.01031735809644507, abs=1e-4))
    assert(field.data[32,32] ==  approx(0.0, abs=1e-4))

