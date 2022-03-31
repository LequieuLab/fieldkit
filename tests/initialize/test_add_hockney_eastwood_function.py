import fieldkit as fk
from pytest import approx
import numpy as np

def test_2d():
    npw = (32,16)
    
    #h = np.array([[1,0,0],[0,1,0],[0,0,1]]) 
    h = np.array([[2,0],[0,1]]) 
    field = fk.Field(npw = npw, h = h)
    
    # check P=5
    fk.add_hockney_eastwood_function(field, center=(0.0,0.0), P=5, height=1)
    assert(np.sum(field.data) ==  approx(1.0, abs=1e-4)) # confirm that integral = 1
    assert(field.data[0,0] ==  approx(0.35875108506944436, abs=1e-4))
    assert(field.data[0,1] == field.data[1,0] == field.data[-1,0] == approx(0.11854383680555554, abs=1e-4))
    assert(field.data[0,2] ==  approx(0.0015597873263888886, abs=1e-4))
    assert(field.data[0,3] ==  approx(0.0, abs=1e-4))

    # check P=3
    fk.add_hockney_eastwood_function(field, center=(0.5,0.5), P=3, height=1)
    assert(np.sum(field.data) ==  approx(2.0, abs=1e-4)) # confirm that integral = 2
    assert(field.data[8,8] ==  approx(0.5625, abs=1e-4)) 
    assert(field.data[7,8] == field.data[9,8] == field.data[8,7] == field.data[8,9] ==  approx(0.09375, abs=1e-4)) 
    assert(field.data[6,8] == field.data[10,8] == field.data[8,6] == field.data[8,10] ==  approx(0.00, abs=1e-4)) 
    assert(field.data[7,7] == field.data[9,7] == field.data[7,9] == field.data[9,9] ==  approx(0.015625, abs=1e-4)) 
    
    # check P=2
    fk.add_hockney_eastwood_function(field, center=(0.45,0), P=2, height=1)
    assert(np.sum(field.data) ==  approx(3.0, abs=1e-4)) # confirm that integral = 3
    assert(field.data[8,0] ==  approx(0.2, abs=1e-4)) 
    assert(field.data[7,0] ==  approx(0.8, abs=1e-4)) 
    assert(field.data[7,1] == field.data[7,-1] ==  approx(0.0, abs=1e-4)) 

    fk.write_to_file("fields.dat",[field])
    
if __name__ == "__main__":
    test_2d()
