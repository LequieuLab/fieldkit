import fieldkit as fk
from pytest import approx
import numpy as np

def test():

    filename = '../density_slab_cl.dat'
    fields = fk.read_from_file(filename)
    state = fk.binodal_from_slab(fields, field_idx_for_interfaces=0, plot=False)
    print(state) 
    # TODO: ADD TESTS
       
if __name__ == "__main__":
    test()
