import fieldkit as fk
from pytest import approx
import numpy as np

def test_hockney_eastwood_3d():
    #npw = (32,32,32)
    npw = (16,16,16)
    
    gsdfile = "./initialize/micelle.gsd"
    #trjfile = "./initialize/lam.lammpstrj"
    #psffile = "./initialize/lam.psf"
    frame_index = -1 # use last frame
    P = 2
    fields = fk.particle_to_field_gsd(gsdfile, frame_index, npw, P, normalize=True)

    fk.write_to_VTK("fields.vtk",fields)
    
    # FIXME add tests
    assert(True)


if __name__ == "__main__":
    test_hockney_eastwood_3d()

