import fieldkit as fk
from pytest import approx
import numpy as np

def test_hockney_eastwood_3d():
    #npw = (32,32,32)
    npw = (16,16,16)
    
    trjfile = "frame.lammpstrj"
    psffile = "frame.psf"
    frame_index = -1 # use last frame
    P = 2
    fields = fk.particle_to_field_hockney_eastwood(trjfile,psffile, frame_index, npw, P)

    fk.write_to_VTK("fields.vtk",fields)
    
    # FIXME add tests


if __name__ == "__main__":
    test_hockney_eastwood_3d()

