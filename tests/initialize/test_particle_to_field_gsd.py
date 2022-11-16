import fieldkit as fk
from pytest import approx
import numpy as np
import json

def test_micelle_3d_by_species():
    npw = (16,16,16)
    
    gsdfile = "./initialize/micelle/micelle.gsd"
    frame_index = -1 # use last frame
    P = 2
    fields = fk.particle_to_field_gsd(gsdfile, frame_index, npw, P, normalize=True)

    fk.write_to_file("fields_species.dat",fields)
    fk.write_to_VTK("fields_species.vtk",fields)
    
    # FIXME add tests
    assert(True)

def test_micelle_3d_by_molecule():
    npw = (16,16,16)
    
    gsdfile = "./initialize/micelle/micelle.gsd"
    frame_index = -1 # use last frame
    P = 2
    types = np.zeros(100000,dtype=int) # hardcoded 1e5 atoms in micelle.gsd
    types[14999:] = 1 # first 1207 atoms are polymers (moltype=0), rest are solvent (moltype=1)
    fields = fk.particle_to_field_gsd(gsdfile, frame_index, npw, P, types = types, normalize=True)

    fk.write_to_file("fields_molecule.dat",fields)
    fk.write_to_VTK("fields_molecule.vtk",fields)
    
    # FIXME add tests
    assert(True)



def test_lam_3d():
    #npw = (32,32,32)
    npw = (16,16,16)
    
    gsdfile = "./initialize/lamella/final.gsd"
    frame_index = -1 # use last frame
    P = 4
    fields = fk.particle_to_field_gsd(gsdfile, frame_index, npw, P, normalize=False)

    fk.write_to_file("fields.dat",fields)
    fk.write_to_VTK("fields.vtk",fields)
    
    #read json file
    jsonfile = "./initialize/lamella/polymer.json"
    with open(jsonfile, "r") as rf:
      json_data = json.load(rf)
 
    # add tests
    # check that total sum = natoms
    assert(np.isclose(json_data['config']['natoms'], np.sum(fields[0].data + fields[1].data))) 

    # check cell vs h
    h = fields[0].h
    assert(np.isclose(json_data['config']['cell'][0],h[0][0])) 
    assert(np.isclose(json_data['config']['cell'][1],h[1][1])) 
    assert(np.isclose(json_data['config']['cell'][2],h[2][2])) 
    assert(np.isclose(json_data['config']['cell'][3],h[0][1])) 
    assert(np.isclose(json_data['config']['cell'][4],h[0][2])) 
    assert(np.isclose(json_data['config']['cell'][5],h[1][2])) 

    # add other tests

if __name__ == "__main__":
    test_micelle_3d_by_molecule()
    #test_micelle_3d_by_species()
    #test_lam_3d()

