import fieldkit as fk
from pytest import approx
import numpy as np

def test_2d():
    npw = (64,64)
    h = np.array([[3,0],[0,1]]) 
    
    field = fk.Field(npw_Nd = npw, h=h)

    fk.add_gaussian(field, center=(0.0,0.0), sigma = 0.2, height=100)
    fk.add_gaussian(field, center=(0.5,0.5), sigma = 0.2, height=100)
    fk.add_gaussian(field, center=(2.0,0.5), sigma = 0.2, height=100)

    fk.write_to_file("fields.dat",[field])
   
    threshold = 0.25
    print(fk.calc_domain_stats_mesh(field,threshold,plotMesh=False,outputMesh=False,add_periodic_domains=False, applyPBC=True))



def test_3d():
    npw = (32,32,32)
    h = np.array([[3,0,0],[0,1,0],[0,0,2]]) 
    
    field = fk.Field(npw_Nd = npw, h=h)

    fk.add_gaussian(field, center=(0.0,0.0,0.0), sigma = 0.2, height=100)
    fk.add_gaussian(field, center=(0.5,0.5,0.0), sigma = 0.2, height=100)
    fk.add_gaussian(field, center=(2.0,0.5,0.0), sigma = 0.2, height=100)

    fk.write_to_VTK("fields.vtk",[field])
    
    threshold = 0.1
    print(fk.calc_domain_stats_mesh(field,threshold,plotMesh=False,outputMesh=False,add_periodic_domains=False, applyPBC=True))
    
if __name__ == "__main__":
    test_3d()
