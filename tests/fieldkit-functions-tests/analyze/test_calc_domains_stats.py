import fieldkit as fk
from pytest import approx
import numpy as np

def test_2d_cubic():
    npw = (32,32)
    h = np.array([[1,0],[0,1]]) 
    field = fk.Field(npw_Nd = npw, h=h)
    
    sig = 0.15
    r = np.sqrt(2*np.log(2)) * sig # = half width at half max
    c = [[0,0], [0.5,0.5]]
    fk.add_gaussian(field, center=c[0], sigma=sig, height=1)
    fk.add_gaussian(field, center=c[1], sigma=sig, height=1)
    field.data /= np.max(field.data) # normalize so max == 0 

    fk.write_to_file("fields.dat",[field])
   
    threshold = 0.5
    stats = fk.calc_domain_stats(field,threshold,plotMesh=True,outputMesh=False,add_periodic_domains=False, applyPBC=True)
    
    assert(stats['ndomains'] == 2) 
    for i in range(stats['ndomains']):
        assert(stats['center'][i] == approx(c[i],1e-2))
        assert(stats['surface_area'][i] == approx(2*np.pi*r,1e-2))
        assert(stats['volume'][i] == approx(np.pi*r**2,1e-2))
        assert(stats['IQ'][i] == approx(1.0,1e-2))

def test_2d_rectangular():
    npw = (32,64)
    h = np.array([[1,0],[0,2]]) 
    
    field = fk.Field(npw_Nd = npw, h=h)
    
    sig = 0.20
    r = np.sqrt(2*np.log(2)) * sig # = half width at half max
    c = [[0,0],[0.5,1]]
    fk.add_gaussian(field, center=c[0], sigma=sig, height=1)
    fk.add_gaussian(field, center=c[1], sigma=sig, height=1)
    field.data /= np.max(field.data) # normalize so max == 0 

    fk.write_to_file("fields.dat",[field])
   
    threshold = 0.5
    stats = fk.calc_domain_stats(field,threshold,plotMesh=True,outputMesh=False,add_periodic_domains=False, applyPBC=True)
    print(stats) 
    
    assert(stats['ndomains'] == 2) 
    for i in range(stats['ndomains']):
        assert(stats['center'][i] == approx(c[i],1e-2))
        assert(stats['surface_area'][i] == approx(2*np.pi*r,1e-2))
        assert(stats['volume'][i] == approx(np.pi*r**2,1e-2))
        assert(stats['IQ'][i] == approx(1.0,1e-2))


def test_3d_cubic():
    npw = (32,32,32)
    h = np.array([[2,0,0],[0,2,0],[0,0,2]]) 
    
    field = fk.Field(npw_Nd = npw, h=h)

    sig = 0.20
    r = np.sqrt(2*np.log(2)) * sig # = half width at half max
    c = [[0,0,0], [1,1,1]]
    fk.add_gaussian(field, center=c[0], sigma=sig, height=1)
    fk.add_gaussian(field, center=c[1], sigma=sig, height=1)
    field.data /= np.max(field.data) # normalize so max == 0 

    fk.write_to_VTK("fields.vtk",[field])
    
    threshold = 0.5
    stats = fk.calc_domain_stats(field,threshold,plotMesh=True,outputMesh=False,add_periodic_domains=False, applyPBC=True)
    print(stats)
    
    assert(stats['ndomains'] == 2) 
    for i in range(stats['ndomains']):
        assert(stats['center'][i] == approx(c[i],1e-2))
        assert(stats['surface_area'][i] == approx(4*np.pi*r**2,abs=1e-1))
        assert(stats['volume'][i] == approx(4/3*np.pi*r**3,rel=1e-1))
        assert(stats['IQ'][i] == approx(1.0,3e-2))

def test_3d_orthorhombic():
    npw = (32,32,64)
    h = np.array([[2,0,0],[0,2,0],[0,0,4]]) 
    
    field = fk.Field(npw_Nd = npw, h=h)

    sig = 0.20
    r = np.sqrt(2*np.log(2)) * sig # = half width at half max
    c = [[0,0,0], [1,1,1]]
    fk.add_gaussian(field, center=c[0], sigma=sig, height=1)
    fk.add_gaussian(field, center=c[1], sigma=sig, height=1)
    field.data /= np.max(field.data) # normalize so max == 0 

    fk.write_to_VTK("fields.vtk",[field])
    
    threshold = 0.5
    stats = fk.calc_domain_stats(field,threshold,plotMesh=True,outputMesh=False,add_periodic_domains=False, applyPBC=True)
    print(stats)
    
    assert(stats['ndomains'] == 2) 
    for i in range(stats['ndomains']):
        assert(stats['center'][i] == approx(c[i],1e-2))
        assert(stats['surface_area'][i] == approx(4*np.pi*r**2,abs=1e-1))
        assert(stats['volume'][i] == approx(4/3*np.pi*r**3,rel=1e-1))
        assert(stats['IQ'][i] == approx(1.0,3e-2))


   
if __name__ == "__main__":
    test_3d_orthorhombic()
