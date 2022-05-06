import fieldkit as fk
from pytest import approx
import numpy as np

def test_1d_P2():
  
  P = 2
  ngridpts = np.array([32])
  
  # xparticle = 0
  xparticle = np.array([0])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P)
  assert(np.sum(weights) == 1)
  assert(np.all(grid_indicies == np.array([[0],[1]])))
  assert(np.all(np.isclose(weights,[1.0,0.0])))

  # xparticle = 0.5
  xparticle = np.array([0.5])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P)
  assert(np.sum(weights) == 1)
  assert(np.all(grid_indicies == np.array([[0],[1]])))
  assert(np.all(np.isclose(weights,[0.5,0.5])))

  # xparticle = 15
  xparticle = np.array([15])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P)
  assert(np.sum(weights) == 1)
  assert(np.all(grid_indicies == np.array([[15],[16]])))
  assert(np.all(np.isclose(weights,[1.0,0.0])))

  # xparticle = 31.8
  xparticle = np.array([31.8])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P)
  assert(np.sum(weights) == 1)
  assert(np.all(grid_indicies == np.array([[31],[0]])))
  assert(np.all(np.isclose(weights,[0.2,0.8])))

def test_1d_P3():
  
  P = 3
  ngridpts = np.array([32])
  
  # xparticle = 0
  xparticle = np.array([0])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P)
  assert(np.sum(weights) == 1)
  assert(np.all(grid_indicies == np.array([[31],[0],[1]])))
  assert(np.all(np.isclose(weights,[0.125,0.75,0.125])))

  # xparticle = 0.5
  xparticle = np.array([0.5])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P)
  assert(np.sum(weights) == 1)
  assert(np.all(grid_indicies == np.array([[31],[0],[1]])))
  assert(np.all(np.isclose(weights,[0.0,0.5,0.5])))

  # xparticle = 15
  xparticle = np.array([15])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P)
  assert(np.sum(weights) == 1)
  assert(np.all(grid_indicies == np.array([[14],[15],[16]])))
  assert(np.all(np.isclose(weights,[0.125,0.75,0.125])))

  # xparticle = 31.8
  xparticle = np.array([31.8])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P)
  assert(np.sum(weights) == 1)
  assert(np.all(grid_indicies == np.array([[31],[0],[1]])))
  assert(np.all(np.isclose(weights,[0.245,0.71,0.045])))

def test_1d_P4():
  
  P = 4
  ngridpts = np.array([32])
  
  # xparticle = 0
  xparticle = np.array([0])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P)
  assert(np.isclose(np.sum(weights),1))
  assert(np.all(grid_indicies == np.array([[31],[0],[1],[2]])))
  assert(np.all(np.isclose(weights,[0.166666,0.666666,0.166666,0.0],atol=1e-4)))

  # xparticle = 0.5
  xparticle = np.array([0.5])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P)
  assert(np.isclose(np.sum(weights),1))
  assert(np.all(grid_indicies == np.array([[31],[0],[1],[2]])))
  assert(np.all(np.isclose(weights,[0.02083333, 0.47916667, 0.47916667, 0.02083333])))

  # xparticle = 15
  xparticle = np.array([15])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P)
  assert(np.isclose(np.sum(weights),1))
  assert(np.all(grid_indicies == np.array([[14],[15],[16],[17]])))
  assert(np.all(np.isclose(weights,[0.166666,0.666666,0.166666,0.0],atol=1e-4)))

  # xparticle = 31.8
  xparticle = np.array([31.8])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P)
  assert(np.isclose(np.sum(weights),1))
  assert(np.all(grid_indicies == np.array([[30],[31],[0],[1]])))
  assert(np.all(np.isclose(weights,[0.00133333, 0.28266667, 0.63066667, 0.08533333])))

def test_2d_P2():
  
  P = 2
  dim = 2
  ngridpts = np.array([32,16])
  
  # xparticle = 0,0
  xparticle = np.array([0,0])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P**dim)
  assert(np.isclose(np.sum(weights),1))
  assert(np.all(grid_indicies == np.array([[0,0],[0,1],[1,0],[1,1]])))
  assert(np.all(np.isclose(weights,[1,0,0,0],atol=1e-4)))

  # xparticle = 0.5,0.5
  xparticle = np.array([0.5,0.5])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P**dim)
  assert(np.isclose(np.sum(weights),1))
  assert(np.all(grid_indicies == np.array([[0,0],[0,1],[1,0],[1,1]])))
  assert(np.all(np.isclose(weights,[0.25,0.25,0.25,0.25],atol=1e-4)))

  # xparticle = 0.0,0.5
  xparticle = np.array([0.0,0.5])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P**dim)
  assert(np.isclose(np.sum(weights),1))
  assert(np.all(grid_indicies == np.array([[0,0],[0,1],[1,0],[1,1]])))
  assert(np.all(np.isclose(weights,[0.5,0.5,0.0,0.0],atol=1e-4)))

  # xparticle = 0.5,0.0
  xparticle = np.array([0.5,0.0])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P**dim)
  assert(np.isclose(np.sum(weights),1))
  assert(np.all(grid_indicies == np.array([[0,0],[0,1],[1,0],[1,1]])))
  assert(np.all(np.isclose(weights,[0.5,0.0,0.5,0.0],atol=1e-4)))

  # xparticle = 15,15
  xparticle = np.array([15,15])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P**dim)
  assert(np.isclose(np.sum(weights),1))
  assert(np.all(grid_indicies == np.array([[15,15],[15,0],[16,15],[16,0]])))
  assert(np.all(np.isclose(weights,[1,0,0,0],atol=1e-4)))

  # xparticle = 15.5,15.5
  xparticle = np.array([15.5,15.5])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P**dim)
  assert(np.isclose(np.sum(weights),1))
  assert(np.all(grid_indicies == np.array([[15,15],[15,0],[16,15],[16,0]])))
  assert(np.all(np.isclose(weights,[0.25,0.25,0.25,0.25],atol=1e-4)))

def test_2d_P3():
  
  P = 3
  dim = 2
  ngridpts = np.array([32,16])
  
  # xparticle = 0,0
  xparticle = np.array([0,0])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P**dim)
  assert(np.isclose(np.sum(weights),1))
  assert(np.all(grid_indicies == np.array([[31,15],[31,0],[31,1],[0,15],[0,0],[0,1],[1,15],[1,0],[1,1]])))
  assert(np.all(np.isclose(weights,[0.015625, 0.09375 , 0.015625, 0.09375 , 0.5625  , 0.09375 , 0.015625, 0.09375 , 0.015625], atol=1e-4)))

  # xparticle = 0.5,0.5
  xparticle = np.array([0.5,0.5])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P**dim)
  assert(np.isclose(np.sum(weights),1))
  assert(np.all(grid_indicies == np.array([[31,15],[31,0],[31,1],[0,15],[0,0],[0,1],[1,15],[1,0],[1,1]])))
  assert(np.all(np.isclose(weights,[0.  , 0.  , 0.  , 0.  , 0.25, 0.25, 0.  , 0.25, 0.25],atol=1e-4)))

  # xparticle = 0.0,0.5
  xparticle = np.array([0.0,0.5])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P**dim)
  assert(np.isclose(np.sum(weights),1))
  assert(np.all(grid_indicies == np.array([[31,15],[31,0],[31,1],[0,15],[0,0],[0,1],[1,15],[1,0],[1,1]])))
  assert(np.all(np.isclose(weights,[0.    , 0.0625, 0.0625, 0.    , 0.375 , 0.375 , 0.    , 0.0625,0.0625],atol=1e-4)))

  # xparticle = 0.5,0.0
  xparticle = np.array([0.5,0.0])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P**dim)
  assert(np.isclose(np.sum(weights),1))
  assert(np.all(grid_indicies == np.array([[31,15],[31,0],[31,1],[0,15],[0,0],[0,1],[1,15],[1,0],[1,1]])))
  assert(np.all(np.isclose(weights,[0.    , 0.    , 0.    , 0.0625, 0.375 , 0.0625, 0.0625, 0.375 ,0.0625],atol=1e-4)))

  # xparticle = 15,15
  xparticle = np.array([15,15])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P**dim)
  assert(np.isclose(np.sum(weights),1))
  assert(np.all(grid_indicies == np.array([[14,14],[14,15],[14,0],[15,14],[15,15],[15,0],[16,14],[16,15],[16,0]])))
  assert(np.all(np.isclose(weights,[0.015625, 0.09375 , 0.015625, 0.09375 , 0.5625  , 0.09375 , 0.015625, 0.09375 , 0.015625], atol=1e-4)))

  # xparticle = 15.5,15.5
  xparticle = np.array([15.5,15.5])
  grid_indicies, weights = fk.hockney_eastwood_function_jit(xparticle,ngridpts,P=P)
  assert(len(grid_indicies) == len(weights) == P**dim)
  assert(np.isclose(np.sum(weights),1))
  assert(np.all(grid_indicies == np.array([[15,15],[15,0],[15,1],[16,15],[16,0],[16,1],[17,15],[17,0],[17,1]])))
  assert(np.all(np.isclose(weights,[0.25, 0.25, 0.  , 0.25, 0.25, 0.  , 0.  , 0.  , 0.  ],atol=1e-4)))


if __name__ == "__main__":
  test_1d_P2()
  test_1d_P3()
  test_1d_P4()

  test_2d_P2()
  test_2d_P3()
