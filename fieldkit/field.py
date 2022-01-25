import numpy as np

class Field:
    """ Class for Field object used to model 1D, 2D, 3D structures.
    
    Attributes:
        npw_Nd: number of grid points in each dimension
        npw_total: total number of grid points for all dimensions
        dim: dimension of field
        data: stored values of field at each grid point
        h: a cell tensor, in which the columns are box vectors that describe the coordinate system
        is_real_space: tells if Fourier transform has been used
        coords: stores x,y,z coordinates of each grid point
        
    """
    def __init__(self,npw_Nd=None, data=None, h=None):
        """ Defines default values for class attributes. """
        self.npw_Nd = None        
        self.npw_total = None     
        self.dim = len(npw_Nd)             
        self.data = None        
        self.h = None             
        self.is_real_space = True 
        self.coords = None       

        #manipulate attributes based on input values
        if np.all(npw_Nd != None):
            self.npw_Nd = tuple(npw_Nd)
            self.npw_total = np.prod(npw_Nd)
        
        # initialize data
        if np.all(data != None):
            self.set_data(data)
        elif np.all(npw_Nd != None): 
            self.data= np.zeros(npw_Nd) # default to zeros
        
        # initialize h
        if np.all(h != None):
            self.set_h(h)
        else: 
            self.set_h(np.eye(self.dim)) # default to identity

    def hvoxel(self):
      # TODO: remove all uses of hvoxel. I added this to support the domain analysis code
      if self.dim == 2:
        hvoxel = np.array([self.coords[1,0],self.coords[0,1]])
      elif self.dim == 3:
        hvoxel = np.array([self.coords[1,0,0],self.coords[0,1,0],self.coords[0,0,1]])
      return hvoxel
        
    def is_orthorhombic(self):
        # sanity checks
        #assert(self.h != None), "Calling is_orthorhombic but h is not initialized"
        assert(not np.all(self.h == 0)), f"Calling is_orthorhombic but h is all zeros. {h}"

        # check that upper and lower diag are zero
        diag = np.diag(np.diag(self.h)) # this is a 2d matrix with all offdiag = 0
        nondiag = self.h - diag # 2d matrix with diagonal elements = 0
         
        zero_nondiag = np.all(nondiag ==0.0)
        if zero_nondiag:
            return True 
        else:
            return False
        
    def volume(self):
      return np.linalg.det(self.h)

    def set_h(self,h):
        """ Checks that h is square and sets coords. """
        assert(self.dim*self.dim == h.size)
        for i in range(len(h.shape)):
            assert(h.shape[i] == self.dim)
        self.h = h 
        
        # if set_h is called, update coords 
        self.coords = self.CoordsFromH()
        #if not np.all(self.coords != None):
        #    self.coords = self.CoordsFromH()

    def set_data(self,data):
        """Sets data based on size."""
        if self.npw_total != None:
            # if already init with size, check that size matches
            assert(data.size == self.npw_total)
            assert(data.shape == self.npw_Nd)
        else:
            # if not yet init with size, then infer size from data
            self.npw_total = data.size
            self.npw_Nd = data.shape
            self.dim = len(data.shape)
            assert (self.dim >= 1 and self.dim <=3)
    
        # now set data
        self.data = data

   
    def CoordsFromH(self):
        """Compute coordinates of field from current "h" box vector."""
        coords = np.zeros(list(self.npw_Nd)+[self.dim])
        for idx in np.ndindex(self.npw_Nd):
            x = np.array(idx) / self.npw_Nd #0-1
            r = np.matmul(self.h,x)
            coords[idx] = r
        return coords;

    def set_coords(self, coords):
        """User can manually set coords."""
        self.coords = coords



