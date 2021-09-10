
from .field import *

from skimage import measure

# Need to increase recursion limit for burning algorithm
import sys
sys.setrecursionlimit(800000)

def calc_domain_stats(field, density_threshold, plotMesh=False,outputMesh=False,add_periodic_domains=False, applyPBC=True):
    ''' calculate properties of domains using a mesh
        Adapted from domaintools.py
    '''
    
    dim = field.dim
    h = field.h
    
    ndomains, domainID, image_flags, domainBorder = identify_discrete_domains(field, density_threshold)

    com = np.zeros((ndomains, dim))
    surface_area = np.zeros(ndomains)
    volume = np.zeros(ndomains)
    volume_nomesh = np.zeros(ndomains)
    IQ = np.zeros(ndomains)

    #for each domain
    for idomain in range(0,ndomains):
        # calc center of domain
        com[idomain,:] = _calc_domain_center(idomain+1, field, domainID, image_flags, units='coord')
        

        if dim == 2:
            # mesh domain
            contours,density_centered = _mesh_single_domain(field, idomain+1, density_threshold, domainID, image_flags,domainBorder, wrap_before_mesh=applyPBC)
            assert (len(contours) == 1), "The contour should only be one curve, if not the area and volume calculations will be completely wrong!"

            # get surface area (perimeter) and volume (area)
            surface_area[idomain] = _contour_perimeter(contours[0])
            volume[idomain] = _contour_area(contours[0])

            if plotMesh: 
                # draw surface behind the mesh
                _plot_contours_2D(field,contours,filename="mesh.{}.png".format(idomain+1),surface=density_centered)
        if dim == 3: 
            # mesh domain
            verts, faces, normals, values,density_centered = _mesh_single_domain(field, idomain+1, density_threshold, domainID, image_flags,domainBorder,wrap_before_mesh=applyPBC)

            # get surface area, volume and isoperimetric quotient
            surface_area[idomain] = measure.mesh_surface_area(verts, faces)
            volume[idomain] = _mesh_volume(verts,faces)
            if plotMesh: 
                _plot_mesh_3D(field, verts,faces, filename="mesh.{}.png".format(idomain+1))
            if outputMesh:
                _write_mesh(verts,faces,fileprefix="mesh.{}.".format(idomain+1))

        IQ[idomain] = _calc_IQ(dim,surface_area[idomain], volume[idomain])

        # should work for 2d and 3d
        volume_nomesh[idomain] = _voxel_volume(field, idomain+1, domainID) # get volume from voxels
        
        # sanity check to make sure that volume is reasonable
        if not np.isclose(volume_nomesh[idomain], volume[idomain], rtol=0.2):
            print(f"WARNING: the volume computed via mesh is >20% different from volume from voxels {volume[idomain]} != {volume_nomesh[idomain]}. Check the mesh to make sure the domain mesh is fully closed.")

    if add_periodic_domains:
        for idomain in range(1,ndomains+1):  
            extracom = _pbc_domain_locs(idomain,com[idomain-1])
            if extracom:
                com = np.concatenate((com,extracom))
                extra_num = len(extracom)
                IQ = np.concatenate((IQ,[IQ[idomain-1]]*extra_num)) 
                surface_area = np.concatenate((surface_area,[surface_area[idomain-1]]*extra_num)) 
                volume = np.concatenate((volume,[volume[idomain-1]]*extra_num)) 


    stats = {'ndomains': ndomains, 'center': com, 'surface_area': surface_area, 'volume': volume, 'volume_nomesh':volume_nomesh, 'IQ':IQ}

    return stats

def _calc_domain_center(idomain, field, domainID, image_flags, units='box'):
    ''' given a domain index, apply PBC and return the center of mass
        Can return result in 'box' units (0 to Nx) or in 'coord' units (0 to boxl)
    '''
    h = field.h
    

    isdomain = (domainID == idomain)
    N = np.sum(isdomain)
    indicies = np.transpose(np.nonzero(isdomain))
    
    dim = len(domainID.shape) 
    Nx = domainID.shape
    coords = np.zeros((N,dim))
  
    #TODO could I do this without for loop? (will be faster)
    for i in range(N):
        index = tuple(indicies[i])
  
        if units == "box":
           coord = index + image_flags[index] * Nx
        elif units == "coord":
           # this was orig (for orthorhombic boxes
           #coord = self.__coords[index] + self.__image_flags[index] * self.__boxl
  
           # new for non-orthorhombic boxes (check that this works for orthorhombic though!)
           #shift = np.array(np.mat(self.__hvoxel.T) * np.mat(image_flags[index] * Nx).T).T
           #shift = np.array(np.mat(h.T) * np.mat(image_flags[index]).T).T # depreciated np.mat
           shift = np.dot(h, image_flags[index])
           coord = field.coords[index] + shift
        else:
            raise ValueError("Invalid units entry of \'%s\'" % units)
        
        coords[i] = coord
    
    # now average in order to get center of the domain (each point weighted evenly)
    return np.average(coords,axis=0)
  
def _mesh_single_domain(field, idomain, density_threshold, domainID, image_flags,domainBorder, wrap_before_mesh=True):
    '''
    Function to:
    1) apply PBC to the domains so that an entire domain is continuous (ie not split across boundaries)
    2) Grab a little margin around each domain (the domain's "border") so that marching cubes can interpolate. The border is computed in identifyAndIndexDomains().
    3) Mesh the domain using marching cubes
    '''

    Nx = field.npw_Nd
    h = field.h
    dim = field.dim

    hvoxel = field.hvoxel() # TODO: remove all uses of hvoxel
    
    isdomain = (domainID == idomain)
    #isborder = (self.__borderID == idomain)
    isborder = np.zeros(Nx,dtype=np.bool)
    # convert to tuple to correctly set indicies of isborder
    isborder[tuple(domainBorder[idomain-1])] = True
    
    alldensity = field.data

    # center box and properties around center of mass (so that domains don't cross pbc)
    # np.roll is the key function here
    # if domains percolate then this will break 
    com_box = _calc_domain_center(idomain,field, domainID, image_flags, units='box')
    #com_coord = self.calcDomainCOM(idomain,units='coord')
    #coords_tmp = np.copy(self.__coords)
    for i in range(dim):
        shift = int(0.5*Nx[i] - com_box[i])
        isdomain = np.roll(isdomain,shift,axis=i)
        isborder = np.roll(isborder,shift,axis=i)
        #coords_tmp = np.roll(coords_tmp,shift,axis=i)
        alldensity = np.roll(alldensity,shift,axis=i)

           
    # isolate the domain of interest
    isdomain_or_isborder = isdomain + isborder # since both bool, sum is the union of the two fields
    mydensity = np.zeros(Nx)
    mydensity[isdomain_or_isborder] = alldensity[isdomain_or_isborder]
    
    # mesh! (using scikit-image)
    if dim == 2:
        # calculate contours in 'box' units
        contours = measure.find_contours(mydensity, density_threshold) 

        # convert 'box' units to 'coords' units (this is key for non-orthorhombic cells)
        for i,c in enumerate(contours):
            #contours[i] = np.array((np.mat(hvoxel).T * np.mat(c).T).T) #old depreciated use of np.mat
            contours[i] = np.dot(c, hvoxel)

        return contours, alldensity
    elif dim == 3:
        #from skimage import measure

        #verts, faces, normals, values = measure.marching_cubes_lewiner(mydensity, self.__density_threshold, spacing = self.__gridspacing)
        # do not use spacing=self.__gridspacing, let marching cubes calculate verticies in 'box' units (0,Nx) 
        verts, faces, normals, values = measure.marching_cubes(mydensity, density_threshold)

        # convert 'box' units to 'coords' units (this is key for non-orthorhombic cells)
        for i,v in enumerate(verts):
            #verts[i] = np.array((np.mat(hvoxel).T * np.mat(v).T).T) # depreciated np.mat
            verts[i] = np.dot(hvoxel,v)
            n = normals[i]
            #normals[i] = np.array((np.mat(hvoxel).T * np.mat(n).T).T) # depreciated np.mat
            normals[i] = np.dot(hvoxel,n) 


        return verts, faces, normals, values, alldensity
    else:
        raise ValueError("Meshing makes no sense in 1 dimension!")

def _contour_perimeter(contour):
    '''calculate perimeter of contour by suming up the line-segment lengths
    '''
    assert (np.all(contour[0] == contour[-1])), "Contour must be closed! (1st point == last point)"

    #TODO vectorize this for loop
    p = 0.0
    n=contour.shape[0]
    for i in range(n-1):
       v = contour[i+1] - contour[i] 
       p += np.sqrt(np.square(v).sum())
    return p

def _contour_area(contour):
    ''' Calculate area of shape enclosed in contour
        similar to calculating mesh volume
        use trick from http://geomalgorithms.com/a01-_area.html
    '''
    assert (np.all(contour[0] == contour[-1])), "Contour must be closed! (1st point == last point)"
    
    #TODO vectorize this for loop
    area = 0.0
    n=contour.shape[0]
    for i in range(n-1):
        area += np.cross(contour[i],contour[i+1])
    return 0.5*np.abs(area)


def _mesh_volume(verts, faces):
    '''calculate volume of a mesh, using cross product trick
    '''
    actual_verts = verts[faces]
    v0 = actual_verts[:,0,:]
    v1 = actual_verts[:,1,:]
    v2 = actual_verts[:,2,:]
   
    # TODO: dont do the volume rescaling here, instead change the actual position of "verts" in getDomainStats my scaling each vert position by h (or something along these lines)

    # introduce factor to scale the volume if non-orthorhombic box
    # this is because the mesh is generated assuming a
    #if self.__orthorhombic:
    #    factor=1.0
    #else:
    #    factor = self.__volvoxel / np.prod(self.__gridspacing)
    factor=1.0 # FIXME, assumes cell is orthorhombic. Uncomment lines above to fix

    # 1/6 \sum v0 \cdot (v1 x v2)
    return factor * 1.0/6.0 * np.abs( (v0*np.cross(v1,v2)).sum(axis=1).sum() )

def _voxel_volume(field, idomain, regionID):
    ''' Get volume of idomain using voxels
    '''
    #v_voxel = np.prod(self.__gridspacing) # volume of single voxel
    v_voxel = np.linalg.det(field.h) / field.npw_total
    #v_voxel = self.__volvoxel
    n_voxel = np.sum(regionID == idomain) # number of voxels in ith domain
    return v_voxel*n_voxel

def _calc_IQ(dim, area, vol):
   '''returns isoperimetric coefficient. 1 for perfect circle or sphere, less for other shapes
      note that in 2d "area" is actually perimeter, and "vol" is actually area
      This difference didn't seem to warrant a completely different method though
   '''
   if dim == 2:
       return 4.0*np.pi*vol / (area * area)
   elif dim == 3:
       return 36.0*np.pi * vol*vol / (area * area * area)

def _plot_contours_2D(field, contours, surface=None, filename=None):
    ''' Plot a mesh from marching squares
    '''
    import matplotlib.pyplot as plt
    Nx = field.npw_Nd
    hvoxel = field.hvoxel()

    # Display the image and plot all contours found
    fig, ax = plt.subplots()

    ax.set_aspect(1)

    if surface is not None:
        x = np.arange(Nx[0])
        y = np.arange(Nx[1])
        xx,yy = np.meshgrid(x,y)
        
        # nice one-liner to rotate all of xx and yy using hvoxel
        xxrot,yyrot = np.einsum('ji, mni -> jmn', hvoxel.T, np.dstack([xx, yy]))
        
        # using pcolormesh allows us to use non-orthorhombic boxes 
        im=ax.pcolormesh(xxrot,yyrot,surface.T,shading='auto')
        fig.colorbar(im,ax=ax)
        
        # imshow only worked for orthorhombic boxes
        #ax.imshow(surface.T, interpolation='nearest')

    for n, contour in enumerate(contours):
        ax.plot(contour[:, 0], contour[:, 1], linewidth=2, color='k',ls='--',marker='o')

    #ax.axis('image')
    #ax.set_xticks([])
    #ax.set_yticks([])
    if not filename:
        plt.show()
    else:
        plt.savefig(filename)
    plt.close()

def _plot_mesh_3D(field, verts, faces, filename=None):
    ''' Plot a mesh from marching cubes
    '''
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    assert(field.is_orthorhombic())
    assert(field.dim == 3)
    boxl = np.diag(field.h)

    # Display resulting triangular mesh using Matplotlib. This can also be done
    # with mayavi (see skimage.measure.marching_cubes_lewiner docstring).
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Fancy indexing: `verts[faces]` to generate a collection of triangles
    mesh = Poly3DCollection(verts[faces])
    mesh.set_edgecolor('k')
    ax.add_collection3d(mesh)

    ax.set_xlim(0, boxl[0])  
    ax.set_ylim(0, boxl[1])  
    ax.set_zlim(0, boxl[2])  

    plt.tight_layout()
    if not filename:
        plt.show()
    else:
        plt.savefig(filename)

    plt.close()

def _write_Mesh(verts,faces,fileprefix="mesh."):
   '''save mesh to a file'''
   np.savetxt(fileprefix + "verts.dat",verts,header='Autogenerated mesh file. Contains x y z positions of each vertex' )
   np.savetxt(fileprefix + "faces.dat",faces, header='Autogenerated mesh file. Contains vertex indicies of each triangle in mesh')
    


 

def identify_discrete_domains(field, density_threshold):
    ''' Identify all discrete domains present within a field given a density threshold

        Uses the "burning algorithm"
        Adapted from domaintools.py

        Also sets - image_flags (which PBC a domain belongs to) and 
                  - isborder (whether a grid is adjacent to a domain)
        Though these are not currently used outside of this function

    Args:
        field: field to identify domains of
        density_threshold: threshold used to define domains

    Returns:
        ndomains: the number of discrete domains 
        domainID: a ndarray of shape equal to field.data containing the index of that domain. If domainID == 0, it is the continuous domain. Points with domainID == i, correspond to the ith domain
        imageflags: imageflags in each dimension
    '''

    isdomain_array = field.data > density_threshold 

    Nx = field.npw_Nd
    dim = field.dim

    # if domainID == -1, it has not been visited
    domainID = np.full(Nx,-1, dtype=np.int32)
    # image_flags are only for the domains themselves, the image flags of the border are not needed
    image_flags = np.zeros(list(Nx) + [dim])
    
    domainBorder = [[]]

    domain_number = 1;

    #this is where the recursive magic happens
    for i in np.ndindex(Nx):
      if (domainID[i] == -1):
        if (isdomain_array[i]):
          current_image_flag = np.zeros(dim)
          _spread_domain(i, domain_number, isdomain_array,current_image_flag, domainID, image_flags, domainBorder, Nx);
          domainBorder.append([])
          domain_number += 1;
        else:
          # note - dont assign borders here, this is acomplished inside of spread_domain()
          domainID[i] = 0;
          image_flags[i]= np.zeros(dim)
    
    # now cleaning up
    ndomains = domain_number-1;
    
    # remove last element from lists (should be empty)
    assert (domainBorder[-1] == [])
    del domainBorder[-1] 
    
    # check that lengths of domain structs are correct
    assert (len(domainBorder) == ndomains)

    # convert border and imageflag lists to numpy arrays
    for i in range(ndomains):
        domainBorder[i] = np.array(domainBorder[i]).transpose()
    
    return ndomains, domainID, image_flags, domainBorder
    

def _spread_domain(coord_center, domain_number, isdomain_array,current_image_flag, domainID, image_flags, domainBorder, Nx):
    ''' recursive function:
        given a point, find the neighbors of that point, 
        for each neighbor, send back into function
        
        WARNING: domainID and image_flags must be passed by reference
    '''

    domainID[coord_center] = domain_number;
    image_flags[coord_center] = current_image_flag

    neighbors,neigh_image_flags = _get_neighbors(coord_center, current_image_flag, Nx);

    for i in range(len(neighbors)):
        neighbor = neighbors[i]
        image_flag = tuple(neigh_image_flags[i])
        if (domainID[neighbor] == -1):
            if (isdomain_array[neighbor]):
              _spread_domain(neighbor, domain_number, isdomain_array, image_flag, domainID, image_flags, domainBorder, Nx);
            else:
              domainID[neighbor] = 0;
              

        if domainID[neighbor] == 0:
          
          # only append to list if neighbor isn't in there already
          if neighbor not in domainBorder[domain_number-1]:
              # must have neighbors that are domain (since spread domain is only called 
              #   if coord_center is a domain). Therefore, it's a border
              domainBorder[domain_number-1].append(neighbor)

              # set image flags of non-domain adjacent to domain according to the domain
              # basically, I need the border to have the correct image flags
              # NOTE: image flags of borders aren't used anymore
              #self.__domainBorderImageFlags[domain_number-1].append(image_flag)

   
def _get_neighbors(coord_center,center_image_flag, Nx):
   ''' given a coord (tuple), return 
        1) the neighbors of that coord (also tuple) AND 
        2) the image_flag (which PBC) that neighbor corresponds to
   '''
   dim = len(coord_center)

   neighbors = [];
   neigh_image_flags = np.tile(center_image_flag, (2*dim,1))
   for i in range(dim):
      coord_neigh = np.copy(coord_center)
      coord_neigh[i] -= 1;
      _apply_pbc(coord_neigh, neigh_image_flags[2*i], Nx);
      neighbors.append(tuple(coord_neigh))

      coord_neigh = np.copy(coord_center)
      coord_neigh[i] += 1
      _apply_pbc(coord_neigh, neigh_image_flags[2*i+1], Nx)
      neighbors.append(tuple(coord_neigh))

   return neighbors, neigh_image_flags


def _apply_pbc(coord,image_flag, Nx):
  ''' apply periodic boundary conditions in index space. helper function for get_neighbors'''
  dim = len(coord)
  for i in range(dim):
     if coord[i] >= Nx[i]: 
         coord[i] = 0
         image_flag[i] += 1
     if coord[i] < 0:             
         coord[i] = Nx[i] - 1
         image_flag[i] -= 1


def _pbc_domain_locs(field, idomain,regionID, image_flags, local_com):
    '''This function returns the locations of the other domains on the periodic boundary.  
    for example for a domain with its center on the corner of the box, it would return all
    the other box corners'''
    
    dim = field.dim
    hvoxel = field.hvoxel() # TODO: remove all uses of hvoxel
    Nx = field.npw_Nd
      

    extra_com = []
    domain = (regionID == idomain)
    local_flags = image_flags[domain]

    # unique_flags contains a minimal list of the PBC that we want to add new domains to
    unique_flags = set([])
    for i in range(np.shape(local_flags)[0]):
        unique_flags.add(tuple(local_flags[i]))
    if dim == 2:
        unique_flags.remove((0,0))#remove duplicate com
    elif dim == 3:
        unique_flags.remove((0,0,0))#remove duplicate com


    #pdb.set_trace()
    for flag in unique_flags:
        flag = np.array(flag)

        # old - trenton 
        #new_com = -1*flag*self.__boxl+local_com
        #find the location of the extra periodic com by adding the box length times the flag to the current com 

        # new - without assuming orthorhombic box
        shift = np.array(hvoxel.T * (flag * Nx).T).T
        shift = np.reshape(shift,(dim,))
        new_com = local_com - shift #note minus shift to follow trentons convention

        extra_com.append(new_com)
        num_extra = len(extra_com)
    return extra_com

