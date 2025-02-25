'''
Functions to read/write fields to files
'''


import numpy as np
import logging
import sys
import os
import h5py

from .field import *
from glob import glob


def read_from_file(filename): 
  """ Reads from textfile and output a list of Field objects.
  Adapted from iotools.py.
  
  Args:
      filename: String for name of file to be read through this function.
          
  Returns:
      field_list: A list of Fields objects.
  """
 
  # Check whether the file is gzipped and handle that seamlessly
  f = open(filename,'r')
  # Now parse the header to obtain the grid dimension and other format information
  version=int(f.readline().strip().split()[3])
  nfields=int(f.readline().strip().split()[3])

  dim=int(f.readline().strip().split()[3])
  npw=f.readline().strip().split()[4:4+dim]
  npw=[int(i) for i in npw] # Convert all list entries to int
  M = np.prod(npw)
  kspacedata=True
  complexdata=True
  
  flagsline=f.readline().strip().split()
  if flagsline[4] == "0":
    kspacedata=False
  if flagsline[9] == "0":
    complexdata=False
    
  f.readline() # Skip the last header line
  logging.info("  * Format version = {}".format(version))
  logging.info("  * Number of fields = {}".format(nfields))
  logging.info("  * Number of spatial dimensions = {}".format(dim))
  if dim == 3:
    logging.info("  * Grid dimension = {}x{}x{}".format(*npw))
  elif dim == 2:
    logging.info("  * Grid dimension = {}x{}".format(*npw))
  logging.info("  * K-space data? ",kspacedata)
  logging.info("  * Complex data? ",complexdata)
  if kspacedata:
    sys.stderr.write("\nError: k-space data is not supported\n")
    sys.exit(1)

  # Get the grid coordinates and the fields
  AllData = np.loadtxt(f)
  coords_flat = AllData[:,:dim]
  fields_flat = AllData[:,dim:]
  # note data is stored along COLUMNS 
  
  # if complex, combine columns into complex-valued data
  if complexdata:
    assert(2*nfields == fields_flat.shape[1])
    tmp = np.zeros((M,nfields),dtype=complex)
    for ifield in range(nfields):
      tmp[:,ifield].real = fields_flat[:,2*ifield]
      tmp[:,ifield].imag = fields_flat[:,2*ifield + 1]
    fields_flat = tmp 
    
  
  # calculate h (cell tensor) from coords
  coords = np.reshape(coords_flat ,list(npw) + [dim]) # reshape coords for easier usage
  if dim == 1:
      hvoxel = np.array([coords[1]])
  elif dim == 2:
      hvoxel = np.array([coords[1,0],coords[0,1]])
  elif dim == 3:
      hvoxel = np.array([coords[1,0,0],coords[0,1,0],coords[0,0,1]])
  h = hvoxel * npw

  # now create field objects in a list
  field_list = []
  for ifield in range(nfields):
    # reshape data before creating field
    data = np.reshape(fields_flat[:,ifield] , npw)
    field_list.append(Field(h=h,npw=npw,data=data))

  return field_list

def write_to_file(filename, fields): 
  """ Creates a text file based on the fields variable.
  Adapted from iotools.py.
  
  Args:
      filename: String for name of file to be written through this function.
      fields: A list of Fields objects.
          
  Returns:
      A text file based on fields.
  
  Raises:
      TypeError: if fields is not a list, than make it a one element list.
  """
  
  try:
    nfields = len(fields)
  except TypeError:
    fields = [fields]
    nfields = len(fields)

  print(f"Writing {nfields} fields to {filename}")

  # confirm that all fields have same dimension and coords
  for i in range(nfields):
    if (i !=0 ):
      assert (fields[i].dim == fields[i-1].dim)
      assert (fields[i].npw == fields[i-1].npw)
      assert (np.all(fields[i].h == fields[i-1].h))
      assert (np.all(fields[i].is_real_space == fields[i-1].is_real_space))
      assert (np.all(fields[i].data.dtype == fields[i-1].data.dtype))
      assert (np.all(fields[i].coords == fields[i-1].coords))

  # set values
  dim = fields[0].dim 
  npw = fields[0].npw
  h = fields[0].h
  coords = fields[0].coords
  assert(fields[0].is_real_space), "Can only currently write real space fields"
  is_real_space = True
  if fields[0].data.dtype == complex:
    complexdata = True
  else:
    complexdata = False
  kspacedata = (not is_real_space)
  M = np.prod(npw)


  # begin writing
  fout = open(filename,'w')
  fout.write("# Format version 3\n")
  fout.write("# nfields = {0}\n".format(nfields))
  fout.write("# NDim = {0}\n".format(dim))
  fout.write("# PW grid = ")
  for i in range(dim):
      fout.write("{0} ".format(npw[i]))
  fout.write("\n# k-space data = {0} , complex data = {1}\n".format(int(kspacedata),int(complexdata)))

  # Output coords and fields
  if dim == 1:
    fout.write("# Columns: x fielddata\n")
    for ix in range(0,npw[0]):
      fout.write("%12.10g " % coords[ix][0])
      for n in range(nfields):
        if complexdata:
          fout.write("%20.10g " % fields[n].data[ix].real)
          fout.write("%20.10g " % fields[n].data[ix].imag)
        else:
          fout.write("%20.10g " % fields[n].data[ix])
      fout.write("\n")

  elif dim == 2:
    fout.write("# Columns: x y fielddata\n")
    for ix in range(0,npw[0]):
      for iy in range(0,npw[1]):
        fout.write("%12.10g %12.10g " % (coords[ix,iy,0], coords[ix,iy,1]))
        for n in range(nfields):
          if complexdata:
            fout.write("%20.10g " % (fields[n].data[ix,iy].real))
            fout.write("%20.10g " % (fields[n].data[ix,iy].imag))
          else:
            fout.write("%20.10g " % (fields[n].data[ix,iy]))
        fout.write("\n")
      fout.write("\n")

  elif dim == 3:
    fout.write("# Columns: x y z fielddata\n")
    for ix in range(0,npw[0]):
      for iy in range(0,npw[1]):
        for iz in range(0,npw[2]):
          fout.write("%12.10g %12.10g %12.10g" % (coords[ix,iy,iz,0], coords[ix,iy,iz,1],coords[ix,iy,iz,2]))
          for n in range(nfields):
            if complexdata:
              fout.write("%20.10g " % (fields[n].data[ix,iy,iz].real))
              fout.write("%20.10g " % (fields[n].data[ix,iy,iz].imag))
            else:
              fout.write("%20.10g " % (fields[n].data[ix,iy,iz]))
          fout.write("\n")
        fout.write("\n")

  fout.close()
  
def get_PolyFTS_to_VTK_IdxMap(M,Nx):
    idx=np.zeros(M,dtype=np.uint64)
    ndim = len(Nx)
    if ndim == 1:
        for ix in range(Nx[0]):
            idx[ix] = ix
    elif ndim == 2:
        #looks good
        m=0
        for iy in range(Nx[1]):
          for ix in range(Nx[0]):
            idx[m] = ix*Nx[1] + iy
            m+=1
    elif ndim == 3:
        m=0
        for iz in range(Nx[2]):
          for iy in range(Nx[1]):
            for ix in range(Nx[0]):
                idx[m] = ix*Nx[1]*Nx[2] + iy*Nx[2] + iz
                #idx[m] = iz*Nx[0]*Nx[1] + iy*Nx[0] + ix
                m+=1
    return idx

  
def write_to_VTK(filename, fields):
    """ Creates a VTK file based on the fields variable.
    Adapted from FTS-tools/plot/PolyFTS_to_VTK.py.

    Args:
        filename: String for name of file to be written through this function.
        fields: A list of Fields objects.
  
    Returns:
        A VTK file based on fields.
    
    """
   
    #Check if vtk file is already in current directory
    #if os.path.isfile(filename):
    #    print("File '{}' already exist, file will be overwritten".format(filename))
    
    
    #Check if fields in Field object are compatible for VTK
    for i in range(1,len(fields)):
        assert(fields[i-1].npw == fields[i].npw)
        assert((fields[i-1].h == fields[i].h).all())

    #parameters
    dim = fields[0].dim
    npw = fields[0].npw
    orthorhombic = True
    binary = False
    M = fields[0].npw_total
    AllCoords = fields[0].coords
    
    #Create AllFields 
    fielddata = []
    for f in range(len(fields)):
        data = fields[f].data.flatten().reshape(fields[0].npw_total, 1).real
        fielddata.append(data)
    
    AllFields = np.concatenate(fielddata, axis=1)
    
    nfields = AllFields.shape[-1]
   
    #determine spacing
    if orthorhombic == True:
        if dim == 1:
            spacing = (AllCoords[1][0]) #2d array (x, field_index)
        if dim == 2:
            spacing = (AllCoords[1,0][0], AllCoords[0,1][1]) #3d array (x,y, field_index)
        if dim == 3:
            #spacing = (AllCoords[1,0,0][2], AllCoords[0,1,0][1], AllCoords[0,0,1][0])
            spacing = (max(AllCoords[1,0,0]), max(AllCoords[0,1,0]), max(AllCoords[0,0,1])) #4d
            if np.any(np.isclose(spacing,0.0)):
              raise RuntimeError (f"Error! One of the grid spacings is close to zero: {spacing}")
            #spacing = (AllCoords[1,0,0][2], AllCoords[0,1,0][1], AllCoords[0,0,1][0])
        logging.info("Mesh spacing       {}".format(spacing))
    
    #Manipulating AllFields and AllCoords 
    AllCoords = np.ravel(AllCoords)
    AllCoords = np.reshape(AllCoords,(M,dim ))
    AllFields = np.ravel(AllFields)
    AllFields = np.reshape(AllFields,(M,nfields))

    #Generate mapping from PolyFTS order to VTK order 
    IdxMap = get_PolyFTS_to_VTK_IdxMap(M,npw)
    # Remap field samples from PolyFTS order to VTK order 
    AllCoords = AllCoords[IdxMap,:]
    AllFields = AllFields[IdxMap,:]

    AllCoords = AllCoords.T 
    AllFields = AllFields.T

    #Write VTK file
    o = open(filename,"w")
    
    if orthorhombic:
        o.write("# vtk DataFile Version 3.0\n")
        o.write("PolyFTS field data\n")
        o.write("ASCII\n")
        o.write("DATASET STRUCTURED_POINTS\n")
        if dim == 1:
            o.write("DIMENSIONS {} 1 1\n".format(*npw))
            o.write("ORIGIN 0\n")
            o.write("SPACING {} 0 0\n".format(*spacing))
        elif dim == 2:
            o.write("DIMENSIONS {} {} 1\n".format(*npw))
            o.write("ORIGIN 0 0 0\n")
            o.write("SPACING {} {} 0\n".format(*spacing))
        elif dim == 3:
            o.write("DIMENSIONS {} {} {}\n".format(*npw))
            o.write("ORIGIN 0 0 0\n")
            o.write("SPACING {} {} {}\n".format(*spacing))
        o.write("POINT_DATA {0}\n".format(M))
        for i in range(nfields):
            o.write("SCALARS field{0} float 1\n".format(i))
            o.write("LOOKUP_TABLE default\n")
            o.close();o = open(filename,"ab")
            np.savetxt(o, AllFields[i], fmt="%.11f")
            o.close();o = open(filename,"a")
    else:
      o.write("# vtk DataFile Version 3.0\n")
      o.write("PolyFTS field data\n")
      if binary:
          o.write("BINARY\n")
      else:
          o.write("ASCII\n")
      o.write("DATASET STRUCTURED_GRID\n")
      if dim == 1:
          o.write("DIMENSIONS {} 1 1\n".format(*npw))
      elif dim == 2:
          o.write("DIMENSIONS {} {} 1\n".format(*npw))
      elif dim == 3:
          o.write("DIMENSIONS {} {} {}\n".format(*npw))
      o.write("POINTS {} double\n".format(M))
      # Remap the mesh coordinates to VTK order (Change in the future)
      AllCoords = AllCoords[:,IdxMap]
      #np.savetxt('coords.dat',AllCoords.transpose()) # Debug
      o.close(); o = open(filename,'ab')
      if binary:
          AllCoords.transpose().tofile(o)
      else:
          np.savetxt(o, AllCoords.transpose(), fmt="%0.11f")
      o.close(); o = open(filename,'a')

      o.write("\nPOINT_DATA {0}\n".format(M))
      for i in range(nfields):
      #for i in range(1):
          # write as scalar
          #o.write("SCALARS field{0} float 1\n".format(i)) #STARTED TO SEG FAULT WHEN TRYING TO READ SCALARS
          #o.write("LOOKUP_TABLE default\n")
          #np.savetxt(o, AllFields[i], fmt="%14.11f")
          #np.savetxt(o, np.vstack([AllCoords,AllFields[i]]).transpose(), fmt="%14.11f")

          # write as vector
          o.write("VECTORS field{0} double\n".format(i)) #writing as vector fixed the seg fault
          tmp=np.zeros((3,M))
          tmp[0,:] = AllFields[i]
          o.close(); o = open(filename,'ab')
          if binary:
              tmp.transpose().tofile(o)
          else:
              np.savetxt(o, tmp.transpose(), fmt="%14.11f")
          o.close(); o = open(filename,'a')
    o.close()
    
def write_to_HDF5(filename, fields):
  """ Creates a HDF5 file based on the fields variable.
  
  Args:
      filename: String for name of file to be written through this function.
      fields: A list of Fields objects.
          
  Returns:
      A HDF5 file based on fields.
  
  Raises:
      TypeError: if fields is not a list, than make it a one element list.
  """
  # check if file is HDF5
  if not filename.endswith(".h5"):
    print("Must be an HDF5 file")
    sys.exit(1)

  try:
    nfields = len(fields)
  except TypeError:
    fields = [fields]
    nfields = len(fields)

  print(f"Writing {nfields} fields to {filename}")

  # confirm that all fields have same dimension and coords
  for i in range(nfields):
    if (i !=0 ):
      assert (fields[i].dim == fields[i-1].dim)
      assert (fields[i].npw == fields[i-1].npw)
      assert (np.all(fields[i].h == fields[i-1].h))
      assert (np.all(fields[i].is_real_space == fields[i-1].is_real_space))
      assert (np.all(fields[i].data.dtype == fields[i-1].data.dtype))
      assert (np.all(fields[i].coords == fields[i-1].coords))

  # set values
  dim = fields[0].dim 
  npw = fields[0].npw
  h = fields[0].h
  coords = fields[0].coords
  assert(fields[0].is_real_space), "Can only currently write real space fields"
  is_real_space = True
  if fields[0].data.dtype in [np.complex64, np.complex128]:
    complexdata = True
  else:
    complexdata = False
  kspacedata = (not is_real_space)
  M = np.prod(npw)

  # begin writing
  fout = h5py.File(filename, 'w')
  group = fout.create_group("fields")

  # write attributes
  group.attrs["cell_tensor"] = h
  group.attrs.create("complexdata", complexdata, dtype="u1")
  group.attrs.create("dim", dim, dtype=np.intc)
  group.attrs.create("kspacedata", kspacedata, dtype="u1")
  group.attrs["npw"] = npw

  # write datasets
  for i in range(nfields):
    if complexdata:
        data_to_write = fields[i].data
    else:
        data_to_write = np.empty(npw, dtype=complex)
        data_to_write.real = fields[i].data
    dset = group.create_dataset(f"field{i}", npw, data=data_to_write)

  # close file
  fout.close()

def read_from_HDF5(filename):
  """ Reads from HDF5 file and output a list of Field objects.
  
  Args:
      filename: String for name of file to be read through this function.
          
  Returns:
      field_list: A list of Fields objects.
  """
  # check if file is HDF5
  if not filename.endswith(".h5") and not filename.endswith(".hdf5"):
    print("Must be an HDF5 file")
    sys.exit(1)

  # open HDF5 file and group
  f = h5py.File(filename, 'r')
  group = f["fields"]

  # get values from group attributes
  nfields = len(group)
  h = group.attrs["cell_tensor"]
  complexdata = group.attrs["complexdata"]
  dim = group.attrs["dim"]
  kspacedata = group.attrs["kspacedata"]
  npw = group.attrs["npw"]

  # log attributes
  logging.info("  * Number of fields = {}".format(nfields))
  logging.info("  * Number of spatial dimensions = {}".format(dim))
  if dim == 3:
    logging.info("  * Grid dimension = {}x{}x{}".format(*npw))
  elif dim == 2:
    logging.info("  * Grid dimension = {}x{}".format(*npw))
  logging.info("  * K-space data? ",kspacedata)
  logging.info("  * Complex data? ",complexdata)
  if kspacedata:
    sys.stderr.write("\nError: k-space data is not supported\n")
    sys.exit(1)

  # read data to field_list
  field_list = []
  for i in range(nfields):
    data_to_read = group[f"field{i}"][:]
    if not complexdata:
        data_to_read = data_to_read.real
    field_list.append(Field(h=h, npw=npw, data=data_to_read))

  # close file
  f.close()

  return field_list

def write_to_cube_files(prefix, fields):
    """ Creates a cube file for each Field object in the fields list. Because 
        Field objects lack atoms, each voxel is considered an atom with its
        charge set to its value in the Field object.
    
    Args:
        prefix: String for the prefix of files written through this function.
        fields: A list of Field objects.
  
    Returns:
        A cube file for each Field object in the fields list.
    
    """

    # Check if the fields are compatible with the cube format.

    for i in range(len(fields)):
        assert(fields[i].dim == 3), "Field objects are not 3-dimensional"
        assert(fields[i].is_orthorhombic()), "Field objects are not orthorhombic"
    for i in range(1,len(fields)):
        assert(fields[i-1].npw == fields[i].npw), "Field objects have inconsistent npw"
        assert((fields[i-1].h == fields[i].h).all()), "Field objects have inconsistent cell tensors"

    # Get parameters from the first Field object.

    npw_x = fields[0].npw[0]
    npw_y = fields[0].npw[1]
    npw_z = fields[0].npw[2]
    nvoxels = npw_x * npw_y * npw_z
    xvoxel_len = fields[0].h[0][0] / npw_x
    yvoxel_len = fields[0].h[1][1] / npw_y
    zvoxel_len = fields[0].h[2][2] / npw_z

    # Loop through each Field object in the fields list.

    for i in range(len(fields)):
        data = fields[i].data.real
        coords = fields[i].coords

        # Write the cube file for the given Field object. This follows 
        # https://paulbourke.net/dataformats/cube for formatting.

        with open(f"{prefix}_field{i}.cube", "w") as o:

            # Write the header.

            o.write(f"GAUSSIAN CUBE FILE WRITTEN BY FIELDKIT\n")
            o.write(f"OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n")
            o.write(f"{nvoxels:>8}    {0:.6f}    {0:.6f}    {0:.6f}\n")
            o.write(f"{npw_x:>8}    {xvoxel_len:.6f}    {0:.6f}    {0:.6f}\n")
            o.write(f"{npw_y:>8}    {0:.6f}    {yvoxel_len:.6f}    {0:.6f}\n")
            o.write(f"{npw_z:>8}    {0:.6f}    {0:.6f}    {zvoxel_len:.6f}\n")
            for ix in range(npw_x):
                for iy in range(npw_y):
                    for iz in range(npw_z):
                        x = coords[ix][iy][iz][0]
                        y = coords[ix][iy][iz][1]
                        z = coords[ix][iy][iz][2]
                        value = data[ix][iy][iz]
                        o.write(f"{1:>8}    {value:.6f}    {x:.6f}    {y:.6f}    {z:.6f}\n")

            # Write the volumetric data.

            for ix in range(npw_x):
                for iy in range(npw_y):
                    for iz in range(npw_z):
                        o.write(f"{data[ix][iy][iz]:.5E} ")
                        if iz % 6 == 5:
                            o.write("\n")
                    o.write("\n")

def read_from_cube_files(prefix):
    """ Reads cube file(s) and converts them to a list of Field objects.

    Args:
        prefix: String for the prefix of the cube file(s) read by this function.

    Returns:
        A list of Field objects.

    """

    # Loop through each cube file.

    field_list = []
    filepaths = sorted(glob(f"{prefix}_field*.cube"))
    for filepath in filepaths:
        with open(filepath, "r") as o:
            h = np.zeros((3, 3))
            npw = []
            values = []
            for i, line in enumerate(o):
                if i <= 2:  # Skip the first 3 lines.
                    continue
                tokens = line.split()
                if 3 <= i <= 5:  # Extract the npw and cell tensor.
                    npw_j = int(tokens[0])
                    npw.append(npw_j)
                    h[i-3][i-3] = float(tokens[i-2]) * npw_j
                else:
                    if len(tokens) == 6:  # Stop when volumetric data reached.
                        break
                    value = float(tokens[1])
                    values.append(value)

        # Rehape the values given the npw.

        values = np.asarray(values)
        data = np.reshape(values, npw)

        # Create the Field object and append it to the field list.

        field_list.append(Field(h=h,npw=npw,data=data))

    # Return the list of Field objects.

    return field_list

