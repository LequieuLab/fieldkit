'''
Functions to plot fields 
'''

from .field import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot(fields,dpi=100,show=True,filename=None,imag=False):
  """Plot fields using matplotlib
    
    Args: 
      fields: a list of Field objects
      dpi: dpi (resolution) of specified image
      show: whether or not to show the plot
      filename: output filename for plotfile
      imag: whether to plot 'real' or 'imag' component of fields

    Returns:
      none

  """
  try:
    nfields = len(fields)
  except TypeError:
    fields = [fields]
    nfields = len(fields)
  

  if nfields == 1:
    nrows = 1; ncols = 1
  elif nfields == 2:
    nrows = 1; ncols = 2
  elif nfields == 3:
    nrows = 1; ncols = 3
  elif nfields == 4:
    nrows = 2; ncols = 2
  else:
    # try to factor
    #nrows = 6
    #while (nfields % nrows != 0): 
    #  nrows -= 1
    #ncols = int(nfields / nrows)

    # force square
    nrows = int(np.ceil(nfields **0.5))
    ncols = int(np.ceil(nfields **0.5))
    print(f"Automatically set {nrows = }, {ncols = }")
  #else:
  #  raise RuntimeError("nfields > 4 not currently supported")
 
  fig = plt.figure(figsize=(ncols*3.33,nrows*3.33),dpi=dpi)

  for i in np.ndindex(nrows,ncols):
    ifield = i[0]*ncols + i[1]
    if ifield >= nfields:
      continue
    field = fields[ifield]

    # create new axis for plot
    if fields[0].dim <= 2:
      ax = fig.add_subplot(nrows,ncols,ifield+1)
    else:
      ax = fig.add_subplot(nrows,ncols,ifield+1, projection='3d')

    # grab either the real or imag part of fields (depending on imag input arg)
    if field.is_complex():
      if imag == False:
        print(f"Note: not plotting imaginary part of field {ifield}")
        data = field.data.real
      else:
        print(f"Note: only plotting imaginary part of field {ifield}")
        data = field.data.imag
    else:
      if (imag==True):
        print("Warning: imag set to 'True' but fields are purely real")
      data = field.data

    # set title
    if nfields != 1:
      ax.set_title(f'Field {ifield}')

    if field.dim == 1:
        if nfields != 1:
          ax.set_title(f'Field {ifield}')
        ax.set_xlabel('x')
        ax.set_ylabel('field value')
        ax.plot(field.coords,data)

    if field.dim == 2:
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.axis('equal')
        pc = ax.pcolormesh(field.coords[:,:,0],field.coords[:,:,1], data,shading='auto',cmap='coolwarm')
        cb = fig.colorbar(pc,ax=ax)
        cb.set_label('field value')

    if field.dim == 3:
        from mpl_toolkits.mplot3d import Axes3D
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z',labelpad=0)
        scatter = ax.scatter(field.coords[:,:,:,0], field.coords[:,:,:,1],field.coords[:,:,:,2],c=data,cmap='coolwarm')
        
        # from https://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
        #divider = make_axes_locatable(ax)
        #cax = divider.append_axes("right", size="5%", pad=0.05)
        #cb = fig.colorbar(scatter,cax=cax)

        cb = fig.colorbar(scatter,ax=ax,pad=0.2,shrink=0.5)
        #cb.set_label('field value')
  plt.tight_layout()
  if filename:
    plt.savefig(filename)

  if show:
    plt.show()

  plt.close()


