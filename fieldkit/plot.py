
from .field import *
import matplotlib.pyplot as plt

def plot(fields,dpi=100):
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
    nrows = 1; ncols = 2
  elif nfields == 4:
    nrows = 2; ncols = 2
  else:
    raise RuntimeError("nfields > 4 not currently supported")
 
  fig = plt.figure(figsize=(ncols*3.33,nrows*3.33),dpi=dpi)

  # TODO: need to handle multiple fields by creating multiple axes
  axes = fig.subplots(nrows,ncols,squeeze=False)
 
  for i,ax in np.ndenumerate(axes):
    ifield = i[0]*ncols + i[1]*nrows
    field = fields[ifield]

    if field.dim == 1:
        ax.plot(field.coords,field.data,label=f'Field {i}' )
        ax.legend()

    if field.dim == 2:
        if nfields != 1:
          ax.set_title(f'Field {ifield}')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.axis('equal')
        pc = ax.pcolormesh(field.coords[:,:,0],field.coords[:,:,1], field.data.real,shading='auto',cmap='coolwarm')
        fig.colorbar(pc,ax=ax)

    if field.dim == 3:
        pass

  plt.tight_layout()
  plt.show()


