
from .field import *
import matplotlib.pyplot as plt

def plot(fields,dpi=300):
  try:
    nfields = len(fields)
  except TypeError:
    fields = [fields]
    nfields = len(fields)
  
  fig = plt.figure(figsize=(3.33,3.33),dpi=dpi)

  # TODO: need to handle multiple fields by creating multiple axes
  ax = fig.subplots(1,1)
 
  for i,field in enumerate(fields):
    if field.dim == 1:
        ax.plot(field.coords,field.data,label=f'Field {i}' )
        ax.legend()
    if field.dim == 2:
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.axis('equal')
        # note x,y are reversed
        pc = ax.pcolormesh(field.coords[:,:,1],field.coords[:,:,0], field.data,shading='auto',cmap='coolwarm')
        fig.colorbar(pc,ax=ax)
  plt.show()


