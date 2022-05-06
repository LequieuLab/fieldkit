#!/usr/bin/env python3
from fk_setpath import *
import fieldkit as fk 
import argparse

parser = argparse.ArgumentParser(description='Change resolution of a field')
parser.add_argument('-gsd','--gsdfile',required=True,help='Input gsd trajectory file')
parser.add_argument('-res','--resolution',nargs=3,required=True,type=int,help='New resolution of fields ')
parser.add_argument('-f','--frames',default=[-1],nargs='+',type=int,help='frame indicies to use')
parser.add_argument('-P','--assignment_order',default=4,type=int,help='order of Hockney-Eastwood assignment function')
parser.add_argument('-o','--outfile',default='particle_density.dat',type=str,help='name of output filename')
parser.add_argument('--no_normalize',action='store_true',help='do not normalize the fields by rho0')
args = parser.parse_args()
#print(args)

P = args.assignment_order
npw = args.resolution
#vtkfile = ''.join(args.outfile.split('.')[0:-1]) + '.vtk'

if args.no_normalize:
  normalize = False
else:
  normalize = True

fields = fk.particle_to_field_gsd(args.gsdfile, args.frames, npw, P, normalize=normalize, use_jit=True)

# write to file
fk.write_to_file(args.outfile,fields)

# write VTK
#fk.write_to_VTK(vtkfile,fields)
