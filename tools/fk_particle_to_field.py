#!/usr/bin/env python3
from fk_setpath import *
import fieldkit as fk 
import argparse


parser = argparse.ArgumentParser(description='Change resolution of a field')
parser.add_argument('-trj','--trjfile',required=True,help='Input trajectory file')
parser.add_argument('-psf','--psffile',required=True,help='Input psf tolology file')
parser.add_argument('-res','--resolution',nargs=3,required=True,type=int,help='New resolution of fields ')
parser.add_argument('-f','--frame',default=-1,type=int,help='frame index to use')
parser.add_argument('-P','--assignment_order',default=4,type=int,help='frame index to use')
parser.add_argument('-o','--outfile',default='particle_density.dat',type=str,help='frame index to use')
args = parser.parse_args()
#print(args)

trjfile = args.trjfile
psffile = args.psffile
frame_index = args.frame 
P = args.assignment_order
npw = args.resolution
vtkfile = ''.join(args.outfile.split('.')[0:-1]) + '.vtk'

fields = fk.particle_to_field_hockney_eastwood(trjfile,psffile, frame_index, npw, P)

# write to file
fk.write_to_file(args.outfile,fields)

# write VTK
fk.write_to_VTK(vtkfile,fields)

