#!/usr/bin/env python3
from fk_setpath import *
import fieldkit as fk 
import argparse


parser = argparse.ArgumentParser(description='Modify resolution of field so that voxels are cubic')
parser.add_argument('-i','--infile',required=True,help='Input filename of field')
parser.add_argument('-o','--outfile',required=True,help='Output filename of field')
parser.add_argument('-tol','--tol',required=True,type=float,help='error tolerance for how cubic voxels will be')
args = parser.parse_args()
#print(args)

# read file
fields = fk.read_from_file(args.infile)

# change resolution
new_fields = fk.cubic_voxels(fields, args.tol)

# write to file
fk.write_to_file(args.outfile,new_fields)

