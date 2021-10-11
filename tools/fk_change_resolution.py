#!/usr/bin/env python3
from fk_setpath import *
import fieldkit as fk 
import argparse


parser = argparse.ArgumentParser(description='Change resolution of a field')
parser.add_argument('-i','--infile',required=True,help='Input filename of field')
parser.add_argument('-o','--outfile',required=True,help='Output filename of field')
parser.add_argument('-res','--resolution',nargs='+',required=True,type=int,help='New resolution of fields ')
args = parser.parse_args()
#print(args)

# read file
fields = fk.read_from_file(args.infile)

# confirm that dimention of new resolution provided matches that in infile
dim = fields[0].dim
dim_new = len(args.resolution)
assert(dim == dim_new), f"Input file is {dim}d, but new resolution is {dim_new}d"

# change resolution
new_fields = fk.change_resolution(fields, args.resolution)

# write to file
fk.write_to_file(args.outfile,new_fields)

