#!/usr/bin/env python3
from fk_setpath import *
import fieldkit as fk 
import argparse


parser = argparse.ArgumentParser(description='Change resolution of a field')
parser.add_argument('-i','--infile',required=True,help='Input filename of fields')
parser.add_argument('-o','--outfile',required=True,help='Output filename of fields')
parser.add_argument('-d','--dim_new',required=True,type=int,help='New dimension of fields ')
parser.add_argument('-npw','--npw_new',nargs='+',required=True,type=int,help='New resolution of fields ')
parser.add_argument('-boxl','--boxl_new',nargs='+',required=True,type=float,help='New box length of expanded fields ')
args = parser.parse_args()
#print(args)

# read file
fields = fk.read_from_file(args.infile)

# expand dimension
new_fields = fk.expand_dimension(fields, args.dim_new, args.npw_new, args.boxl_new)

# write to file
fk.write_to_file(args.outfile,new_fields)

