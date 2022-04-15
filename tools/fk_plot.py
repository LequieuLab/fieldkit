#!/usr/bin/env python3
from fk_setpath import *
import fieldkit as fk 
import argparse


parser = argparse.ArgumentParser(description='Plot a field using matplotlib')
parser.add_argument('infile',help='Input filename of field')
parser.add_argument('-o','--outfile',default=None, help='Output filename of plot')
parser.add_argument('--show',default=True,action='store_true', help='show plot via interactive window?')
parser.add_argument('--dpi',default=100,type=float,help='resolution (dpi) of plot')
args = parser.parse_args()
#print(args)

# read file
fields = fk.read_from_file(args.infile)

# plot
fk.plot(fields, dpi=args.dpi, filename=args.outfile, show=args.show)

