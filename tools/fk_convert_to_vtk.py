#!/usr/bin/env python3
from fk_setpath import *
import fieldkit as fk 
import argparse

import os
import re
def atoi(text):
    return int(text) if text.isdigit() else text
def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]



parser = argparse.ArgumentParser(description='Convert field .dat files to VTK')
parser.add_argument('infiles',nargs='*', help='Input filenames of fields')
args = parser.parse_args()
#print(args)

# parser
args.infiles.sort(key=natural_keys)
print(f"Files to convert: {args.infiles}")

# figure out maximum frame number in inputfiles (used for zero padding)
# OPTION: could just set this to be some big value (i.e 20)
maxdigit=0
for infile in args.infiles:
    head,tail = os.path.split(infile)
    ndigit=sum( c.isdigit() for c in tail)  
    if ndigit > maxdigit: maxdigit = ndigit

# for each input file, generate a VTK file
for infile in args.infiles:

    # Generate the output file name automatically (padded with zeros)
    head,tail = os.path.split(infile)
    outfile,ext = os.path.splitext(tail)
    outfile_notdigits=''.join([c for c in outfile if not c.isdigit()])
    outfile_digits=''.join([c for c in outfile if c.isdigit()])
    zeropadded=str("{val:0>{width}}".format(val=outfile_digits,width=maxdigit))
    outfile = outfile_notdigits + zeropadded  + ".vtk"
   
    print(f"Converting {infile} -> {outfile}")
 
    # read file
    fields = fk.read_from_file(infile)

    # write to VTK
    fk.write_to_VTK(outfile,fields)

