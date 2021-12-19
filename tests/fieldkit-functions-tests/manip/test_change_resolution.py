#!/usr/bin/env python3 

import fieldkit as fk 
from pytest import approx

def test_32x32_to_64x64():
    filename = "density.dat"
    new_resolution = (64, 64)

    fields = fk.read_from_file(filename)
    new_fields = fk.change_resolution(fields, new_resolution)

    assert(new_fields[0].npw_Nd == (64,64))

    assert(approx(new_fields[0].data[0,0], 1e-06) == fields[0].data[0,0])
    assert(approx(new_fields[0].data[0,20], 1e-06) == fields[0].data[0,10])
    assert(approx(new_fields[0].data[10,10], 1e-06) == fields[0].data[5,5]) 
    assert(approx(new_fields[0].data[40,2], 1e-06) == fields[0].data[20,1])
    assert(approx(new_fields[0].data[62,62], 1e-06) == fields[0].data[31,31])

def test_32x32_to_64x32():
    filename = "density.dat"
    new_resolution = (64, 32)

    fields = fk.read_from_file(filename)
    new_fields = fk.change_resolution(fields, new_resolution)

    assert(new_fields[0].npw_Nd == (64,32))

    assert(approx(new_fields[0].data[0,0], 1e-06) == fields[0].data[0,0])
    assert(approx(new_fields[0].data[0,10], 1e-06) == fields[0].data[0,10])
    assert(approx(new_fields[0].data[10,5], 1e-06) == fields[0].data[5,5]) 
    assert(approx(new_fields[0].data[40,1], 1e-06) == fields[0].data[20,1])
    assert(approx(new_fields[0].data[62,31], 1e-06) == fields[0].data[31,31])

