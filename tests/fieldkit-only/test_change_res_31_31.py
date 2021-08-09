import fieldkit as fk 
from pytest import approx

def test_all():
    filename = "density.dat"
    new_resolution = (31, 31)

    fields = fk.read_from_file(filename)
    new_fields = fk.change_resolution(fields, new_resolution)
    fk.write_to_file("change_res_31_31.dat", new_fields)

    assert(new_fields[0].npw_Nd == (31,31))

    assert(approx(new_fields[0].data[0,0], 1e-06) == fields[0].data[0,0])
    assert(approx(new_fields[0].data[0,20], 1e-06) == fields[0].data[0,10])
    assert(approx(new_fields[0].data[10,10], 1e-06) == fields[0].data[5,5]) 
    assert(approx(new_fields[0].data[40,2], 1e-06) == fields[0].data[20,1])
    assert(approx(new_fields[0].data[62,62], 1e-06) == fields[0].data[31,31])

