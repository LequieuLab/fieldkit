import fieldkit as fk

def test_all():
    dat_filename = "density.dat"
    h5_filename = "density.h5"

    dat_fields = fk.read_from_file(dat_filename)
    h5_fields = fk.read_from_HDF5(h5_filename)

    assert(len(dat_fields) ==  2)
    assert(dat_fields[0].npw == (32, 32))

    assert(dat_fields[0].data[0,0] == 0.9867631917)
    assert(dat_fields[0].data[0,10] == 0.0067130977)
    assert(dat_fields[0].data[5,5] == 0.1642596991)
    assert(dat_fields[0].data[20,1] == 0.0413537591)
    assert(dat_fields[0].data[31,31] == 0.9816982147)

    assert(len(h5_fields) ==  2)
    assert(h5_fields[0].npw == (32, 32))

    assert(h5_fields[0].data[0,0] == 0.9867631917)
    assert(h5_fields[0].data[0,10] == 0.0067130977)
    assert(h5_fields[0].data[5,5] == 0.1642596991)
    assert(h5_fields[0].data[20,1] == 0.0413537591)
    assert(h5_fields[0].data[31,31] == 0.9816982147)


    #tests the read_from_file function in 2d
