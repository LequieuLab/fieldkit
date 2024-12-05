import fieldkit as fk

def test_all():
    dat_filename = "density_1D.dat"
    h5_filename = "density_1D.h5"

    dat_fields = fk.read_from_file(dat_filename)
    h5_fields = fk.read_from_HDF5(h5_filename)

    assert(len(dat_fields) ==  2)
    assert(dat_fields[0].npw == (64,))

    assert(dat_fields[0].data[0] == 0.0036036893)
    assert(dat_fields[0].data[5] == 0.0096679716)
    assert(dat_fields[0].data[20] == 0.7067403875)
    assert(dat_fields[0].data[32] == 0.9831727311)
    assert(dat_fields[0].data[63] == 0.0034669962)

    assert(len(h5_fields) ==  2)
    assert(h5_fields[0].npw == (64,))

    assert(h5_fields[0].data[0] == 0.0036036893)
    assert(h5_fields[0].data[5] == 0.0096679716)
    assert(h5_fields[0].data[20] == 0.7067403875)
    assert(h5_fields[0].data[32] == 0.9831727311)
    assert(h5_fields[0].data[63] == 0.0034669962)

    #tests the read_from_file function in 1d
