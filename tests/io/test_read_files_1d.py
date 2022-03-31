import fieldkit as fk

def test_all():
    filename = "density_1D.dat"

    fields = fk.read_from_file(filename)

    assert(len(fields) ==  2)
    assert(fields[0].npw == (64,))

    assert(fields[0].data[0] == 0.0036036893)
    assert(fields[0].data[5] == 0.0096679716)
    assert(fields[0].data[20] == 0.7067403875)
    assert(fields[0].data[32] == 0.9831727311)
    assert(fields[0].data[63] == 0.0034669962)

    #tests the read_from_file function in 1d
