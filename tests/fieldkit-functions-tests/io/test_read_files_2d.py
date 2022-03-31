import fieldkit as fk

def test_all():
    filename = "density.dat"

    fields = fk.read_from_file(filename)

    assert(len(fields) ==  2)
    assert(fields[0].npw == (32, 32))

    assert(fields[0].data[0,0] == 0.9867631917)
    assert(fields[0].data[0,10] == 0.0067130977)
    assert(fields[0].data[5,5] == 0.1642596991)
    assert(fields[0].data[20,1] == 0.0413537591)
    assert(fields[0].data[31,31] == 0.9816982147)

    #tests the read_from_file function in 2d
