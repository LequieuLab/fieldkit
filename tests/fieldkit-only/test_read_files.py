import fieldkit as fk

def test_all():
    filename = "density.dat"

    fields = fk.read_from_file(filename)

    assert(len(fields) ==  2)
    assert(fields[0].npw_Nd == (32, 32))

    #tests the read_from_file function
