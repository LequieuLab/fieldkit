import fieldkit as fk

def test_all():
    filename = "density.dat"

    fields = fk.read_from_file(filename)
    new_fields = fk.roll(fields, [0.5,0.5])

    assert(new_fields[0].data[0,0] == fields[0].data[16,16])
    assert(new_fields[0].data[0,10] == fields[0].data[16,26])
    assert(new_fields[0].data[5,5] == fields[0].data[21,21])
    assert(new_fields[0].data[20,1] == fields[0].data[4, 17])
    assert(new_fields[0].data[15, 15] == fields[0].data[31,31])
