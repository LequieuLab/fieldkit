import fieldkit as fk

def test_all():
    filename = "density.dat"

    fields = fk.read_from_file(filename)
    new_fields = fk.roll(fields, [0.25,0.5])

    assert(new_fields[0].data[8,16] == fields[0].data[0,0])
    assert(new_fields[0].data[8,26] == fields[0].data[0,10])
    assert(new_fields[0].data[13,21] == fields[0].data[5,5])
    assert(new_fields[0].data[28,17] == fields[0].data[20, 1])
    assert(new_fields[0].data[7, 15] == fields[0].data[31,31])
