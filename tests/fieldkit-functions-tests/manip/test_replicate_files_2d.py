import fieldkit as fk

def test_all():
    filename="density.dat"
    nreplicates = (4,4)

    fields = fk.read_from_file(filename)
    fields_new = fk.replicate_fields(fields, nreplicates)
 
    assert(fields_new[0].npw == (128,128))


    assert(fields_new[0].data[0,32] == fields[0].data[0,0])
    assert(fields_new[0].data[0,42] == fields[0].data[0,10])
    assert(fields_new[0].data[5,37] == fields[0].data[5,5])
    assert(fields_new[0].data[20,65] == fields[0].data[20,1])
    assert(fields_new[0].data[127,127]==fields[0].data[31,31])


    #test replicate_fields function in 2d fields
