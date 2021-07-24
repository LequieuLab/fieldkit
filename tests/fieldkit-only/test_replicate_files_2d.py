import fieldkit as fk

def test_all():
    filename="density.dat"
    nreplicates = (4,4)

    fields = fk.read_from_file(filename)
    fields_new = fk.replicate_fields(fields, nreplicates)
 
    assert(fields_new[0].npw_Nd == (128,128))

    #test replicate_fields function in 2d fields
