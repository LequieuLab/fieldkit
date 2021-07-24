import fieldkit as fk

def test_all():
    filename = "density.dat"
    new_resolution = (64, 32)

    fields = fk.read_from_file(filename);
    new_fields = fk.change_resolution(fields, new_resolution)
 
    assert(new_fields[0].npw_Nd == (64, 32))

    #test change_resolution function for 2d fields
