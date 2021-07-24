import fieldkit as fk

def test_all():
    import os
    
    filename = "density.dat"
    
    fields = fk.read_from_file(filename)
    fk.write_to_file("tmp.dat", fields)
    fk.write_to_VTK("tmp.vtk", fields)

    assert(os.path.exists("tmp.dat"))
    assert(os.path.exists("tmp.vtk"))

    #test the write_to_file and write_to_VTK functions
