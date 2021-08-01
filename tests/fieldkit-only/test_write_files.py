import fieldkit as fk

def test_all():
    import os
    
    filename = "density.dat"
    
    fields = fk.read_from_file(filename)
    fk.write_to_file("tmp.dat", fields)
    fk.write_to_VTK("tmp.vtk", fields)

    #test if the files are created
    assert(os.path.exists("tmp.dat"))
    assert(os.path.exists("tmp.vtk"))

    #test if the files are not empty
    assert(os.path.getsize("tmp.dat") != 0)
    assert(os.path.getsize("tmp.vtk") != 0)

    #test the write_to_file and write_to_VTK functions
