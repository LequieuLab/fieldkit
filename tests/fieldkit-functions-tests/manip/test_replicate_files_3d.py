import fieldkit as fk

def test_all():
    import numpy as np

    npw = (32, 32, 32)
    h = np.eye(3) * 2
    nreplicates = (2, 2, 2)

    fields = fk.initialize_phase("A15", npw, h)
    fields_new = fk.replicate_fields(fields, nreplicates)

    assert(fields_new[0].npw_Nd == (64, 64, 64))
    
    assert(fields_new[0].data[0,32,0] == fields[0].data[0,0,0])
    assert(fields_new[0].data[0,63,0] == fields[0].data[0,0,0])
    assert(fields_new[0].data[5,37,0] == fields[0].data[5,5,0])
    assert(fields_new[0].data[20,33,0] == fields[0].data[20,1,0])
    assert(fields_new[0].data[63,63,0] == fields[0].data[31,31,0])


    #test replicate_fields function for 3d fields
