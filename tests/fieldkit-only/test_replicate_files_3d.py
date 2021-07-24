import fieldkit as fk

def test_all():
    import numpy as np

    npw = (32, 32, 32)
    h = np.eye(3) * 2
    nreplicates = (2, 2, 2)

    fields = fk.initialize_phase("A15", npw, h)
    fields_new = fk.replicate_fields(fields, nreplicates)
 
    assert(fields_new[0].npw_Nd == (64, 64, 64))

    #test replicate_fields function for 3d fields
