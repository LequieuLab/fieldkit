import fieldkit as fk

def test_all():
    import numpy as np

    npw = (32, 32, 32)
    h = np.eye(3) * 2
    new_resolution = (64, 64, 64)

    fields = fk.initialize_phase("A15", npw, h)
    field = fields[0]
    new_fields = fk.change_resolution(fields, new_resolution)
 
    assert(new_fields[0].npw_Nd == (64, 64, 64))

    #test change_resolution for 3d phases
