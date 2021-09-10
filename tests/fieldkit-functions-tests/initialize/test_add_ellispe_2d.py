import fieldkit as fk
from pytest import approx

def test_units_scaled():
    npw = (64,64)
    
    field1 = fk.Field(npw_Nd = npw)
    field2 = fk.Field(npw_Nd = npw)

    fk.add_ellipse(field1, center=(0.0,0.0), axis_lengths=(0.5, 0.5), height=1)
    fk.add_ellipse(field2, center=(0.0,0.0), axis_lengths=(0.3, 0.2), height=1)

    assert(approx(field1.data[0,0], 1e-02) == 1.00)
    assert(approx(field1.data[0,63], 1e-02) == 1.00)
    assert(approx(field1.data[63,0], 1e-02) == 1.00)
    assert(approx(field1.data[63,63], 1e-02) == 1.00)
    assert(approx(field1.data[32,32], 1e01) == 0.00)

    assert(approx(field2.data[0,0], 1e-02) == 1.00)
    assert(approx(field2.data[0,63], 1e-02) == 1.00)
    assert(approx(field2.data[63,0], 1e-02) == 1.00)
    assert(approx(field2.data[63,63], 1e-02) == 1.00)
    assert(approx(field2.data[19,12], 1e-01) == 0.022)

#TODO test unscaled units
