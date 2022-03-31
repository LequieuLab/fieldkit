import fieldkit as fk
from pytest import approx

def test_all():
    npw = [64]

    field1 = fk.Field(npw = npw)
    field2 = fk.Field(npw = npw)
    field3 = fk.Field(npw = npw)

    fk.add_ellipse(field1, center=[0.0], axis_lengths=[0.5])
    fk.add_ellipse(field2, center=[0.5], axis_lengths=[0.5])
    fk.add_ellipse(field3, center=[0.5], axis_lengths=[0.25])

    assert(approx(field1.data[0],1e-02) == 1.00)
    assert(approx(field1.data[63], 1e-02) == 1.00)
    assert(approx(field1.data[32], 1e-03) == 0.867)

    assert(approx(field2.data[0],1e-02) == 0.867)
    assert(approx(field2.data[63], 1e-02) == 0.867)
    assert(approx(field2.data[32], 1e-02) == 1.00)

    assert(approx(field3.data[0],1e-02) == 0.00)
    assert(approx(field3.data[63], 1e-02) == 0.00)
    assert(approx(field3.data[32], 1e-02) == 1.00)


