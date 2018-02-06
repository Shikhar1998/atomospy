from atomospy import c, doppler
from sympy import S

def test_doppler():
    assert doppler(0) == 0
    assert doppler(c) == S.ComplexInfinity
