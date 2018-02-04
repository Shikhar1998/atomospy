from atomospy import c, gamma, lorentz_transform
from sympy import S, zoo, nan
from sympy.physics.vector import *

def test_gamma():
    assert gamma(0) == S(1)
    assert gamma(c) == S.ComplexInfinity

def test_lorentz_transform():
    N = ReferenceFrame("N")
    L = [1, 2, 3, 4]
    # Testing one dimensional lorentz transformation
    v = c*N.x + 0*N.y + 0*N.z
    assert lorentz_transform(N, L, v) == [zoo, 2.0000, 3.0000, zoo]
    # Testing two dimensional lorentz transformation
    v = c*N.x + 2*N.y + 0*N.z
    assert lorentz_transform(N, L, v) == [zoo, -6.0000, 3.0000, zoo]
    v = c*N.x + c*N.y + 0*N.z
    assert lorentz_transform(N, L, v) == [zoo, zoo, 3.0000, nan]
    # Testing three dimensional lorentz transformation
    v = 1*N.x + 2*N.y + 4*N.z
    assert lorentz_transform(N, L, v) == [-3.0000, -6.0000, -13.000, 12.000]
    v = c*N.x + c*N.y + c*N.z
    assert lorentz_transform(N, L, v) == [zoo, zoo, zoo, nan]
