from atomospy import c, gamma, lorentz_transform, lorentz_1d
from sympy import S, zoo, nan, Matrix

def test_gamma():
    assert gamma(0) == S(1)
    assert gamma(c) == S.ComplexInfinity

def test_lorentz_1d():
    L = [1, 2]
    v = c/2
    assert lorentz_1d(L, v, evaluate=True, precesion=25) == \
           Matrix([[1.0000], [2.0000], [-1.6678e-9]])

def test_lorentz_transform():
    L = [1, 2, 3]
    v = 0
    assert lorentz_transform(L, v, 0, evaluate=True, precision=25) == \
           Matrix([[1.00000000000000], [2.00000000000000], [3.00000000000000]])
    v = c/4
    assert lorentz_transform(L, v, 90, evaluate=True, precision=15) == \
           Matrix([[100746818.386124], [-201010089.577463], [2.99999999888263]])
    v = c/2
    assert lorentz_transform(L, v, 45, evaluate=True, precision=15) == \
           Matrix([[-236231354.403674], [-382641686.711420], [2.99999999628555]])
    v = 3*c/4
    assert lorentz_transform(L, v, 135, evaluate=True, precision=15) == \
           Matrix([[671894145.587836], [-59607595.5064823], [3.00000000204979]])
    v = c-0.0000001
    assert lorentz_transform(L, v, 180, evaluate=True, precision=25) == \
           Matrix([[20854337462797735.66744330], [27917497404806660.16951075], [116235962.3706579454679635]])
