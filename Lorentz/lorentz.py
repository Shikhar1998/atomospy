"""
Implementation of lorentz transformation
"""
from sympy import S, Pow, Mul, Add, evalf
from constants import c
from sympy.physics.vector import *

def gamma(v):
    """
    This is the rescaling factor for a particle
    in motion along the x-direction
    """
    return S(1)/Pow(1 - v**2/c**2, S(1)/2)

def lorentz_transform(N, L, v, evaluate=True, precision=5):
    """
    Reference: https://en.wikipedia.org/wiki/Lorentz_transformation

    Paramters:

    N: Input Reference Frame
    L: space time array (x, y, z, t) in frame F
    v: vector, velocity of frame F* with respect to F 
    L_: space time array (x*, y*, z*, t*) in frame F*

    from atomospy.py import lorentz_transform
    from sympy.physics.vector import *
    >>> N = ReferenceFrame("N")
    >>> L = [1, 2, 3, 4]
    >>> v = c*N.x + N.y + N.z
    >>> lorentz_transform(N, L, v)
    [zoo, -2.0, -1.0, zoo]
    """
    
    g_x = gamma(N.x.dot(v))
    g_y = gamma(N.y.dot(v))
    g_z = gamma(N.z.dot(v))

    g = [g_x, g_y, g_z]
    v_val = [N.x.dot(v), N.y.dot(v), N.z.dot(v)]

    L_ = [0, 0, 0, 0]
    for i in range(3):
            L_[i] = Mul(g[i], L[i] - v_val[i] * L[3]) 
            L_[3] = Add(L_[3], Mul(g[i], L[3] - v_val[i] * L[i]/ c**2))
    if evaluate is True:
        for j in range(4):
            L_[j] = L_[j].evalf(precision)
    return L_
