"""
Implementation of lorentz transformation
"""
from sympy import S, Pow
from .constants import c

def _gamma(v):
    return S(1)/Pow(1 - v**2/c**2, 1/2)

def lorentz_space(L, v):
    """
    Reference: https://en.wikipedia.org/wiki/Lorentz_transformation

    Paramters:

    L: space time array (x, y, z, t) in frame F
    v: velocity of frame F* with respect to F
    """
    g = gamma(v)
    f_ = L
    

    
