"""
Implementation of Relativistic Doppler effect
"""
from atomospy import c
from sympy import S, sqrt

def beta(v):
    return v/S(c)

def doppler(v_source):
    """
    If the sign of v_source is positive then the source
    is moving towards the observer and viceversa.
    """
    if v_source>c:
        raise(ValueError("v_source cannot be greater than c"))
    return v_source*sqrt((1+beta(v_source))/(1-beta(v_source)))
