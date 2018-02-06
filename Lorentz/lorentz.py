"""
Implementation of Lorentz Transformation
"""

from sympy import S, Pow, Mul, Add, evalf, sin, cos, Matrix, sqrt
#from atomospy import c
c = 299792458
from sympy.physics.vector import *

def gamma(v):
    """
    This is the rescaling factor for a particle
    in motion along the x, y or z direction.
    """

    return S(1)/Pow(1 - v**2/c**2, S(1)/2)

def lorentz_1d(L, v, evaluate=True, precesion=5):
    """
    Parameters:
    ==========

    L: 1-d space time array (x, t) in frame F
    v: vector, velocity of frame F* with respect to F
    L_: space time array (x*, t*) in frame F*
    evaluate: Calculate the equivalent sympy expression
    precesion: Precesion desired for the answer

    This function calculates the one dimensional Lorentz
    Transform for a single point.

    Examples
    ========

    from atomospy.py import lorentz_transform
    >>> L = [1, 3]
    >>> v = c/2
    >>> lorentz_1d(L, v, True, 15)
    Matrix([[1.0000], [3.0000], [-1.6678e-9]])
    """

    if len(L)<2:
         raise(TypeError("Incorrect input dimension size"))
    L.insert(2, 0)
    return lorentz_transform(L, v, 0)

def lorentz_transform(L, v, theta, evaluate=True, precision=5):
    """
    Reference: https://en.wikipedia.org/wiki/Lorentz_transformation

    Paramters:
    =========

    L: space time array (x, y, t) in frame F
    v: vector, velocity of frame F* with respect to F 
    L_: space time array (x*, y*, t*) in frame F*

    Examples
    ========

    from atomospy.py import lorentz_transform
    >>> L = [1, 2, 3]
    >>> v = c/2
    >>> theta = 0
    >>> lorentz_transform(L, v, theta, True, 15)
    Matrix([[3.4162e+8], [-2.9243e+8], [3.0000]])
    """

    if v>c:
        raise(ValueError("The magnitude can never be greater than speed of light"))
    if len(L)<3:
        raise(TypeError("Incorrect input dimension size"))

    g = gamma(v)

    array_1 = [Mul(g, cos(theta)**2) + sin(theta)**2, Mul(Mul(sin(theta), cos(theta)), (g - 1)), -Mul(Mul(v, g), cos(theta))]
    array_2 = [Mul(Mul(sin(theta), cos(theta)), (g - 1)), Mul(g, sin(theta)**2) + cos(theta)**2, -Mul(Mul(v, g), sin(theta))]
    array_3 = [-Mul(Mul(v, g), cos(theta))/ (c**2), -Mul(Mul(v, g), sin(theta))/ (c**2), g]
    transformation_matrix = Matrix([array_1, array_2, array_3])
    L = Matrix(L)
    L_ = transformation_matrix * L
    if evaluate is True:
        L_ = L_.evalf(precision)
    return L_
L = [1, 2]
v = c/2
print lorentz_1d(L, v, evaluate=True, precesion=25)
