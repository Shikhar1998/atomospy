"Plotting Lorentz transformation results"

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from sympy.utilities.decorator import doctest_depends_on
from sympy import S

@doctest_depends_on(modules=('matplotlib'))
def plot_lorentz(L, L_, style = 'ro'):
    """
    For plotting the four dimensional plots the
    fourth dimension in a three dimension

    Default style = 'ro'
    """

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if S.ComplexInfinity in L or S.ComplexInfinity in L_:
        raise("ComplexError: Cannot plot with complex quantities")
    ax.scatter([L[0], L_[0]], [L[1], L_[1]],
             [L[2], L_[2]], [L[2], L_[2]], cmap=plt.hot())
    plt.show()
