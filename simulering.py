import math
import numpy as np

## CONSTANTS
g = 9.81  # m/s
I_0 = 1
m = 0.0027  # kg
r = 0.001  # m


## FUNCTIONS
def first_derivative(f, x, h=0.000001):
    """
    Computes first_derivative of function f at point x
    https://en.wikipedia.org/wiki/Numerical_differentiation#Finite_difference_formulas
    """
    return ( f(x + h) - f(x - h) ) / h**2


def second_derivative(f, x, h=0.0000001):
    """
    https://www.scss.tcd.ie/~dahyotr/CS7ET01/01112007.pdf
    """
    return ( f(x + h) - 2 * f(x) + f(x - h) ) / h**2


def alpha(x, f):
    return math.atan(abs(first_derivative(f, x)))


def accel_at_x(angle):
    """
    function to be used in eulers method
    """
    denom = 1 + I_0 / (m * r**2)
    numer = g * math.sin(angle)

    return numer / denom


def x_comp(val, rads):
    return val * math.cos(rads)


def krumning(f, x):
    numer = 1 + first_derivative(f, x)**2
    return abs(second_derivative(f, x)) / (numer**(3/2))


def bane(x):
    return x**2


def main()
    step_h = 0.01  # s
    upper_t = 5  # s

    time_ = [i for i in range(0, upper_t / step_h)]

    pos_x = [0] * len(time_)
    speed_x = [0] * len(time_)
    accel_tangent = [0] * len(time_)

    for n in nme_[:-1]:  # n from 0 to len-1
        angle_rad = alpha(bane, pos_x[n])
        accel_tangent[n+1] = 


if __name__ == '__main__':
    main()