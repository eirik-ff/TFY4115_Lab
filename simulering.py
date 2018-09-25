import math
import numpy as np
from matplotlib import pyplot as plt

## CONSTANTS
g = 9.81  # m/s
c = 2/5  # for solid sphere, I_0 / mr^2



## FUNCTIONS
def first_derivative(f, x, h=0.00001):
    """
    Computes first_derivative of function f at point x
    https://en.wikipedia.org/wiki/Numerical_differentiation#Finite_difference_formulas
    """
    return (f(x + h) - f(x)) / h


def second_derivative(f, x, h=0.00001):
    """
    https://www.scss.tcd.ie/~dahyotr/CS7ET01/01112007.pdf
    """
    return ( f(x + h) - 2 * f(x) + f(x - h) ) / h**2


def alpha(f, x):
    return math.atan(-1*(first_derivative(f, x)))


def accel_at_x(angle):
    """
    function to be used in eulers method
    """
    denom = 1 + c
    numer = g * math.sin(angle)

    return numer / denom


def x_comp(val, rads):
    return val * math.cos(rads)


def y_comp(val, rads):
    return val * math.sin(rads)


def krumning(f, x):
    denom = 1 + first_derivative(f, x)**2
    return abs(second_derivative(f, x)) / (denom**(3/2))


def bane(x):
    return 3433.2281869052126*x**15 + -35425.80967983638*x**14 + 164086.51268783736*x**13 + -449939.19577019487*x**12 + 810832.0563168445*x**11 + -1007456.1154303564*x**10 + 880102.071263024*x**9 + -540502.2906810011*x**8 + 228926.39127614404*x**7 + -63888.210540241824*x**6 + 10617.08175653279*x**5 + -759.7060900777494*x**4 + -33.33911458989251*x**3 + 8.555075814607894*x**2 + -1.7645318855075836*x**1 + 0.6532652470202632*x**0


def main():
    step_h = 0.001  # s
    UPPER_T = 1.1  # s, how long to simulate

    time_ = [i for i in range(0, int(UPPER_T / step_h))]

    angle = [0.0] * len(time_)
    pos_x = [0.0] * len(time_)
    pos_y = [0.0] * len(time_)
    speed_tangential = [0.0] * len(time_)
    accel_tangential = [0.0] * len(time_)

    angle[0] = alpha(bane, pos_x[0])
    pos_y[0] = bane(pos_x[0])

    # calc points for curve
    t = [i/100.0 for i in range(0, 140, 1)]
    b = list()
    for x in t:
        b.append(bane(x))

    # simulation loop
    for n in time_[:-1]:  # n from 0 to len-1
        angle[n+1] = alpha(bane, pos_x[n])
        accel_tangential[n+1] = accel_at_x(angle[n])
        speed_tangential[n+1] = speed_tangential[n] + step_h * accel_tangential[n]
        pos_x[n+1] = pos_x[n] + step_h * x_comp(speed_tangential[n], angle[n])
        pos_y[n+1] = pos_y[n] - step_h * y_comp(speed_tangential[n], angle[n])

    # write simulation to file
    f = open("simulation.txt", "w")
    print("Time\tAngle\tAccel\tSpeed\tX    \tY    ", file=f)
    for vals in zip(time_, angle, accel_tangential, speed_tangential, pos_x, pos_y):
        print(str(vals[0]) + "\t" + "\t".join("%.3f" % i for i in vals[1:]), file=f)
    f.close()
    print("Wrote to simulation.txt")

    # plotting
    print("Plotting...")
    fig = plt.figure()
    plt.title("Simulation results")
    plt.xlabel("Time step")
    plt.grid()

    plt.plot(time_, accel_tangential, "g")
    plt.plot(time_, speed_tangential, "b")
    plt.plot(time_, pos_x, "c")
    plt.plot(time_, pos_y, "m")
    # plt.plot(pos_x, pos_y, 'r--')
    # plt.plot(t, b)  # bane plot

    legends = ("Accel", "Speed", "X", "Y")
    plt.legend(legends, loc="best")

    plt.show()



if __name__ == '__main__':
    main()