import math
import numpy as np
from matplotlib import pyplot as plt
from iptrack import *

import tkinter as tk
from tkinter import filedialog

## CONSTANTS
g = 9.81  # m/s
c = 2/5  # for solid sphere, I_0 / mr^2



## FUNCTIONS
def first_derivative(f, x, h=0.0001):
    """
    Computes first_derivative of function f at point x
    https://en.wikipedia.org/wiki/Numerical_differentiation#Finite_difference_formulas
    """
    return (f(x + h) - f(x)) / h


def second_derivative(f, x, h=0.001):
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


def normal_g_force(v_tangential, krumning):
    return 1 + v_tangential**2 * krumning / g


def bane(x, coeffs):
    # len - 1 because the last term is a constant, not x**1, and i starts at 0
    return sum(c * x**(len(coeffs) - 1 - i) for i, c in enumerate(coeffs))


def write_csv(*args, filename, header, ext="csv"):
    """
    *args: lists to be written
    filename: filename without extension
    header: list, file header with desc. of which values are where. should 
            correspond with *args
    ext: extension without dot
    """
    print("Writing to file...")
    f = open(filename + "." + ext, "w")
    print("#", ",".join(header), file=f)
    for vals in zip(*args):
        print(",".join("%.4f" % i for i in vals), file=f)
    f.close()
    print("Wrote to", filename + "." + ext)


def main():
    print("Enter file with data from Tracker")
    root = tk.Tk()
    root.withdraw()
    filename = filedialog.askopenfilename()
    root.destroy()  # removes window instance

    coeffs = iptrack(filename)
    func = lambda x: bane(x, coeffs)

    step_h = 0.001  # s, eulers method step
    UPPER_T = 1.1  # s, how long to simulate

    time_ = [i for i in range(0, int(UPPER_T / step_h))]

    angle = [0.0] * len(time_)
    pos_x = [0.0] * len(time_)
    pos_y = [0.0] * len(time_)
    speed_tangential = [0.0] * len(time_)
    accel_tangential = [0.0] * len(time_)
    normal = [0.0] * len(time_)
    curvature = [0.0] * len(time_)
    dfdx = [0.0] * len(time_)  # first derivative
    d2fdx2 = [0.0] * len(time_)  # second derivative

    angle[0] = alpha(func, pos_x[0])
    pos_y[0] = func(pos_x[0])

    # calc points for curve
    t = [i * step_h for i in range(0, int(UPPER_T / step_h))]
    b = list()
    for x in t:
        b.append(func(x))

    # simulation loop
    for n in time_[:-1]:  # n from 0 to len-1
        angle[n+1] = alpha(func, pos_x[n])
        accel_tangential[n+1] = accel_at_x(angle[n])
        speed_tangential[n+1] = speed_tangential[n] + step_h * accel_tangential[n]
        pos_x[n+1] = pos_x[n] + step_h * x_comp(speed_tangential[n], angle[n])
        pos_y[n+1] = pos_y[n] - step_h * y_comp(speed_tangential[n], angle[n])

        curvature[n+1] = krumning(func, pos_x[n])
        normal[n+1] = normal_g_force(speed_tangential[n], krumning(func, pos_x[n]))

        dfdx[n] = first_derivative(func, pos_x[n])
        d2fdx2[n] = second_derivative(func, pos_x[n])



    ans = input("Do you want to write simulation results to file? (y/N) ").lower()
    if ans.startswith("y"):
        # write simulation to file
        write_csv(time_, angle, pos_x, pos_y, speed_tangential, accel_tangential,\
                normal, curvature, dfdx, d2fdx2, filename="krumbane1_sim",
                header=["time","angle","x","y","speed_tangential",\
                "accel_tangential","normal_force","curvature","dfdx","d2fdx2"])

    # plotting
    print("Plotting...")
    fig = plt.figure("krumbane1")
    plt.title("Simulation results")
    plt.xlabel("Time step")
    plt.grid()

    plt.plot(time_, accel_tangential, "g")
    plt.plot(time_, curvature, "y")
    plt.plot(time_, speed_tangential, "b")
    plt.plot(time_, pos_x, "c")
    plt.plot(time_, pos_y, "m")
    plt.plot(time_, normal, "k--")
    # plt.plot(time_, dfdx, "m")
    # plt.plot(time_, d2fdx2, "c")

    # plt.plot(pos_x, pos_y, 'r--')
    # plt.plot(t, b)  # bane plot

    legends = (
                "Accel", 
                "Curvature",
                "Speed", 
                "X", 
                "Y", 
                "Normal G-force",
                # "df/dx",
                # "d^2f/dx^2",
                ""
            )
    plt.legend(legends, loc="best")

    plt.show()



if __name__ == '__main__':
    main()