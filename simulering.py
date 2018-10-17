import os
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
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
    """
    in radians
    """
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
    filepath = filename + "." + ext

    if filepath in os.listdir(os.getcwd()):
        ans = input("Overwrite old data? (Y/n) ").lower()
        if ans == "n":
            count = 0
            for f in os.listdir(os.getcwd()):
                if filename in f:
                    count += 1

            filepath = filename + "_" + str(count) + "." + ext

    print("Writing to file", filepath)
    f = open(filepath, "w")
    print("#", ",".join(header), file=f)
    for vals in zip(*args):
        print(",".join("%.4f" % i for i in vals), file=f)
    f.close()
    print("Wrote to", filepath)


def plot_g_kraft_og_bane_vs_x(pos_x, pos_y, normal, fig):
    ax1 = fig.add_subplot(111)
    plt.xlabel("x-posisjon $x$ [m]")
    plt.grid()

    color = "r"
    ax1.set_ylabel("G-krefter $n_{gf}$ [1]", color=color)
    ax1.plot(pos_x, normal, color + "-")

    ax2 = ax1.twinx()
    color = "g"
    ax2.set_ylabel("y-posisjon $y$ [m]", color=color)
    ax2.plot(pos_x, pos_y, color + "-.")

    fig.tight_layout()


def main():
    print("Enter file with data from Tracker")
    # root = tk.Tk()
    # root.withdraw()
    # filepath = filedialog.askopenfilename()
    # if filepath == "":
    #     print("No file chosen, exiting...")
    #     return
    # root.destroy()  # removes window instance
    filepath = r"C:/Users/eirik/googledrive/1_Skole/1_Universitetet/_3_Semester/TFY4115_Fysikk/Lab/Programfiler/test2.txt"
    filename = os.path.basename(filepath).split(".")[0]

    data = np.loadtxt(filepath, skiprows=2)  # t x y v
    coeffs = iptrack_data(data)
    track = lambda x: bane(x, coeffs)  # banen

    # finne bunnpunktet
    y_data = data[:,2]
    y_min = np.argmin(y_data)

    # eksperimentell data 
    curvature_exp = krumning(track, data[:,1])
    speed_exp = data[:,3]
    normal_g_force_exp = 1 + speed_exp**2 * curvature_exp / g

    step_h = 0.001  # s, eulers method step
    UPPER_T = 1.1  # s, how long to simulate

    # find upper time limit
    with open(filepath, "r") as f:
        lines = f.readlines()
        last = lines[-1]
        t, x, y, v = last.split()
        UPPER_T = float(t)

    time_ = [i for i in range(0, int(UPPER_T / step_h))]
    time_plot = [i*step_h for i in time_]  # x values in seconds

    angle = [0.0] * len(time_)
    pos_x = [0.0] * len(time_)
    pos_y = [0.0] * len(time_)
    speed_tangential = [0.0] * len(time_)
    accel_tangential = [0.0] * len(time_)
    normal = [0.0] * len(time_)
    curvature = [0.0] * len(time_)
    dfdx = [0.0] * len(time_)  # first derivative
    d2fdx2 = [0.0] * len(time_)  # second derivative

    pos_x[0] = data[:,1][0]
    angle[0] = alpha(track, pos_x[0])
    pos_y[0] = track(pos_x[0])
    speed_tangential[0] = data[:,3][0]

    # calc points for curve
    t = [i * step_h for i in range(0, int(UPPER_T / step_h))]
    b = list()
    for x in t:
        b.append(track(x))

    # simulation loop
    for n in time_[:-1]:  # n from 0 to len-1
        # simulering 
        angle[n+1] = alpha(track, pos_x[n])
        accel_tangential[n+1] = accel_at_x(angle[n])
        speed_tangential[n+1] = speed_tangential[n] + step_h * accel_tangential[n]
        pos_x[n+1] = pos_x[n] + step_h * x_comp(speed_tangential[n], angle[n])
        pos_y[n+1] = pos_y[n] - step_h * y_comp(speed_tangential[n], angle[n])

        curvature[n+1] = krumning(track, pos_x[n])
        normal[n+1] = normal_g_force(speed_tangential[n], krumning(track, pos_x[n]))

        dfdx[n] = first_derivative(track, pos_x[n])
        d2fdx2[n] = second_derivative(track, pos_x[n])



    # ans = input("Do you want to write simulation results to file? (y/N) ").lower()
    ans = "n"
    if ans == "y":
        # write simulation to file
        write_csv(time_, angle, pos_x, pos_y, speed_tangential, accel_tangential,\
                normal, curvature, dfdx, d2fdx2, filename=filename + "_sim",
                header=["time","angle","x","y","speed_tangential",\
                "accel_tangential","normal_force","curvature","dfdx","d2fdx2"])

    # plotting
    print("Plotting...")
    fig = plt.figure(filename, figsize=(6,4), dpi=180)
    # plot_g_kraft_og_bane_vs_x(pos_x, pos_y, normal, fig)
    # plot_g_kraft_og_bane_vs_x(pos_x, pos_y, normal, fig)

    # plt.plot(data[:,1], normal_g_force_exp)
    # plt.axvline(x=data[y_min][1], color="k", linewidth=0.7)
    plt.plot(data[:,0], data[:,3])
    plt.plot(time_plot, speed_tangential)

    print("G-kraft i bunnpunkt:", normal_g_force_exp[y_min], normal[10 * y_min])
    print(normal[10*y_min - 10:10*y_min + 11])

    # plt.plot(pos_x, pos_y, 'r--')
    # plt.plot(t, b)  # bane plot

    legends = (
                # "Accel", 
                # "Curvature",
                # "Speed", 
                # "X", 
                # "Y", 
                "Normal G-force",
                # "df/dx",
                # "d^2f/dx^2",
                ""
            )
    # plt.legend(legends, loc="best")

    # ans = input("Save svg figure? (y/N) ").lower()
    ans = "y"
    if ans == "y":
        plt.savefig(filename + ".svg", format="svg")
    plt.show()


if __name__ == '__main__':
    main()