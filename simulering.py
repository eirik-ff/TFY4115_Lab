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
m = 0.0303  # kg


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
    return np.arctan(-1*(first_derivative(f, x)))


def accel_at_x(angle):
    """
    function to be used in eulers method
    """
    denom = 1 + c
    numer = g * np.sin(angle)

    return numer / denom


def x_comp(val, rads):
    return val * np.cos(rads)


def y_comp(val, rads):
    return val * np.sin(rads)


def krumning(f, x):
    denom = 1 + first_derivative(f, x)**2
    return abs(second_derivative(f, x)) / (denom**(3/2))


def normal_g_force(v_tangential, krumning, angle):
    return np.cos(angle) + v_tangential**2 * krumning / g


def friction(alpha, m, a):
    """
    alpha: vinkel (rad)
    m: masse
    a: akselerasjon
    """
    return m * (g * np.sin(alpha) - a)


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


def main():
    print("Enter file with data from Tracker")
    # root = tk.Tk()
    # root.withdraw()
    # filepath = filedialog.askopenfilename()
    # if filepath == "":
    #     print("No file chosen, exiting...")
    #     return
    # root.destroy()  # removes window instance
    filepath = r"./data/fart1.txt"
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
    angles_exp = alpha(track, data[:,1])
    # normal_g_force_exp = 1 + speed_exp**2 * curvature_exp / g
    normal_g_force_exp = normal_g_force(speed_exp, curvature_exp, angles_exp)

    step_h = 0.001  # s, eulers method step
    LOWER_T = 0
    UPPER_T = 1.1  # s, how long to simulate

    # find upper time limit
    with open(filepath, "r") as f:
        lines = f.readlines()
        first = lines[2]
        last = lines[-1]
        t, x, y, v = last.split()
        UPPER_T = float(t)
        t, x, y, v = first.split()
        LOWER_T = float(t)

    time_ = [i for i in range(int(LOWER_T / step_h), int(UPPER_T / step_h))]
    time_plot = [i*step_h for i in time_]  # x values in seconds

    angle = [0.0] * len(time_)
    pos_x = [0.0] * len(time_)
    pos_y = [0.0] * len(time_)
    speed_tangential = [0.0] * len(time_)
    accel_tangential = [0.0] * len(time_)
    normal = [0.0] * len(time_)
    curvature = [0.0] * len(time_)
    friction_sim = [0.0] * len(time_)
    dfdx = [0.0] * len(time_)  # first derivative
    d2fdx2 = [0.0] * len(time_)  # second derivative

    pos_x[0] = data[:,1][0]
    angle[0] = alpha(track, pos_x[0])
    pos_y[0] = track(pos_x[0])
    speed_tangential[0] = data[:,3][0]

    # simulation loop
    # for n in time_[:-1]:  # n from 0 to len-1
    for n in range(len(time_) - 1):
        # simulering 
        angle[n+1] = alpha(track, pos_x[n])
        accel_tangential[n+1] = accel_at_x(angle[n])
        speed_tangential[n+1] = speed_tangential[n] + step_h * accel_tangential[n]
        pos_x[n+1] = pos_x[n] + step_h * x_comp(speed_tangential[n], angle[n])
        pos_y[n+1] = pos_y[n] - step_h * y_comp(speed_tangential[n], angle[n])

        curvature[n+1] = krumning(track, pos_x[n])
        normal[n+1] = normal_g_force(speed_tangential[n], krumning(track, pos_x[n]), angle[n])
        friction_sim[n+1] = m * (g * np.sin(angle[n+1]) + accel_tangential[n+1])

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
    plt.grid()


    # graf som sammenlikner normalkraft simulert vs faktisk
    ax1 = fig.add_subplot(111)
    plt.xlabel("$x$-posisjon [m]")

    color = "r"
    ax1.set_ylabel("G-krefter $n_{g}$ [1]")
    l1, = ax1.plot(pos_x, normal, color + "-")
    l3, = ax1.plot(data[:,1], normal_g_force_exp)  # eksperimentell data
    ax1.set_ylim(0, 3.1)

    ax2 = ax1.twinx()
    color = "g"
    ax2.set_ylabel("$y$-posisjon [m]", color=color)
    ax2.set_ylim(0, 1.0323)
    l2, = ax2.plot(pos_x, pos_y, color + "-.")
    l4  = plt.axvline(x=data[y_min][1], color="k", linewidth=0.7)
    legends = (
                "Simulering",
                "Bane",
                "Eksp. data",
                "Bunnpunkt"
            )
    plt.legend([l1, l2, l3, l4], legends, loc="upper left")
    fig.tight_layout()

    #################################################################

    # # graf som sammenlikner simuleringen og eksperimentell data
    # plt.xlabel("Tid [s]")
    # plt.ylabel("$y$-posisjon [m]")
    # plt.plot(data[:,0], data[:,2])
    # plt.plot(time_plot, pos_y, "r")
    # legends = (
    #             "Eksp. data",
    #             "Simulering"
    #         )
    # plt.legend(legends, loc="best")

    #################################################################

    # plt.xlabel("Tid [s]")
    # # graf som viser friksjonskraften og y-posisjonen
    # ax1 = fig.add_subplot(111)
    # color = "b"
    # ax1.set_ylabel("Friksjons- og normalkraft [N]")
    # l1, = ax1.plot(time_plot, friction_sim, color + "-")
    # color = "r"
    # l5, = ax1.plot(time_plot, [i*m for i in normal], color + "-")

    # ax2 = ax1.twinx()
    # color = "g"
    # ax2.set_ylabel("$y$-posisjon [m]", color=color)
    # l2, = ax2.plot(time_plot, pos_y, color + "-.")

    # fig.tight_layout()
    # l3 = plt.axvline(x=time_plot[pos_y.index(min(pos_y))], color="k", linewidth=0.7)

    # plt.legend([l1, l5, l2, l3], ["Friksjonskraft", "Normalkraft", "y-posisjon", "Bunnpunkt"], loc="lower left")

    ###############################################################

    # # graf til martin, krumning og g-kraft
    # plt.xlabel("$x$-posisjon [m]")
    # upper_bound = 920
    # a, b = -0.2, 3.2
    # ax1 = fig.add_subplot(111)
    # color = "r"
    # ax1.set_ylabel("Krumning [1/m]", color=color)
    # ax1.set_ylim(a, b)
    # l1, = ax1.plot(pos_x[:upper_bound], curvature[:upper_bound], color + "-")

    # ax2 = ax1.twinx()
    # color = "b"
    # ax2.set_ylabel("G-kraft [1]", color=color)
    # ax2.set_ylim(a, b)
    # l2, = ax2.plot(pos_x[:upper_bound], normal[:upper_bound], color + "-")

    # fig.tight_layout()
    # l3 = plt.axvline(x=pos_x[pos_y.index(min(pos_y))], color="k", linewidth=0.7)

    # plt.legend([l1, l2, l3], ["Krumning", "G-kraft", "Bunnpunkt"], loc="lower center")

    ##############################################################

    # ans = input("Save svg figure? (y/N) ").lower()
    ans = "y"
    if ans == "y":
        plt.savefig(filename + ".svg", format="svg")
    plt.show()


if __name__ == '__main__':
    main()