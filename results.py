import numpy as np
import iptrack
import math
import os
from simulering import krumning, normal_g_force, bane, alpha
import matplotlib
import matplotlib.pyplot as plt


g = 9.81  # m/s^2


def parse_data(filepath):
    return np.loadtxt(filepath, skiprows=2)


def get_minimum_point_index(data):
    """
    assumes data is: t x y v
    """
    return np.argmin(data[:,2])


def average(a):
    return sum(a) / len(a)


def standard_deviation(a):
    avg = average(a)
    numer = sum((x - avg)**2 for x in a)
    return math.sqrt(numer / (len(a) - 1))


def standard_error(a):
    s = standard_deviation(a)
    return s / math.sqrt(len(a))


def lpf(x, alpha, x0=None):
    """
    low pass filter, see https://www.embeddedrelated.com/showarticle/779.php 
    """
    y = [0] * len(x)
    yk = x[0] if x0 is None else x0
    for k in range(len(y)):
        yk += alpha * (x[k] - yk)
        y[k] = yk
    return y


def main():
    # laste inn et og et forsøk
    # finne normalkraft i bunnpunkt og large i array
    # finne std avvik og gjennomsnitt

    g_force_at_min_array = []
    curvature_at_min_array = []
    alpha_at_min_array = []
    speed_at_min_array = []

    data_path = r"C:/Users/eirik/googledrive/1_Skole/1_Universitetet/_3_Semester/TFY4115_Fysikk/Lab/Dag-3/data/"

    for f in os.listdir(data_path):
        filepath = data_path + f
        data = np.loadtxt(filepath, skiprows=2)  # t x y v
        coeffs = iptrack.iptrack_data(data)
        track = lambda x: bane(x, coeffs)  # banen

        # finne bunnpunkt
        y_data = data[:,2]
        y_min_index = np.argmin(y_data)
        x_min_val = data[:,1][y_min_index]

        curvature_exp = krumning(track, data[:,1])
        speed_exp = data[:,3]
        normal_g_force_exp = 1 + speed_exp**2 * curvature_exp / g

        g_force_at_min_array.append(normal_g_force_exp[y_min_index])
        curvature_at_min_array.append(1 / curvature_exp[y_min_index])
        alpha_at_min_array.append(alpha(track, x_min_val) * 180 / 3.1415926)
        speed_at_min_array.append(data[:,3][y_min_index])


    print("G-kraft i bunn [1]")
    print("Avg:", average(g_force_at_min_array))
    print("Std.avvik:", standard_deviation(g_force_at_min_array))
    print("Std.feil:", standard_error(g_force_at_min_array))

    print("\nKrumningsradius i bunn [m]")
    print("Avg:", average(curvature_at_min_array))
    print("Std.avvik:", standard_deviation(curvature_at_min_array))
    print("Std.feil:", standard_error(curvature_at_min_array))

    print("\nVinkel i bunn [*]")
    print("Avg:", average(alpha_at_min_array))
    print("Std.avvik:", standard_deviation(alpha_at_min_array))
    print("Std.feil:", standard_error(alpha_at_min_array))

    print("\nFart i bunn [m/s]")
    print("Avg:", average(speed_at_min_array))
    print("Std.avvik:", standard_deviation(speed_at_min_array))
    print("Std.feil:", standard_error(speed_at_min_array))

    fig = plt.figure()
    plt.plot(data[:,0], data[:,3])
    for t in [0.01, 0.05, 0.1, 0.2, 0.5]:
        plt.plot(data[:,0][1:], lpf(np.diff(data[:,3])/np.diff(data[:,0]), (data[:,0][1] - data[:,0][0])/t))
    # plt.plot(data[:,0][1:], np.diff(data[:,3])/np.diff(data[:,0]))
    # plt.plot(data)
    # plt.plot(data[:,1], curvature_exp)
    legends = tuple(str(i) for i in [0.01, 0.05, 0.1, 0.2, 0.5])
    plt.legend(("fart",) + legends, loc="best")
    plt.show()


if __name__ == '__main__':
    main()