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


def main():
    # laste inn et og et forsøk
    # finne normalkraft i bunnpunkt og large i array
    # finne std avvik og gjennomsnitt

    g_force_at_min_array = []
    curvature_at_min_array = []
    alpha_at_min_array = []

    data_path = r"C:/Users/eirik/googledrive/1_Skole/1_Universitetet/_3_Semester/TFY4115_Fysikk/Lab/Dag-3/data/"

    for f in os.listdir(data_path):
        filepath = data_path + f
        data = np.loadtxt(filepath, skiprows=2)  # t x y v
        coeffs = iptrack.iptrack_data(data)
        track = lambda x: bane(x, coeffs)  # banen

        # finne bunnpunkt
        y_data = data[:,2]
        y_min = np.argmin(y_data)
        x_min = data[:,1][y_min]

        curvature_exp = krumning(track, data[:,1])
        speed_exp = data[:,3]
        normal_g_force_exp = 1 + speed_exp**2 * curvature_exp / g

        g_force_at_min_array.append(normal_g_force_exp[y_min])
        curvature_at_min_array.append(1 / curvature_exp[y_min])
        alpha_at_min_array.append(alpha(track, x_min) * 180 / 3.1415926)


    s = np.sqrt(data[:,1]**2 + data[:,2]**2)

    print("g_force_at_min_array")
    print("Avg:", average(g_force_at_min_array))
    print("Std.avvik:", standard_deviation(g_force_at_min_array))
    print("Std.feil:", standard_error(g_force_at_min_array))

    print("\ncurvature_at_min_array")
    print("Avg:", average(curvature_at_min_array))
    print("Std.avvik:", standard_deviation(curvature_at_min_array))
    print("Std.feil:", standard_error(curvature_at_min_array))

    print("\nalpha_at_min_array")
    print("Avg:", average(alpha_at_min_array))
    print("Std.avvik:", standard_deviation(alpha_at_min_array))
    print("Std.feil:", standard_error(alpha_at_min_array))

    fig = plt.figure()
    plt.plot(data[:,1], data[:,2])
    plt.plot(data[:,0], s)
    # plt.plot(data)
    # plt.plot(data[:,1], curvature_exp)
    plt.show()


if __name__ == '__main__':
    main()