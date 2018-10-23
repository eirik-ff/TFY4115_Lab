import numpy as np
import iptrack
import math
import os
from simulering import krumning, normal_g_force, bane, alpha, accel_at_x, x_comp
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
    g_force_at_max_array = []
    curvature_at_min_array = []
    curvature_at_max_array = []
    alpha_at_min_array = []
    alpha_at_max_array = []
    speed_at_min_array = []
    speed_at_max_array = []
    g_force_at_min_sim_array = []
    g_force_at_max_sim_array = []

    data_path = r"C:/Users/eirik/googledrive/1_Skole/1_Universitetet/_3_Semester/TFY4115_Fysikk/Lab/Dag-3/data/"

    for f in os.listdir(data_path):
        filepath = data_path + f
        data = np.loadtxt(filepath, skiprows=2)  # t x y v
        coeffs = iptrack.iptrack_data(data)
        track = lambda x: bane(x, coeffs)  # banen

        ## simulering
        step_h = 0.001  # s, eulers method step
        LOWER_T = 0
        UPPER_T = 1.1  # s, how long to simulate

        # find upper time limit
        with open(filepath, "r") as q:
            lines = q.readlines()
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
        speed_tangential = [0.0] * len(time_)
        accel_tangential = [0.0] * len(time_)
        normal = [0.0] * len(time_)

        # set initial values
        pos_x[0] = data[:,1][0]
        angle[0] = alpha(track, pos_x[0])
        speed_tangential[0] = data[:,3][0]

        for n in range(len(time_) - 1):
            angle[n+1] = alpha(track, pos_x[n])
            accel_tangential[n+1] = accel_at_x(angle[n])
            speed_tangential[n+1] = speed_tangential[n] + step_h * accel_tangential[n]
            pos_x[n+1] = pos_x[n] + step_h * x_comp(speed_tangential[n], angle[n])

            normal[n+1] = normal_g_force(speed_tangential[n], krumning(track, pos_x[n]), angle[n])

        # finne bunnpunkt
        y_data = data[:,2]
        y_min_index = np.argmin(y_data)
        x_min_val = data[:,1][y_min_index]
        t_min_val = data[:,0][y_min_index]

        curvature_exp = krumning(track, data[:,1])
        speed_exp = data[:,3]
        # normal_g_force_exp = 1 + speed_exp**2 * curvature_exp / g
        angles_exp = alpha(track, data[:,1])
        normal_g_force_exp = normal_g_force(speed_exp, curvature_exp, angles_exp)

        y_max_index = np.argmax(normal_g_force_exp)
        x_max_val = data[:,1][y_max_index]
        t_max_val = data[:,0][y_max_index]

        g_force_at_min_array.append(normal_g_force_exp[y_min_index])
        g_force_at_max_array.append(normal_g_force_exp[y_max_index])
        curvature_at_min_array.append(1 / curvature_exp[y_min_index])
        curvature_at_max_array.append(1 / curvature_exp[y_max_index])
        alpha_at_min_array.append(alpha(track, x_min_val) * 180 / 3.1415926)
        alpha_at_max_array.append(alpha(track, x_max_val) * 180 / 3.1415926)
        speed_at_min_array.append(data[:,3][y_min_index])
        speed_at_max_array.append(data[:,3][y_max_index])
        g_force_at_min_sim_array.append(normal[int(t_min_val / step_h)])
        g_force_at_max_sim_array.append(normal[int(t_max_val / step_h)])


    print("G-kraft i bunn (eksp.) [1]")
    print("Avg:", average(g_force_at_min_array))
    print("Std.avvik:", standard_deviation(g_force_at_min_array))
    print("Std.feil:", standard_error(g_force_at_min_array))

    print("\nG-kraft i bunn (sim) [1]")
    print("Avg:", average(g_force_at_min_sim_array))
    print("Std.avvik:", standard_deviation(g_force_at_min_sim_array))
    print("Std.feil:", standard_error(g_force_at_min_sim_array))


    print("\n\nG-kraft i maksimal (eksp) [1]")
    print("Avg:", average(g_force_at_max_array))
    print("Std.avvik:", standard_deviation(g_force_at_max_array))
    print("Std.feil:", standard_error(g_force_at_max_array))

    print("\nG-kraft i maksimal (sim) [1]")
    print("Avg:", average(g_force_at_max_sim_array))
    print("Std.avvik:", standard_deviation(g_force_at_max_sim_array))
    print("Std.feil:", standard_error(g_force_at_max_sim_array))


    print("\n\nKrumningsradius i bunn [m]")
    print("Avg:", average(curvature_at_min_array))
    print("Std.avvik:", standard_deviation(curvature_at_min_array))
    print("Std.feil:", standard_error(curvature_at_min_array))

    print("\nKrumningsradius i maksimal [m]")
    print("Avg:", average(curvature_at_max_array))
    print("Std.avvik:", standard_deviation(curvature_at_max_array))
    print("Std.feil:", standard_error(curvature_at_max_array))


    print("\n\nVinkel i bunn []")
    print("Avg:", average(alpha_at_min_array))
    print("Std.avvik:", standard_deviation(alpha_at_min_array))
    print("Std.feil:", standard_error(alpha_at_min_array))

    print("\nVinkel i maksimal [*]")
    print("Avg:", average(alpha_at_max_array))
    print("Std.avvik:", standard_deviation(alpha_at_max_array))
    print("Std.feil:", standard_error(alpha_at_max_array))


    print("\n\nFart i bunn [m/s]")
    print("Avg:", average(speed_at_min_array))
    print("Std.avvik:", standard_deviation(speed_at_min_array))
    print("Std.feil:", standard_error(speed_at_min_array))

    print("\nFart i maksimal [m/s]")
    print("Avg:", average(speed_at_max_array))
    print("Std.avvik:", standard_deviation(speed_at_max_array))
    print("Std.feil:", standard_error(speed_at_max_array))



if __name__ == '__main__':
    main()