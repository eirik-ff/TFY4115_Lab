# iptrack - interpolate track
#
# SYNTAX
# p=iptrack(filename)
#
# INPUT
# filename: data file containing exported tracking data on the standard
# Tracker export format
#
# mass_A
# t x   y
# 0.0   -1.0686477620876644 42.80071293284619
# 0.04  -0.714777136706708  42.62727536827738
# ...
#
# OUTPUT
# p=iptrack(filename) returns the coefficients of a polynomial of degree 15
# that is the least square fit to the data y(x). Coefficients are given in
# descending powers.

import numpy as np

# returns the coefficients of a polynomial of degree 15
def iptrack(filename):
    data=np.loadtxt(filename,skiprows=2)
    return np.polyfit(data[:,1],data[:,2],15)


def main():
    data = "krumbane.txt"
    coeffs = iptrack(data)

    # generates the polynomial expression
    poly = ""
    for n, c in enumerate(coeffs):
        poly += "+ {} x^{}".format(c, 15 - n)

    # prints the polynomial excluding the first + sign
    print(poly[1:])


if __name__ == '__main__':
    main()