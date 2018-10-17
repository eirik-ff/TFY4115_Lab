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
import sys

# returns the coefficients of a polynomial of degree 15
def iptrack(filename):
    data=np.loadtxt(filename,skiprows=2)
    return np.polyfit(data[:,1],data[:,2],15)


def iptrack_data(data):
    return np.polyfit(data[:,1],data[:,2],15)


def print_expr(coeffs, var="x", py=False):
    """
    coeffs: list of coeffs from iptrack function
    var: variable in expression
    py: python friendly with * and **
    """
    term = "+ {0} {1}^{2} "
    if py:
        term = "+ {0}*{1}**{2} "

    # generates the polynomial expression
    poly = ""
    for n, c in enumerate(coeffs):
        poly += term.format(c, var, 15 - n)

    # prints the polynomial excluding the first + sign
    print(poly[1:])


def main():
    data = "test2.txt"
    coeffs = iptrack(data)
    print_expr(coeffs, py=True)    


if __name__ == '__main__':
    main()
