# This file converts gmsh transfinite specifications into different formats. The gmsh specification
# for a given edge (and hence of fixed length) is based on the number of elements and a geometric
# progression ratio. Seldom is this useful. Many a times, one of the following is more useful
# - Start and end element size
# - Start element size and number of elements
#
# This small script provides forward and inverse functions for these conversions. The gmsh
# specification will be called format 0, start/end size as format 1 and start size along with
# number of elements as format 2.

import numpy as np
from numpy.polynomial import Polynomial

def format_0to1(length, n_cells, geom_ratio):
    # Output start size and end size
    # Total length is the sum of geometric progression
    # Use that to calculate starting value, that will be the min if ratio > 1 and max otherwise
    if geom_ratio == 1:
        start = length/n_cells
    else:
        start = length*(geom_ratio-1)/(geom_ratio**n_cells-1)
    end = start*geom_ratio**(n_cells-1)
    print("Given:\n\tLength {}\n\tNumber of cells {}\n\tGeometric ratio {}".format(
        length, n_cells, geom_ratio
    ))
    print("Output:\n\tStart size {}\n\tEnd size {}".format(start, end))
    return (start, end)

def format_1to0(length, start_size, end_size):
    # This is the inverse of format_0to1
    # The ratio of sum of geometric progression and the starting term is a function of the ratio
    # of the end to start value and the geometric ratio. Use this to get r
    # S/a = \frac{r^{n-1}-1/r}{1-1/r}
    # where r^{n-1} is the ratio of end size to start size
    geom_ratio = (start_size - length)/(end_size - length)
    if geom_ratio == 1:
        n_cells = length/start_size
    else:
        n_cells = np.log(end_size/start_size)/np.log(geom_ratio) + 1
    print("Given:\n\tLength {}\n\tStart size {}\n\tEnd size {}".format(
        length, start_size, end_size
    ))
    print("Output:\n\tNumber of cells {}\n\tGeometric ratio {}".format(n_cells, geom_ratio))
    return (n_cells, geom_ratio)

def format_0to2(length, n_cells, geom_ratio):
    # Output size start size and n_cells
    # n_cells is already available, use format_0to1 for start size
    format_0to1(length, n_cells, geom_ratio)

def format_2to0(length, start_size, n_cells):
    # inverse of format_0to2
    # S/a = \frac{r^n - 1}{r-1}
    # S, a and n are known, r is required. A nonlinear polynomial equation.
    coeffs = np.zeros(n_cells+1)
    coeffs[-1] = 1 # r^n
    coeffs[0] = -1+length/start_size # constant
    coeffs[1] = -length/start_size # r^1
    # https://numpy.org/doc/stable/reference/routines.polynomials.html
    p = Polynomial(coeffs)
    r = p.roots()
    real_roots = r.real[abs(r.imag)<1e-4]
    real_pos_roots = real_roots[real_roots>0]
    print("Given:\n\tLength {}\n\tStart size {}\n\tNumber of cells {}".format(
        length, start_size, n_cells
    ))
    print("Real positive roots for geometric ratio:\n\t{}".format(real_pos_roots))
    print("Scale ratios (geom_ratio^(n_cells-1)):\n\t{}".format(real_pos_roots**(n_cells-1)))

format_2to0(40e-3, 1e-6, 80)

