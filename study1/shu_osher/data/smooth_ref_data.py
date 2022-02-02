import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, splev
from scipy.signal import savgol_filter

ref_data = np.genfromtxt("rho_ref_data_monotonic.txt")
x_ref = ref_data[:,0]
rho_ref = ref_data[:,1]

# bi-variate spline interpolation
# https://stackoverflow.com/questions/14344099/smooth-spline-representation-of-an-arbitrary-contour-flength-x-y
dist_ord = 2 # norm type for distance calculation
t = np.zeros_like(x_ref)
t[1:] = (
    np.abs((x_ref[1:] - x_ref[:-1]))**dist_ord +
    np.abs((rho_ref[1:] - rho_ref[:-1]))**dist_ord
)**(1.0/dist_ord)
t = np.cumsum(t)
t /= t[-1]
nt = np.linspace(0,1,6400)
spl_ord = 2 # spline degree
x_spl = make_interp_spline(t, x_ref, spl_ord)
rho_spl = make_interp_spline(t, rho_ref, spl_ord)
x_interp = splev(nt, x_spl)
rho_interp = splev(nt, rho_spl)

# smoothing data
# https://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way
window_size = 7
pol_ord = 3
rho_smooth = savgol_filter(rho_ref, window_size, pol_ord)

plt.plot(x_ref, rho_ref, "bo-", label="Digitizer", markersize=2)
plt.plot(x_interp, rho_interp, "r-", label="Interpolated, d={}, p={}".format(dist_ord, spl_ord))
plt.plot(x_ref, rho_smooth, "g-", label="Filtered, w={}, p={}".format(window_size, pol_ord))
plt.grid()
plt.legend()
plt.show()
