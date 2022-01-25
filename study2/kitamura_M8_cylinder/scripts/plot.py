import argparse
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams["font.size"] = 10

parser = argparse.ArgumentParser(
    description = "A script to plot extracted data for Kitamura's case."
)
parser.add_argument("res_dir", help="Result directory (absolute or relative).")
parser.add_argument("counter", help="Output counter.", type=int)
args = parser.parse_args()

res_dir = args.res_dir
if res_dir[-1] != "/":
    res_dir += "/"
counter = args.counter

surface_data = np.genfromtxt(
    "{}surface_data_{}.csv".format(res_dir, counter), delimiter=",", names=True
)

r = 0.02
p_inf = 370.7
p0 = p_inf*84.9 # this ratio is given in Kitamura (2010)
p = surface_data["p"]
theta = np.arcsin(surface_data["Points1"]/r)*180/np.pi

fig, ax = plt.subplots(1,1)
ax.plot(theta, p/p0, "r-")
ax.grid()
ax.set_xlabel(r"$\theta$ [degrees]")
ax.set_ylabel(r"$\frac{p}{p_0}$")
plt.show()
