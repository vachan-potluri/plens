import argparse
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams["font.size"] = 10
# plt.rcParams["figure.figsize"] = 7,6

## Some important notes
# cp is p/rho u^2
# Although the paper says cp is taken as p/(0.5 rho u^2), the values don't match when this
# expression is used. I have checked the freestream pressure and verified that the cp formula
# is missing a factor 2

# Also, the paper plots quantities vs wall coordinate, instead of x-coordinate

p_inf = 9.99473
rho_inf = 0.0004827
u_inf = 2401.824
L = 0.438912
wedge_angle = 18*np.pi/180

parser = argparse.ArgumentParser(
    description = "A script to plot extracted data for Hung wedge 18 case."
)
parser.add_argument("res_dir", help="Result directory (absolute or relative).")
parser.add_argument("counter", help="Output counter.", type=int)
args = parser.parse_args()

res_dir = args.res_dir
if res_dir[-1] != "/":
    res_dir += "/"
counter = args.counter
plate_data = np.genfromtxt(
    "{}plate_data_{}.csv".format(res_dir, counter), delimiter=",", names=True
)
ramp_data = np.genfromtxt(
    "{}ramp_data_{}.csv".format(res_dir, counter), delimiter=",", names=True
)

holden_cp = np.genfromtxt("../data/holden_cp.csv", delimiter=",")
hung_cp = np.genfromtxt("../data/hung_cp.csv", delimiter=",")
# harvey_St = np.genfromtxt("../data/St_harvey.csv", delimiter=",")
# candler_St = np.genfromtxt("../data/St_gnoffo.csv", delimiter=",")

x_plate = plate_data["Points0"]/L
x_ramp = ramp_data["Points0"]/L
s_wall = np.concatenate([x_plate, 1+(x_ramp-1)/np.cos(wedge_angle)])

cp_plate = 2*(plate_data["p"])/(rho_inf*u_inf**2)
cp_ramp = 2*(ramp_data["p"])/(rho_inf*u_inf**2)
cp = np.concatenate([cp_plate, cp_ramp])
cp_mask = np.isfinite(cp)

# St_plate = -2*plate_data["qy"]/(rho_inf*u_inf**3)
# St_ramp = -2*(ramp_data["qy"]*np.cos(np.pi/6) - ramp_data["qx"]*np.sin(np.pi/6))/(rho_inf*u_inf**3)
# St = np.concatenate([St_plate, St_ramp])
# St_mask = np.isfinite(St)

# txy_plate = plate_data["txy"]
# txy_ramp = ramp_data["txy"]
# txy = np.concatenate([txy_plate, txy_ramp])
# txy_mask = np.isfinite(txy)

fig, ax = plt.subplots(1,1)
ax.plot(holden_cp[:,0], 2*holden_cp[:,1], "bo", markersize=4, label="Holden \& Moselle\n(1970, experiment)")
ax.plot(hung_cp[:,0], 2*hung_cp[:,1], "gx", markersize=4, label="Hung \& MacCormack\n(1976, simulation)")
ax.plot(s_wall[cp_mask], cp[cp_mask], "r-", label="PLENS")
ax.set_xlabel(r"$s_{\textrm{wall}}/L$")
ax.set_ylabel(r"$\displaystyle\frac{p}{\frac{1}{2}\rho_\infty u_\infty^2}$", rotation=0, labelpad=20)
ax.grid()
ax.legend(loc="best")
fig.tight_layout()
for fmt in ["png", "pdf"]:
    full_figname = "{}cp_comparison_{}.{}".format(res_dir, counter, fmt)
    fig.savefig(full_figname, format=fmt)
    print("Written file {}".format(full_figname))
plt.show()

# del fig, ax
# fig, ax = plt.subplots(1,1)
# ax.plot(harvey_St[:,0], harvey_St[:,1], "bo", markersize=4, label="Harvey et al (2001)\nexperiment")
# ax.plot(candler_St[:,0], candler_St[:,1], "gx", markersize=4, label="Gnoffo (2001)\nsimulation")
# ax.plot(x[St_mask], St[St_mask], "r-", label="PLENS")
# ax.set_xlabel(r"$x/L$")
# ax.set_ylabel(r"$\displaystyle\frac{2q^{\prime\prime}_w}{\rho_\infty u_\infty^3}$", rotation=0, labelpad=20)
# ax.grid()
# ax.legend(loc="best")
# fig.tight_layout()
# for fmt in ["png", "pdf"]:
#     full_figname = "{}St_comparison_{}.{}".format(res_dir, counter, fmt)
#     fig.savefig(full_figname, format=fmt)
#     print("Written file {}".format(full_figname))
# plt.show()



"""
#del fig, ax
fig, ax = plt.subplots(1,1)
ax.plot(x[txy_mask], txy[txy_mask], "r-", label="PLENS")
# Gnoffo's results on unadaptive 272 by 96 grid
x_sep_gnoffo = 0.52
x_att_gnoffo = 1.31
txy_max = np.max(txy[txy_mask])
txy_min = np.min(txy[txy_mask])
txy_range = np.array([txy_min, -txy_min])
ax.plot(x_sep_gnoffo*np.ones(2), txy_range, "g--")
ax.plot(x_att_gnoffo*np.ones(2), txy_range, "g--", label="Gnoffo (2001) separation\n& attachment locations")
ax.legend()
ax.minorticks_on()
ax.grid(which="major", linestyle='-', linewidth='0.5')
ax.grid(which="minor", linestyle='-', linewidth='0.25')
ax.set_xlabel(r"$x/L$")
ax.set_ylabel(r"$\tau_{xy}$ [Pa]")
fig.tight_layout()
fig.savefig("{}txy_comparison_{}.png".format(res_dir, counter), format="png")
fig.savefig("{}txy_comparison_{}.pdf".format(res_dir, counter), format="pdf")
plt.show()
"""
