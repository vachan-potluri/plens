import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
    "font.size": 12,
    "axes.formatter.limits": [-2,2]
})
import argparse

parser = argparse.ArgumentParser(
    description = "A script to plot extracted data for Deepak et al (2013) Mach 9 wedge 18 case."
)
parser.add_argument("res_dir", help="Result directory (absolute or relative).")
parser.add_argument("counter", help="Output counter.", type=int)
args = parser.parse_args()

res_dir = args.res_dir
if res_dir[-1] != "/":
    res_dir += "/"
counter = args.counter

p_inf = 730
u_inf = 2280
rho_inf = 16e-3
T_inf = 160
T_w = 300
Pr = 0.69 # Prandtl number
cp = 1.4*286.7/(1.4-1) # specific heat
L = 8.5e-2 # plate length
theta = 18*np.pi/180 # wedge angle
plate_data = np.genfromtxt(
    "{}plate_data_{}.csv".format(res_dir, counter), delimiter=",", names=True
)
ramp_data = np.genfromtxt(
    "{}ramp_data_{}.csv".format(res_dir, counter), delimiter=",", names=True
)
deepak_p = np.genfromtxt("../data/deepak_p.csv", delimiter=",") # s/L, p/p_inf
olejniczak_p = np.genfromtxt("../data/olejniczak_p.csv", delimiter=",")
mallinson_p = np.genfromtxt("../data/mallinson_p.csv", delimiter=",")
deepak_St = np.genfromtxt("../data/deepak_St.csv", delimiter=",") # s/L, St
olejniczak_St = np.genfromtxt("../data/olejniczak_St.csv", delimiter=",")
mallinson_St = np.genfromtxt("../data/mallinson_St.csv", delimiter=",")

x_plate = plate_data["Points0"]/L
x_ramp = ramp_data["Points0"]/L
s_ramp = 1 + (x_ramp-1)/np.cos(theta)
x = np.concatenate([x_plate, x_ramp])
s = np.concatenate([x_plate, s_ramp])

p_plate = plate_data["p"]
p_ramp = ramp_data["p"]
p = np.concatenate([p_plate, p_ramp])
p_mask = np.isfinite(p)

q_plate = -plate_data["qy"]
q_ramp = -ramp_data["qy"]*np.cos(theta) + ramp_data["qx"]*np.sin(theta)
q = np.concatenate([q_plate, q_ramp])
# See eqs (1-2) in Deepak et al (2013)
St = q/(rho_inf*u_inf*(cp*(T_inf-T_w) + np.sqrt(Pr)*u_inf**2/2))
St_mask = np.isfinite(St)

fig, ax = plt.subplots(1,1)
ax.plot(s[p_mask], p[p_mask]/p_inf, "r-", label="PLENS")
ax.plot(deepak_p[:,0], deepak_p[:,1], "b--", label="Deepak et al\n(2013) simulation")
ax.plot(olejniczak_p[:,0], olejniczak_p[:,1], "g:", label="Olejnicak \& Candler\n(1998) simulation")
ax.plot(mallinson_p[:,0], mallinson_p[:,1], "ko", label="Mallinson (1997)\nexperiment")
ax.set_xlabel(r"$s_\textrm{wall}/L$")
ax.set_ylabel(r"$\displaystyle\frac{p}{p_\infty}$", rotation=0, labelpad=10)
ax.grid()
ax.legend(loc="best")
fig.tight_layout()
plt.show()
for fmt in ["png", "pdf"]:
    temp_filename = "{}p_comparison_{}.{}".format(res_dir, counter, fmt)
    fig.savefig(temp_filename, format=fmt)
    print("Saved figure {}".format(temp_filename))

del fig, ax
fig, ax = plt.subplots(1,1)
ax.plot(s[St_mask], St[St_mask], "r-", label="PLENS")
ax.plot(deepak_St[:,0], deepak_St[:,1], "b--", label="Deepak et al\n(2013) simulation")
ax.plot(olejniczak_St[:,0], olejniczak_St[:,1], "g:", label="Olejnicak \& Candler\n(1998) simulation")
ax.plot(mallinson_St[:,0], mallinson_St[:,1], "ko", label="Mallinson (1997)\nexperiment")
ax.set_xlabel(r"$s_\textrm{wall}/L$")
ax.set_ylabel("St", rotation=0)
ax.grid()
ax.legend(loc="best")
fig.tight_layout()
plt.show()
for fmt in ["png", "pdf"]:
    temp_filename = "{}St_comparison_{}.{}".format(res_dir, counter, fmt)
    fig.savefig(temp_filename, format=fmt)
    print("Saved figure {}".format(temp_filename))
