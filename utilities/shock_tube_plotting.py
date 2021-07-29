# directory for the shock tube test
test_dir = "../test1-2/"

# directory where results are stored within the test directory
res_dir = "result/"

# name of file containing the line plot data from paraview
dfname = "line_data.csv"

# Test data
# access: data[0/1][test id]
rho_ic = [ [1.0,0.125], [1.0,1.0], [1.0,1.0], [5.99924,5.99242], [1.0,1.0] ]
p_ic = [ [1.0,0.1], [0.4,0.4], [1000.0,0.01], [460.894,46.095], [1000,0.01] ]
u_ic = [ [0.75,0.0], [-2.0,2.0], [0.0,0.0], [19.5975,-6.19633], [-19.59745,-19.59745] ]
dia_locs = [0.5, 0.5, 0.5, 0.4, 0.8]
end_times = [0.2, 0.15, 0.012, 0.035, 0.012]

# test settings
R = 8.314462618/0.029
gamma = 1.4
test_id = 2-1
rhol = rho_ic[test_id][0]
rhor = rho_ic[test_id][1]
pl = p_ic[test_id][0]
pr = p_ic[test_id][1]
Tl = pl/(rhol*R)
Tr = pr/(rhor*R)
ul = u_ic[test_id][0]
ur = u_ic[test_id][1]
time = end_times[test_id]
dia_loc = dia_locs[test_id] # location of initial diaphragm

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 14
plt.rcParams["mathtext.fontset"] = "dejavuserif"
from shockTubeSoln import *

# IMP: change the column heading of "Points:0" to "x"
# python doesn't accept colons in heading
fulldfname = test_dir + res_dir + dfname
print("Data file: {}".format(fulldfname))
num_data = np.genfromtxt(fulldfname, delimiter=",", names=True)
x = num_data["x"]

ex_data = shockTubeSoln(pl,pr,Tl,Tr,ul,ur,x-dia_loc,time)

## Plotting
fig, axes = plt.subplots(2,2)
formats = ["b-", "ro"] # exact and numerical

l1, = axes[0,0].plot(x,ex_data[:,4],formats[0])
l2, = axes[0,0].plot(x,num_data["rho"],formats[1],markersize=2)
axes[0,0].set_ylabel(r"$\rho$")

num_p = (gamma-1)*( num_data["rhoE"] - 0.5*(
        num_data["rhou"]**2/num_data["rho"] +
        num_data["rhov"]**2/num_data["rho"] +
        num_data["rhow"]**2/num_data["rho"]
) )
axes[1,0].plot(x,ex_data[:,2],formats[0])
axes[1,0].plot(x,num_p,formats[1],markersize=2)
axes[1,0].set_ylabel(r"$p$")

num_u = num_data["rhou"]/num_data["rho"]
axes[0,1].plot(x,ex_data[:,3],formats[0])
axes[0,1].plot(x,num_u,formats[1],markersize=2)
axes[0,1].set_ylabel(r"$u$")

for ax in axes.reshape(-1):
        ax.grid()
        ax.set_xlabel(r"$x$")
axes[1,1].remove()

fig.legend([l1,l2], ["exact", "simulation"], loc="center", bbox_to_anchor=[0.75,0.25])
fig.tight_layout(rect=[0,0,1,1])

figname = test_dir + res_dir + "comparison.pdf"
print("Figure name: {}".format(figname))
fig.savefig(figname, format="pdf")
