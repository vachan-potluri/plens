# This script sets-up a case. This script is generally executed after 'generate_geo_file.py'
# 1. Create mesh using gmsh (assuming '.geo' file is available and 'gmsh' is an executable). The
#    mesh is created in 'msh' format.
# 2. Write a file for IC
# 3. Write the input file
# 4. Create result directory

import argparse
import os
import subprocess

parser = argparse.ArgumentParser(
    description = "A script to setup a case. Assumes '.geo' file is available for mesh generation "
        + "using gmsh and 'gmsh' is an executable."
)
parser.add_argument("dest", help="Location where setup is being done")
parser.add_argument(
    "-g",
    "--geo_file",
    help="The name of gmsh 'geo' script to generate mesh (relative to the location given in the "
    + "'dest' argument, '.geo' extension is optional). Default: 'mesh.geo'.",
    default="mesh.geo",
    action="store"
)
parser.add_argument(
    "-f",
    "--flux_scheme",
    help="Flux scheme. Default: 'HLLC'.",
    default="HLLC",
    action="store"
)
parser.add_argument(
    "-k",
    "--rk_order",
    help="RK integration order. Default: 3.",
    type=int,
    default=3,
    action="store"
)
parser.add_argument(
    "-c",
    "--cfl",
    help="CFL number. Default: 0.5.",
    type=float,
    default=0.5,
    action="store"
)
parser.add_argument(
    "-w",
    "--write_frequency",
    help="Write frequency. Default: 100.",
    type=float,
    default=100,
    action="store"
)
parser.add_argument(
    "-r",
    "--res_dir",
    help="Name of the result directory (relative to the location given in the 'dest' argument). "
        + "Default: 'result'",
    default="result",
    action="store"
)
args = parser.parse_args()



# 1. Generate mesh
dest = args.dest
if dest[-1] != "/":
    dest += "/"

geo_file = args.geo_file
if geo_file[-4:] != ".geo":
    geo_file += ".geo"

full_mesh_filename = dest + geo_file

# https://stackoverflow.com/questions/12060863/python-subprocess-call-a-bash-alias
command = "gmsh {} -3 -format msh4".format(full_mesh_filename)
subprocess.run(["/bin/bash", "-i", "-c", command])



# 2. IC file
ic_content = """p
2
1
1
-4


3.857143
2.629369
0
0
10.33333
1 + 0.2*sin(5*x)
0
0
0
1
"""
ic_filename = dest + "ic_shu_osher.dat"
ic_file = open(ic_filename, "w")
ic_file.write(ic_content)
ic_file.close()
print("Written IC in file {}".format(ic_filename))



# 3. Input file
mesh_filename = geo_file[:-4] + ".msh"
input_content = """
subsection mesh
	set type = straight
	set format = msh
	set file name = {}
end

subsection Navier-Stokes
	set gas name = air
	set inviscid = true
	set inviscid surface flux scheme = {}
end

subsection IC
	set type = piecewise function
	set file name = ic_shu_osher.dat
end

subsection BCs
	subsection bid1
		set type = uniform inflow
        set prescribed p = 10.33333
        set prescribed velocity = 2.629369 0 0
        set prescribed T = 0.009344639 
	end
	subsection bid2
	    set type = free
	end
	subsection bid3
		set type = symmetry
	end
end

subsection blender parameters
    set variable = pxrho
end

subsection time integration
    set RK order = {}
    set Courant number = {}
    set end time = 1.8
end

subsection data output
    set write frequency = {}
    set directory = {}
end
""".format(
    mesh_filename,
    args.flux_scheme,
    args.rk_order,
    args.cfl,
    args.write_frequency,
    args.res_dir
)
input_filename = dest + "input.prm"
input_file = open(input_filename, "w")
input_file.write(input_content)
input_file.close()
print("Written input file in {}".format(input_filename))



# 4. Create result directory
full_res_dir = dest + args.res_dir
if os.path.isdir(full_res_dir) == False: subprocess.run(["mkdir", full_res_dir])
