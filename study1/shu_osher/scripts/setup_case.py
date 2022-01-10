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
        + "using gmsh."
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
    help="Name of the result directory. Default: 'result'",
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

subprocess.run(["gmsh", full_mesh_filename, "-3", "-format", "msh4"], env=dict(os.environ), shell=True)
