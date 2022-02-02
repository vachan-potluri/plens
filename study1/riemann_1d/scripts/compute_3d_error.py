# Computes 3d error for a given file as the integration of transverse velocity vector
# This can be scaled by a corresponding x-direction velocity for a normalised error

import argparse
from paraview.simple import *

parser = argparse.ArgumentParser(
    description = "A script to extract line data from pvtu files for shock tube cases. This script "
        + "has to be executed using 'pvpython'."
)
parser.add_argument("res_dir", help="Result directory (absolute or relative).")
parser.add_argument("counter", help="Output counter.", type=int)
parser.add_argument(
    "-b",
    "--base_output_filename",
    help="Base output file name. Default: 'output'.",
    default="output",
    action="store"
)
args = parser.parse_args()

res_dir = args.res_dir
if res_dir[-1] != "/":
    res_dir += "/"

full_op_filename = res_dir + args.base_output_filename + "_{:06d}.pvtu".format(args.counter)
print("Reading output file {}".format(full_op_filename))
result = XMLPartitionedUnstructuredGridReader(FileName=full_op_filename)

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=result)
calculator1.ResultArrayName = 'transverse_vel'
calculator1.Function = 'sqrt(v*v+w*w)'

integrateVariables1 = IntegrateVariables(registrationName='IntegrateVariables1', Input=calculator1)
iv_pdata = paraview.servermanager.Fetch(integrateVariables1).GetPointData()

print(iv_pdata.GetArray("transverse_vel").GetValue(0))
