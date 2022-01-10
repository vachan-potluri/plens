# A script to extract line data from pvtu files for Shu-Osher case

import argparse
from paraview.simple import *

parser = argparse.ArgumentParser(
    description = "A script to extract line data from pvtu files for Shu-Osher case. This script "
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
parser.add_argument(
    "-d",
    "--base_data_filename",
    help="Base data file name. Default: 'line_data'. Data is stored as csv file with 'counter' "
        + "appended to the filename. The file is stored in 'res_dir'.",
    default="line_data",
    action="store"
)
args = parser.parse_args()

res_dir = args.res_dir
if res_dir[-1] != "/":
    res_dir += "/"

full_op_filename = res_dir + args.base_output_filename + "_{:06d}.pvtu".format(args.counter)
print("Reading output file {}".format(full_op_filename))
result = XMLPartitionedUnstructuredGridReader(FileName=full_op_filename)

# create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(
    registrationName='PlotOverLine1',
    Input=result,
    Source='Line'
)

plotOverLine1.Source.Point1 = [-5.0, 0.025, 0.025]
plotOverLine1.Source.Point2 = [5.0, 0.025, 0.025]

# save data
full_data_filename = res_dir + args.base_data_filename + "_{:06d}.csv".format(args.counter)
SaveData(
    full_data_filename,
    proxy=plotOverLine1,
    PointDataArrays=['Subdomain', 'T', 'alpha', 'arc_length', 'k', 'mu', 'p', 'qx', 'qy', 'qz', 'rho', 'rhoE', 'rhou', 'rhov', 'rhow', 'steady_state_error', 'txx', 'txy', 'txz', 'tyy', 'tyz', 'tzz', 'u', 'v', 'vtkValidPointMask', 'w']
)
print("Written data file {}".format(full_data_filename))
