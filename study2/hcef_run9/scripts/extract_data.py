# A script to extract line data from pvtu files for HCEF run 14

import argparse
from paraview.simple import *

parser = argparse.ArgumentParser(
    description = "A script to extract line data from pvtu files for HCEF case. This script "
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
    "-r",
    "--resolution",
    help="Resolution of line source. Default: 100.",
    type=int,
    default=100,
    action="store"
)
args = parser.parse_args()

res_dir = args.res_dir
if res_dir[-1] != "/":
    res_dir += "/"
counter = args.counter
full_data_filename = "{}{}_{:06d}.pvtu".format(res_dir, args.base_output_filename, counter)
result = XMLPartitionedUnstructuredGridReader(
    FileName=full_data_filename
)
print("Read file {}".format(full_data_filename))

# create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(
    registrationName='PlotOverLine1',
    Input=result,
    Source='Line'
)

plotOverLine1.Source.Point1 = [0.01, 0.0325, 1e-8]
plotOverLine1.Source.Point2 = [0.1017, 0.0325, 1e-8]
plotOverLine1.Source.Resolution = args.resolution

# save data
temp_filename = '{}cylinder_data_{}.csv'.format(res_dir, counter)
SaveData(
    temp_filename,
    proxy=plotOverLine1,
    PointDataArrays=['Subdomain', 'T', 'alpha', 'arc_length', 'k', 'mu', 'p', 'qx', 'qy', 'qz', 'rho', 'rhoE', 'rhou', 'rhov', 'rhow', 'steady_state_error', 'txx', 'txy', 'txz', 'tyy', 'tyz', 'tzz', 'u', 'v', 'vtkValidPointMask', 'w']
)
print("Written file {}".format(temp_filename))

# create a new 'Plot Over Line'
plotOverLine2 = PlotOverLine(
    registrationName='PlotOverLine2',
    Input=result,
    Source='Line'
)

plotOverLine2.Source.Point1 = [0.1017, 0.0325, 1e-8]
plotOverLine2.Source.Point2 = [0.22, 0.100801, 1e-8]
plotOverLine1.Source.Resolution = args.resolution

# save data
temp_filename = '{}flare_data_{}.csv'.format(res_dir, counter)
SaveData(
    temp_filename,
    proxy=plotOverLine2,
    PointDataArrays=['Subdomain', 'T', 'alpha', 'arc_length', 'k', 'mu', 'p', 'qx', 'qy', 'qz', 'rho', 'rhoE', 'rhou', 'rhov', 'rhow', 'steady_state_error', 'txx', 'txy', 'txz', 'tyy', 'tyz', 'tzz', 'u', 'v', 'vtkValidPointMask', 'w']
)
print("Written file {}".format(temp_filename))

