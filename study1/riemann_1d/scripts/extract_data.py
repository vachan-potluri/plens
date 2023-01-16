# A script to extract line data from pvtu files for Riemann 1d cases

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
parser.add_argument(
    "-r",
    "--resolution",
    help="Resolution of line source. Default: 6400.",
    type=int,
    default=6400,
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
parser.add_argument(
    "-e",
    "--extent",
    help="The domain extent in x direction. Default: '0 1'.",
    type=float,
    nargs=2,
    default=[0.0, 1.0],
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
    Input=result
)

plotOverLine1.Point1 = [args.extent[0], 0.25, 0.25]
plotOverLine1.Point2 = [args.extent[1], 0.25, 0.25]
plotOverLine1.Resolution = args.resolution

# save data
full_data_filename = res_dir + args.base_data_filename + "_{:06d}.csv".format(args.counter)
SaveData(
    full_data_filename,
    proxy=plotOverLine1
)
print("Written data file {}".format(full_data_filename))
