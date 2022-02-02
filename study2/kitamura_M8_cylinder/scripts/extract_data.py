# A script to extract line data from pvtu files for Kitamura's Mach 8 flow over cylinder

import argparse
from paraview.simple import *

parser = argparse.ArgumentParser(
    description = "A script to extract line data from pvtu files for Kitamura's case. This script "
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

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=result)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]
slice1.SliceType.Origin = [-0.02758819074369967, 0.0, 0.009999999776482582]
slice1.HyperTreeGridSlicer.Origin = [-0.02758819074369967, 0.0, 0.009999999776482582]
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# create a new 'Plot On Intersection Curves'
plotOnIntersectionCurves1 = PlotOnIntersectionCurves(
    registrationName='PlotOnIntersectionCurves1',
    Input=slice1
)
plotOnIntersectionCurves1.SliceType = 'Sphere'
plotOnIntersectionCurves1.SliceType.Center = [0.0, 0.0, 0.01]
plotOnIntersectionCurves1.SliceType.Radius = 0.020000001

# save data
temp_filename = '{}surface_data_{}.csv'.format(res_dir, counter)
SaveData(
    temp_filename,
    proxy=plotOnIntersectionCurves1,
    PointDataArrays=['Subdomain', 'T', 'alpha', 'arc_length', 'k', 'loc_dt', 'mu', 'p', 'qx', 'qy', 'qz', 'rho', 'rhoE', 'rhou', 'rhov', 'rhow', 'steady_state_error', 'txx', 'txy', 'txz', 'tyy', 'tyz', 'tzz', 'u', 'v', 'w'],
    UseScientificNotation=1
)
print("Written file {}".format(temp_filename))
