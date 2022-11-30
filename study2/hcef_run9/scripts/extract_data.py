# A script to extract surface data from pvtu files for HCEF geometry

import argparse
import numpy as np
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
parser.add_argument(
    "-a",
    "--azimuths",
    help="Extract data on multiple azimuthal angles (in degrees). Default: ''. If this argument is "
        + "empty, only the zero-degree azimuthal angle is used for extracting.",
    default=[0.0],
    type=float,
    action="store",
    nargs="*"
)
args = parser.parse_args()



# function to save cylinder and flare data on a given azimuthal angle
def extract_data(result, res_dir, counter, azimuth=0):
    phi = np.pi*azimuth/180
    eps = 1e-6
    r1 = 0.0325 + eps
    r2 = 0.1008 + eps
    x1 = 0.01
    x2 = 0.1017
    x3 = 0.22
    plotOverLine1 = PlotOverLine(
        registrationName='PlotOverLine1',
        Input=result,
    )
    plotOverLine1.Point1 = [x1, r1*np.cos(phi), r1*np.sin(phi)]
    plotOverLine1.Point2 = [x2, r1*np.cos(phi), r1*np.sin(phi)]
    plotOverLine1.Resolution = args.resolution
    temp_filename = '{}cylinder_data_phi{:.0f}_{}.csv'.format(res_dir, azimuth, counter)
    SaveData(
        temp_filename,
        proxy=plotOverLine1,
    )
    print("Written file {}".format(temp_filename))

    plotOverLine2 = PlotOverLine(
        registrationName='PlotOverLine2',
        Input=result,
    )
    plotOverLine2.Point1 = plotOverLine1.Point2
    plotOverLine2.Point2 = [x3, r2*np.cos(phi), r2*np.sin(phi)]
    plotOverLine1.Resolution = args.resolution
    temp_filename = '{}flare_data_phi{:.0f}_{}.csv'.format(res_dir, azimuth, counter)
    SaveData(
        temp_filename,
        proxy=plotOverLine2,
    )
    print("Written file {}".format(temp_filename))



res_dir = args.res_dir
if res_dir[-1] != "/":
    res_dir += "/"
counter = args.counter
full_data_filename = "{}{}_{:06d}.pvtu".format(res_dir, args.base_output_filename, counter)
result = XMLPartitionedUnstructuredGridReader(
    FileName=full_data_filename
)
print("Read file {}".format(full_data_filename))

for azimuth in args.azimuths: extract_data(result, res_dir, counter, azimuth)
