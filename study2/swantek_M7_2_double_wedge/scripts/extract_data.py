# A script to extract surface data from pvtu files for Swantek's double wedge geometry

import argparse
import numpy as np
from paraview.simple import *

parser = argparse.ArgumentParser(
    description = "A script to extract line data from pvtu files for Swantek's double wedge "
        "geometry. This script has to be executed using 'pvpython'. The geometry is assumed to "
        "be similar to that in Komives et al (2014)."
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
    help="Resolution of line source. Default: 500.",
    type=int,
    default=500,
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

# Geometry parameters
l_sym = 5e-3 # length of horizontal symmetry portion, before wedge1
l_wedge1 = 50.8e-3
l_wedge2 = 25.4e-3
theta_wedge1 = 30*np.pi/180
theta_wedge2 = 55*np.pi/180
w = 5e-3 # width in z direction
wedge1_factor = 0.2 # extraction on wedge 1 will start after this much percentage of the length
                    # this is required because at origin, values will be too high

# create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(
    registrationName='PlotOverLine1',
    Input=result,
    Source='Line'
)
plotOverLine1.Source.Point1 = [
        l_sym + wedge1_factor*l_wedge1*np.cos(theta_wedge1),
        wedge1_factor*l_wedge1*np.sin(theta_wedge1),
        0.5*w
]
plotOverLine1.Source.Point2 = [
        l_sym + l_wedge1*np.cos(theta_wedge1),
        l_wedge1*np.sin(theta_wedge1),
        0.5*w
]
plotOverLine1.Source.Resolution = args.resolution

# save data
temp_filename = '{}wedge1_data_{}.csv'.format(res_dir, counter)
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

plotOverLine2.Source.Point1 = plotOverLine1.Source.Point2
plotOverLine2.Source.Point2 =[
        l_sym + l_wedge1*np.cos(theta_wedge1) + l_wedge2*np.cos(theta_wedge2),
        l_wedge1*np.sin(theta_wedge1) + l_wedge2*np.sin(theta_wedge2),
        0.5*w
]
plotOverLine1.Source.Resolution = args.resolution

# save data
temp_filename = '{}wedge2_data_{}.csv'.format(res_dir, counter)
SaveData(
    temp_filename,
    proxy=plotOverLine2,
    PointDataArrays=['Subdomain', 'T', 'alpha', 'arc_length', 'k', 'mu', 'p', 'qx', 'qy', 'qz', 'rho', 'rhoE', 'rhou', 'rhov', 'rhow', 'steady_state_error', 'txx', 'txy', 'txz', 'tyy', 'tyz', 'tzz', 'u', 'v', 'vtkValidPointMask', 'w']
)
print("Written file {}".format(temp_filename))

