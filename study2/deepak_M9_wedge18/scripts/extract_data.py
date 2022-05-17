# trace generated using paraview version 5.9.0

#### import the simple module from the paraview
from paraview.simple import *
import argparse

parser = argparse.ArgumentParser(
    description = "A script to extract line data from pvtu files for Deepak et al (2013) "
        "Mach 9 18 degree wedge case. " +
        "This script has to be executed using 'pvpython'."
)
parser.add_argument("res_dir", help="Result directory (absolute or relative).")
parser.add_argument("counter", help="Output counter.", type=int)
args = parser.parse_args()

res_dir = args.res_dir
if res_dir[-1] != "/":
    res_dir += "/"
counter = args.counter

temp_filename = "{}output_{:06d}.pvtu".format(res_dir, counter)
print("Reading file {}".format(temp_filename))
result = XMLPartitionedUnstructuredGridReader(
    FileName=temp_filename
)

# create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(
    registrationName='PlotOverLine1',
    Input=result,
    Source='Line'
)

z = 1e-3
plotOverLine1.Source.Point1 = [0.01, 0.0, z]
plotOverLine1.Source.Point2 = [0.085, 0.0, z]

# save data
temp_filename = '{}plate_data_{}.csv'.format(res_dir, counter)
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
plotOverLine2.Source.Point2 = [0.17535, 0.0293566, z]

# save data
temp_filename = '{}ramp_data_{}.csv'.format(res_dir, counter)
SaveData(
    temp_filename,
    proxy=plotOverLine2,
    PointDataArrays=['Subdomain', 'T', 'alpha', 'arc_length', 'k', 'mu', 'p', 'qx', 'qy', 'qz', 'rho', 'rhoE', 'rhou', 'rhov', 'rhow', 'steady_state_error', 'txx', 'txy', 'txz', 'tyy', 'tyz', 'tzz', 'u', 'v', 'vtkValidPointMask', 'w']
)
print("Written file {}".format(temp_filename))

