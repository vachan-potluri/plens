
from paraview.simple import *

# find source
dummypvd = XMLPartitionedUnstructuredGridReader(FileName='dummy.pvd')

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=dummypvd)
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
SaveData(
    'surface_data.csv',
    proxy=plotOnIntersectionCurves1,
    PointDataArrays=['Subdomain', 'T', 'alpha', 'arc_length', 'k', 'loc_dt', 'mu', 'p', 'qx', 'qy', 'qz', 'rho', 'rhoE', 'rhou', 'rhov', 'rhow', 'steady_state_error', 'txx', 'txy', 'txz', 'tyy', 'tyz', 'tzz', 'u', 'v', 'w'],
    UseScientificNotation=1
)
