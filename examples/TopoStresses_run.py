import TopoStresses as TS
import numpy as np
import matplotlib.pyplot as plt

h = np.loadtxt('Topo_test.txt') # reads the test topography

nx = 40 # model resolution in the x-direction
ny = 40 # model resolution in the y-direction
nz = 5 # model resolution in the z-direction
xl = 10e3 # size of the model in the x-directoin (in m)
yl = 10e3 # size of the model in the y-directoin (in m)
zl = 10e3 # size of the model in the z-directoin (in m)
rhog = 2800*9.81 # product of rho (density in kg/m^3) by g (in m/s^2)
ym = 1e11 # Young's modulus (in Pa)
pr = 0.25 # Poisson's ratio
ks = -1. # yield stress (in Pa)
exx = 0. # imposed strain in the x-directio
eyy = 0. # imposed strain in the y-direction
ibc = 2 # side bondary conditions
# ibc=0 is free displacements 
# ibc=1 is no-slip
# ibc=2 is free-slip
# ibc=3 is no-tilt

TS.topostresses_setup(nx, ny, nz) # initializes and sends the model resolution through the interface
TS.topostresses_set_xlylzl(xl, yl, zl) # sends the model sizes through the interface
TS.topostresses_set_mechanical_properties(rhog, exx, eyy, ym, pr, ks) # sends the mechanical properties through the interface
TS.topostresses_set_bc(ibc) # sets the boundary conditions through the interface
TS.topostresses_set_h(h.flatten()) #  sends the surface topography through the interface
TS.topostresses_execute() # executes TopoStress
stress = np.zeros((6*nx*ny*nz)) # allocates memory to store the stress tensor in a 1D array
TS.topostresses_copy_stress(stress) # recovers the computed stresses through the interface
TS.topostresses_destroy() # releases memory through the interface

print(stress.reshape((6,nz,ny,nx))[3,-1,:,:].min(),stress.reshape((6,nz,ny,nx))[3,-1,:,:].max())