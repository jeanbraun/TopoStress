# TopoStress Documentation
This is a very brief documentation on **TopoStress** a 3D finite element code to compute the elastic stresses within a crustal block bounded along its top surface by an arbitrary topography

The code was developed by Jean Braun, GFZ-Potsdam

It uses simple bi-linear interpolation functions defined in 8-node elements that have a square base but an arbitrary set of vertical heights, thus allowing to map an arbitrary surface topography known at the nodes of a rectangular mesh and “morph” it into a flat bottom at a prescribed depth

The material is assumed to be elastic with known Young’s modulus and Poisson’s ratio. A simple von Moses criterion is also used to compute elasto-plastic stresses if required ( a negative yield criterion implies linear elastic behavior)

Several boundary conditions can be imposed including:
- free boundaries
- free-slip boundaries
- no-slip boundaries
- no-tilt boundaries
In all cases the vertical displacement os forced to be nil on the bottom layer of nodes

To obtain the solution , a large positive definite system of algebraic equations must be solved. For this a Cholesky factorization and back substitution method is used that was developed by Jean Braun but uses the node numbering method/code developed by Sloan (1992).

## Interface
**TopoStress** is written in Fortran but an interface is provided that allows the code to be called from many other languages, including Python. The interface is also used to link the code to the fastscape-lem package that solves the Stream Power Law using the method developed by Braun and Willett (2013) and implemented using the Xarray-simlab modeling framework developed by Benoit Bovy while at the GFZ-Potsdam.

Basic instructions for using the interface are given in the examples section.

## Installation
**TopoStress** must be compiled using a fortran compiler. We strongly suggest to use the Fortran compiler. This is all that is required to use the code from another fortran code.

To produce a python callable library, **TopoStress** must also be compiled using `f2py`.

We suggest to use the `environment.yml` file to create a `conda` environment  (named `TopoStress`) that contains all required libraries:
`conda env create -f environment.yml `

To perform both operations, the user should go into the `src` directory:
`cd src`
and compile the code using the Makefile:
`make`
This will create a fortran library and a python library in the `lib` directory. A copy of the python library will also be created in the `examples` directory.

## Running examples
All examples are located in the examples directory.

### Using Fortran code
To compile and execute the fortran example, the user should go into the examples directory:
`cd examples`
and compile the fortran example:
`make`
To run the fortran example, the use should run the executable:
`./TopoStresses_run`
This example uses a test topography stored in the file `Topo_Test.txt`. This test topography is given on a 40x40 node grid. The user can supply any other topography at any other resolution but will need to edit the `TopoStresses_run.f90` code to adapt it to the new topography. Other parameters (like the boundary conditions) can be modified too.
 
### Using Python code
The python example is also in the examples directory and is contained in a python file (`TopoStresses_run.py`) and a jupyter notebook (`TopoStresses_run.ipynb`). They both need the python library file `TopoStresses.so` to be in the same directory as the python/jupyter notebook files are.

Following are basic indications on how to use the interface if called from python, assuming that the library has been imported using:
`Import TopoStresses as TS`

First the user must initialize the interface by specifying the resolution of the input topography `(nx,ny)` and the resolution of the model in the z-direction, `nz`:
`TS.topostresses_setup(nx, ny, nz)`

The size of the model run in the x-, y-, and z-directions (in meters)
`TS.topostresses_set_xlylzl(xl, yl, zl)`

Mechanical properties, including `rhog`, the product of rock density (in kg/m^3) by gravitational attraction (in m/s^2), `exx` and `eyy` two components of horizontal strain imposed on the crustal block (no units), `ym`, Young’s modulus (in Pascals), `pr`, Poisson’s ratio (no unit) and `ks` a yield stress (in Pascals, according to von Mises criterion)
`TS.topostresses_set_mechanical_properties(rhog, exx, eyy, ym, pr, ks)`

Boundary conditions that correspond to the following bounty conditions for the side boundaries:
- `ibc=0` => free boundaries
- `ibc=1` => no-slip boundaries
- `ibc=2` => free-slip boundaries
- `ibc=3` => no-tilt bonudaries
`TS.topostresses_set_bc(ibc)`

A 1D array representing the 2D topography (in m)
`TS.topostresses_set_h(h.flatten())`

To run the code, the user needs to execute it using:
`TS.topostresses_execute()`

To get the output (the stress tensor at every point of the grid), the user must copy it through the interface into a double precision array that has the proper dimension (6 times the number of nodes) and for which memory has been allocated: 
`stress = np.zeros((6*nx*ny*nz))`
`TS.topostresses_copy_stress(stress)`
The array contains in its six components:
1. sigma_xx
2. sigma_yy
3. sigma_zz
4. tau_xy
5. tau_zx
6. tau_yz

It is important that the user frees the memory allocated internally by the interface by calling the destroy routine:
`TS.topostresses_destroy()`

### Coupling to FastScape using the Xarray-simlab framework
In the jupyter notebook `TopoStressFramework.ipynb`, a topography created by fastcscape-lem is used to compute the stresses.

In this implementation of **TopoStress**, the user specifies the model parameters as an xarray DataSet (see Xarray-simlab documentation). 

