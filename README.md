
# RBVMS

Application that can simulate incompressible flow
using the residual-based variational multiscale (RBVMS) formulation.
This novel formulation is also suitable as a Large-Eddy Simulation model at turbulent Reynolds numbers.

Novel features include:
* Formulation combining numerical and physical arguments
* * Local Momentum conservation
* * Local Angular momentum conservation
* Wide range of discretisation possible, including
* * Tetrahedral
* * Isogeometric
* Weak imposition of Dirichlet boundary conditions
* Conservative boundary forces extraction
* Symmetry BC -- TBD

For details on the formulation see:

Y. Bazilevs, I. Akkerman,
Large eddy simulation of turbulent Taylor–Couette flow using isogeometric analysis and the residual-based variational multiscale method,
Journal of Computational Physics,
Volume 229, Issue 9, 2010, Pages 3402-3414,ISSN 0021-9991,
https://doi.org/10.1016/j.jcp.2010.01.008.

The parallel application is build on [mfem](https://github.com/mfem/mfem), a
lightweight, general, scalable C++ library for finite element methods.

## Install -- Local

To install RBVMS on a Ubuntu-based machine make sure the appropriate packages are installed.

Before installing the packages make sure the machine is updated, by running the folliwing commands in a terminal:

```
sudo apt-get update
sudo apt-get upgrade
```

Install the required packages by running the following command in a terminal:

```
sudo apt-get install cmake git build-essential openmpi-common libopenmpi-dev libmetis5 libmetis-dev
```

The RBVMS code can be downloaded using the following command:

```
git clone git@github.com:IdoAkkerman/rbvms.git

```
or
```
git clone https://github.com/IdoAkkerman/rbvms.git

```
This will create a directory `rbvms` with the code in it. This code can be compiled using the following commands:

```
cd rbvms
mkdir build
cd build
cmake ..
make -j 16
```

Here the number 16 is the number of files that will be compiled concurrently.
NB: This number you choose should depend on the number of cores you have available.

Depending on the machine compilation might take a few minutes.
After this hypre and MFEM will be installed in the `build/external` directory, and
the rbvms executable can be found in the `bin` directory.

## Install -- Clusters

### Delft Blue

To run on Delf Blue please see,
[DelftBlue.md](DelftBlue.md)

### Reynolds

TBD:
To run on Reynolds please see,
[Reynolds.md](Reynolds.md)

## Test Cases


2D test cases:
1. Poiseuille flow: go to directory `cases/poiseuille`
2. Lid driven cavity: go to directory `cases/lid-driven-cavity`
3. Von Karman vortex street: go to directory `cases/von-karman-re100`
4. Naca foil: go to directory `cases/naca-foil`

3D test cases: tbd
1. Channel flow
2. Von Karman vortex street
3. Naca foil
4. Gresho vortex
5. Taylor-green vortex
6. Turbulent channel
7. Taylor Couette flow

## Run
The code can be run by running
```
$BUILD/bin/rbvms
```

Where `$BUILD` is the directory where `cmake` was run.

To see all the command line options run:
```
$BUILD/bin/rbvms -h
```

### Isogeometric


### GMSH

GMSH is an open-source mesher. On Ubuntu-based machine the package can be installed
using the following command:

```
sudo apt-get install gmsh
```

### Visualization using Visit

The results can be visualized using [Visit] https://visit-dav.github.io/visit-website/index.html.
Which is an open-source tool that can natively visualize NURBS based solutions.

# Licenses

**RBVM** is distributed under the terms of the Apache License (Version 2.0).

All new contributions must be made under the Apache-2.0 licenses.

**MFEM** is distributed under the terms of the BSD-3 license. All new contributions 
must be made under this license. See LICENSE and NOTICE for details.

SPDX-License-Identifier: BSD-3-Clause
LLNL Release Number: LLNL-CODE-806117
DOI: 10.11578/dc.20171025.1248

**HYPRE** is distributed under the terms of both the MIT license and the Apache License (Version 2.0). 
Users may choose either license, at their option.

All new contributions must be made under both the MIT and Apache-2.0 licenses.

See LICENSE-MIT, LICENSE-APACHE, COPYRIGHT, and NOTICE for details.

SPDX-License-Identifier: (Apache-2.0 OR MIT)

LLNL-CODE-778117

