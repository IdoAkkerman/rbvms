
**Lid driven cavity**

WORK IN PROGRESS

This test case involves a low Reynolds flow around a 2D cylinder
resulting in a laminar unsteady flow.

The Reynolds number is Re=100 based on the free-stream and diameter of the cylinder.
The test case is documented in:

Transient flow past a circular cylinder: A benchmark solution
Engelman, Michael S. & Jamnia, Mohammad-Ali
International Journal for Numerical Methods in Fluids 11(7): 985 -1000
https://doi.org/10.1002/fld.1650110706

***Howto run***

There are two meshes available a NURBS mesh in
```
von-karman-nurbs.mesh
```
this is a very coarse mesh. To get a meaning full simulation this mesh needs to be refined 
3 or 4 times at minimum.
The boundaries of the mesh are numbered as follows:
1. Inflow
2. Bottom and top sides
3. Cylinder
4. Outflow

The required functions defining the initial and boundary conditions, the forcing-term and
the diffusion are specified in `von-karman.c`. This file needs to be compiled, by running
```
gcc   -shared -o libvkvs.so -fPIC von-karman.c
```
in the `cases/von-karman-re100` directory.
Each time changes are made to the `.c` file this command needs to be run again.

A NURBS simulation can be started by running
```
bash run-nurbs
```
in the `cases/von-karman-re100` directory.
This script will compile `von-karman.c` for you.
The script assumes the `rbvms` executable islocated in `../../build/rbvms/rbvms`.
If this is incorrect the script can be modified.
It is not necessary to run `rbvms` via the script, the script is mainly a convenient way of communicating
what commands could/shoudl be run and what the command line options are.

The file run-nurbs documents the input parameters needed for the simulation.
The script assumes the simulation can run on 16 cores.
The simulation takes on the order of 1 hour to finish if 16 cores are indeed available.

The simulation generates 4 types of output.

Furthermore the forces 
1. A log file. This includes all solver information such as iteration counts and residuals.
2. An output file that includes time traces of forces on all boundaries and addition info on CFL-number and time-step.
3. A directory of solution files for visualization purposes.
4. A directory of restart files for restarting the simulation.

The results of the simulation can be visualized using Visit. The files in the solution directory can be used for this.

The forces in the output file can be postprocessed using a script:
```
python3 postpro.py
```
This will report unsteady  Lift and Drag coefficients, and Period and Strouhal numbers.
The routine monitor sign changes of the lift to determine the period.
Data is averaged over a number of periods. Validation data from the Engelman paper is also printed for a quick comparison.

The time trace of the Lift and Drag coefficients can also be plotted using a gnuplot routine:
```
gnuplot ClCd.gp
```
this routine should also work remotely if the ssh connection is established with a x-forwarding,
```
ssh -Y user@machine
```
The format of the output file is as follows:

Step Time Time-step CFL Fx@bnd1 Fy@bnd1 Fx@bnd2 Fy@bnd2 Fx@bnd3 Fy@bnd3 Fx@bnd4 Fy@bnd4


**Triangular mesh**
The simulation can also be performed using triangular elements.
But first a triangular mesh needs to be generated using the following command:

```
gmsh von-karman-tri.geo -2 format msh22 -clscale 0.2
```
with the last option `-clscale` the size of the elements can be controlled.
Doubling this number creates elements that are approximately 2 times larger,
resulting in a Degree-of-Freedom (DoF) count that is roughly half.
The distribution of the elements can be controlled by modifying the `geo` file.

To run the triangular mesh run:
```
bash run-tri
```

The run parameters are slightly different than the NURBS case. See `run-tet` for details.

Visualization and postprocessing is identical.
