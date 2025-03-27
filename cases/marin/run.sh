#!/bin/bash

# Set variables
exe=../../build/rbvms-ls/rbvms-ls
mesh=marin.msh

# Clean
#rm -rf solution/* output_*.dat

# Create functions
mpicc -shared -o rbvms-ls.so -fPIC marin.c

# Create mesh -- if required
if [ ! -f $mesh ]; then
   gmsh -3 marin.geo -format msh22 -clscale 0.075
fi

# Actual run
nohup mpirun -n 24 --oversubscribe $exe \
-l rbvms-ls.so -m $mesh --normal-bdr "1" \
-s 45 --dt_vis 0.01 --dt 0.01 \
--newton-tolerance 1e-3 --newton-itermax 10 \
--linear-tolerance 1e-3 --linear-itermax 150 -ri 5 > log &


tail -1111f log
