
#gmsh channel-flap-fsi.geo -2 -clscale 4.0 -format msh22
#gcc   -shared -o channel-flap-fsi.so -fPIC channel-flap-fsi.c
mesh=channel-flap-fsi.msh
exe=../../../build-precice2/rbvms/rbvms-fsi
##exe=../../../build/rbvms/rbvms-fsi
#exe=../../../dbg/rbvms/rbvms
#exe=../../../debug/rbvms/rbvms-fsi

rm output* log*

##valgrind --tool=memcheck \
mpirun -n 1 $exe -m $mesh -o 1 -r 0 \
--dyn-visc 0.01 --density 1.0 \
--strong-bdr "2 4 5" --outflow-bdr 3 --fsi-bdr 6 -l channel-flap-fsi.so \
-s 21 -dt 0.01 --dt_vis 0.01 -tf 2.5 \
-lt 1e-3 -ni 50  -ri 10  -pfs 100 | tee log0
