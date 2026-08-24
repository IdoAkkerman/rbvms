
gcc -shared -o perpendicular-flap-fsi.so -fPIC perpendicular-flap-fsi.c

mesh=mesh-flap-nurbs.mesh
exe=../../../build/rbvms/rbvms-fsi
##exe=../../../build/rbvms/rbvms-fsi
#exe=../../../dbg/rbvms/rbvms
#exe=../../../debug/rbvms/rbvms-fsi

rm output* log*

#valgrind --tool=memcheck \
mpirun -n 22 $exe -m $mesh -o 2 -r 2 \
--dyn-visc 1.0 --density 1.0 \
--strong-bdr "1 4 5" --outflow-bdr "20 21" --fsi-bdr 3 -l perpendicular-flap-fsi.so \
-s 45 -dt 0.001 --dt_vis 0.01 -tf 5.0 \
-lt 1e-5 -nt 1e-5 -ni 30  -ri 20  -pfs 1 | tee log0
