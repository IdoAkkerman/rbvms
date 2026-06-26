
#gcc -shared -o perpendicular-flap-fsi.so -fPIC perpendicular-flap-fsi.c

mesh=von-karman-tri.msh
exe=../../../build/rbvms/rbvms-fsi
##exe=../../../build/rbvms/rbvms-fsi
#exe=../../../dbg/rbvms/rbvms
#exe=../../../debug/rbvms/rbvms-fsi

rm output* log*

##valgrind --tool=memcheck \
mpirun -n 1 $exe -m $mesh -o 1 -r 0 \
--dyn-visc 1e-2  --density 1.0 \
--strong-bdr "1 2" --outflow-bdr 4 --fsi-bdr 3 -l libvkvs.so \
-s 45 -dt 0.1 --dt_vis 0.5 -tf 1000 \
-lt 1e-3 -ni 20  -ri 20  -pfs 1 | tee log0
