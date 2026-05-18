
gmsh channel-flap-fsi.geo -2 -clscale 1.0 -format msh22

gcc   -shared -o channel-flap-fsi.so -fPIC channel-flap-fsi.c
mesh=channel-flap-fsi.msh
exe=../../../build/rbvms/rbvms-fsi
#exe=../../../dbg/rbvms/rbvms


rm output* log*

mpirun -n 1 $exe -m $mesh -o 1 -r 0 \
--strong-bdr "2 4 5 6" --outflow-bdr 3 -l channel-flap-fsi.so \
-s 45 -dt 0.01 --dt_vis 0.05 -tf 2.5 \
-lt 1e-3 -ni 20  -ri 10  | tee log0
