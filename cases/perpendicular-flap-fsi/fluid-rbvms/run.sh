#gmsh channel-flap-fsi.geo -2 -clscale 1.0 -format msh22
gcc -shared -o perpendicular-flap-fsi.so -fPIC perpendicular-flap-fsi.c
mesh=mesh-flap.msh
exe=../../../build/rbvms/rbvms

rm output* log*

mpirun -n 1 $exe -m $mesh -o 1 -r 0 \
--dyn-visc 1.0 --density 1.0 \
--strong-bdr "1 20 3 4 5" --outflow-bdr 21 -l perpendicular-flap-fsi.so \
-s 45 -dt 0.01 --dt_vis 0.1 -tf 5.0 \
-lt 1e-3 -ni 20  -ri 10  | tee log0
