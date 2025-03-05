mpicc -shared -o rbvms-ls.so -fPIC dambreak.c
#gmsh -2 dambreak.geo -format msh22 -clscale 0.025

rm solution/*

mpirun -n 16 /home/ido/data/rbvms-ls2/deb/rbvms-ls/rbvms-ls -l rbvms-ls.so -m dambreak.msh --weak-bdr "1" -s 45 -dt 0.0005 --dt_vis 0.01 \
--newton-tolerance 1e-4 --newton-itermax 20 | tee log
 
