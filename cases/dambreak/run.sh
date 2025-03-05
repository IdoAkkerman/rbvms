mpicc -shared -o rbvms-ls.so -fPIC dambreak.c
#gmsh -2 dambreak.geo -format msh22 -clscale 0.025


mpirun -n 1 /home/ido/data/rbvms-ls2/deb/rbvms-ls/rbvms-ls -l rbvms-ls.so -m dambreak.msh --strong-bdr "1" -s 45 --dt_vis 0.01 \
--newton-tolerance 1e-8 --newton-itermax 30
 
