mpicc -shared -o rbvms-ls.so -fPIC cahouet.c
#gmsh -2 cahouet.geo -format msh22 -clscale 0.25

#rm -rf solution/*

##mpirun -n 16 


mpirun -n 16 --oversubscribe \
/home/ido/data/rbvms-ls2/build/rbvms-ls/rbvms-ls -l rbvms-ls.so -m cahouet.msh \
--weak-bdr "1 3 4" --suction-bdr "2" -s 21 -dt 0.06 --dt_vis 0.5  -tf 1000.0 \
--linear-tolerance 1e-3 --linear-itermax 250 \
--newton-tolerance 1e-3 --newton-itermax 5 --redist-penalty 2.0 -ri 10 -rs | tee log 

#mpirun -n 2 valgrind --tool=memcheck --track-origins=yes \
#/home/ido/data/rbvms-ls2/deb/rbvms-ls/rbvms-ls -l rbvms-ls.so -m dambreak.msh --weak-bdr "1" -s 32 -dt 0.01 --dt_vis 0.01 \
#--linear-tolerance 1e-3 --linear-itermax 150 \
#--newton-tolerance 1e-5 --newton-itermax 2 -dist 3 | tee log 
 
