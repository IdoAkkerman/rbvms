#!/bin/sh
#
#SBATCH --job-name="rbvms-von-karman"
#SBATCH --time=01:00:00
#SBATCH --ntasks=16
#SBATCH --partition=compute
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G

DIR=$HOME/rbvms
EXE=$DIR/build/bin/rbvms

source $DIR/config/delft-blue-modules.sh

srun $EXE \
--mesh von-karman-nurbs.mesh --order 2 --refine 4 \
--lib libvkvs.so --strong-bdr "1 2" --outflow-bdr 4  --weak-bdr 3 \
--ode-solver 45 --dt 0.1 --t-final 1000 \
--dt_vis 0.5 --restart-interval 10 \
--linear-tolerance 1e-3 --newton-itermax 20
