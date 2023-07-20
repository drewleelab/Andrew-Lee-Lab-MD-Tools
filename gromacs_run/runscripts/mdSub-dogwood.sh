#!/bin/bash
#SBATCH --job-name=1csm2L
#SBATCH -p 528_queue
#SBATCH -N 4
#SBATCH --ntasks-per-node=44
#SBATCH --time=48:00:00 # time (D-HH:MM)
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
gmx grompp -f nvt.mdp -o nvt.tpr -c em.gro -p topol.top -r em.gro -maxwarn 2
mpirun $MPI_HOSTS gmx_mpi mdrun -v -deffnm nvt
gmx grompp -f npt.mdp -o npt.tpr -c nvt.gro -p topol.top -t nvt.cpt -r nvt.gro -maxwarn 2
mpirun $MPI_HOSTS gmx_mpi mdrun -v -deffnm npt
gmx grompp -f md.mdp -o md_0_1.tpr -c npt.gro -p topol.top -t npt.cpt -r npt.gro -maxwarn 2
mpirun $MPI_HOSTS gmx_mpi mdrun -v -deffnm md_0_1


mpirun $MPI_HOSTS gmx_mpi mdrun -s md_0_1.tpr -cpi md_0_1.cpt -deffnm md_0_1
