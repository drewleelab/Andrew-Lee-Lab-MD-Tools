#!/bin/bash
#SBATCH --job-name=1csm-I226T-
#SBATCH -p 528_queue
#SBATCH -N 4
#SBATCH --ntasks-per-node=44
#SBATCH --time=48:00:00 # time (D-HH:MM)

# make GMX input files
f="1aki.pdb"
prefix="${f::-4}"
mdp_dir="/proj/kpoplab/users/henry/CM/trajectories/mdp"

module load openmpi_4.0.1/gcc_11.2.0
module load gromacs/2020.3

grep -v "HOH" $f >> "${prefix}_clean.pdb"
gmx pdb2gmx -f "${prefix}_clean.pdb" -o "processed.gro" -water spce -ignh

# solvate
gmx editconf -f "processed.gro" -o "newbox.gro" -c -d 1.0 -bt cubic
gmx solvate -cp "newbox.gro" -cs spc216.gro -o "solv.gro" -p topol.top 

# add ions
gmx grompp -f "${mdp_dir}/ions.mdp" -c solv.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

# energy min
gmx grompp -f "${mdp_dir}/em.mdp" -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

# nvt equil
gmx grompp -f "${mdp_dir}/nvt.mdp" -o nvt.tpr -c em.gro -p topol.top -r em.gro -maxwarn 2
mpirun $MPI_HOSTS gmx_mpi mdrun -v -deffnm nvt

# npt equil
gmx grompp -f "${mdp_dir}/npt.mdp" -o npt.tpr -c nvt.gro -p topol.top -t nvt.cpt -r nvt.gro -maxwarn 2
mpirun $MPI_HOSTS gmx_mpi mdrun -v -deffnm npt

# full MD run
gmx grompp -f "${mdp_dir}/md.mdp" -o md_0_1.tpr -c npt.gro -p topol.top -t npt.cpt -r npt.gro -maxwarn 2
mpirun $MPI_HOSTS gmx_mpi mdrun -v -deffnm md_0_1

# continue run if walltime hit
# mpirun $MPI_HOSTS gmx_mpi mdrun -s md_0_1.tpr -cpi md_0_1.cpt -deffnm md_0_1
