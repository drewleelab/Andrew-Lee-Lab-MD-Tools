
f="1csm-I226T-Y212F.pdb"
prefix="${f::-4}"
mdp_dir="/proj/kpoplab/users/henry/CM/trajectories/mdp"

module load openmpi_4.0.1/gcc_11.2.0
module load gromacs/2020.3

grep -v "HOH" $f >> "${prefix}_clean.pdb"
cp -r "${mdp_dir}/charmm36-mar2019.ff" .
gmx pdb2gmx -f "${prefix}_clean.pdb" -o "processed.gro" -water spce -ignh -ff charmm36-mar2019

# solvate
gmx editconf -f "processed.gro" -o "newbox.gro" -c -d 1.0 -bt cubic
gmx solvate -cp "newbox.gro" -cs spc216.gro -o "solv.gro" -p topol.top 

# add ions
gmx grompp -f "${mdp_dir}/ions.mdp" -c solv.gro -p topol.top -o ions.tpr
echo "13 \n" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

