#!/bin/bash

# script to convert raw GROMACS trajectory into analyze-able format
# usage: sh 1_gromacs_postproc.sh processed_traj_name tpr_file.tpr xtc_file.xtc

module load gromacs/2020.3

traj=$1
tpr=$2
xtc=$3
xtc_prefix=${xtc::-4}

# center traj, remove PBC artifacts, and condense into 1 ns resolution
gmx trjconv -s $tpr -f $xtc -o ${xtc_prefix}_center_nojump.xtc -center -pbc nojump -ur compact -dt 1000  << EOF
1
1
EOF

# cut off equilibration time (50 ns) for analysis
gmx trjconv -s $tpr -f ${xtc_prefix}_center_nojump.xtc -o ${traj}.xtc -center -pbc whole  -b 50000 << EOF
1
1
EOF

# convert to PDB for future use
gmx trjconv -s $tpr -f ${traj}.xtc -o ${traj}.pdb -fit rot+trans << EOF
1
1
EOF





