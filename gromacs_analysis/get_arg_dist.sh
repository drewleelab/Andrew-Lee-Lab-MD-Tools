#!/bin/bash

# script to automatically generate arg 16 - glu 23 distance

# usage: sh get_arg_dist.sh traj.pdb topol.tpr index.ndx

traj=$1
tpr=$2
ndx=$3
snap="snap.pdb"

module load gromacs/2020.3

# need PDB file to get chain IDs for atom selection (XTC doesn't have them)
gmx trjconv -f $traj -s $tpr -o $snap -dump 0 << EOF
1
EOF

# make index groups for dist calculation
gmx make_ndx -f $snap -o $ndx  << EOF

keep 0

r 23 & a OE2 & chain A
r 23 & a OE2 & chain B

r 16 & a NH2 & chain A
r 16 & a NH2 & chain B

1 | 3
name 5 23_A_16_A
2 | 4
name 6 23_B_16_B

q
EOF

# get distance values
gmx distance -s $tpr -f $traj -n $ndx -oall "dist_23_A_16_A.xvg" -select 23_A_16_A
gmx distance -s $tpr -f $traj -n $ndx -oall "dist_23_B_16_B.xvg" -select 23_B_16_B

