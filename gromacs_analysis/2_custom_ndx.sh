#!/bin/bash

# script to automatically generate index file for atom-pair distance calculations for CM protein specifically

# usage: sh custom_ndx.sh traj.xtc topol.tpr index.ndx

traj=$1
tpr=$2
ndx=$3
snap="snap.pdb"

module load gromacs/2020.3

# dump the first frame of trajectory to PDB file, which includes atom index info
gmx trjconv -f $traj -s $tpr -o $snap -dump 0 << EOF
1
EOF

# make custom index file out of PDB and custom list of index numbers (distance_index.ndx)
gmx make_ndx -f $snap -o $ndx  << EOF

keep 0

r 23 & a OE2 & chain A
r 23 & a OE2 & chain B

r 157 & a NH2 & chain A
r 157 & a NH2 & chain B

r 204 & a NH2 & chain A
r 204 & a NH2 & chain B

r 208 & a NZ & chain A
r 208 & a NZ & chain B

r 215 & a OD2 & chain A
r 215 & a OD2 & chain B

r 234 & a OH & chain A
r 234 & a OH & chain B

11 | 3
name 13 234_A_157_A
12 | 4
name 14 234_B_157_B

1 | 3
name 15 23_A_157_A
2 | 4
name 16 23_B_157_B

1 | 11
name 17 23_A_234_A
2 | 12
name 18 23_B_234_B

1 | 7
name 19 23_A_208_A
2 | 8
name 20 23_B_208_B

1 | 5
name 21 23_A_204_A
2 | 6
name 22 23_B_204_B

9 | 6
name 23 215_A_204_B
9 | 8
name 24 215_A_208_B

q
EOF
