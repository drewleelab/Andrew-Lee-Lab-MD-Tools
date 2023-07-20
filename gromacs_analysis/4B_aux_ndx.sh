#!/bin/bash

# script to automatically generate index file for atom-pair distance calculations (auxiliary atoms)

# usage: sh custom_ndx.sh traj.xtc topol.tpr index.ndx

traj=$1
tpr=$2
ndx=$3
snap="snap.pdb"

module load gromacs/2020.3

gmx trjconv -f $traj -s $tpr -o $snap -dump 0 << EOF
1
EOF

gmx make_ndx -f $snap -o $ndx  << EOF

keep 0

r 16 & a NH2 & chain A
r 16 & a NH2 & chain B

r 23 & a OE2 & chain A
r 23 & a OE2 & chain B

r 24 & a OD2 & chain A
r 24 & a OD2 & chain B

r 157 & a NH2 & chain A
r 157 & a NH2 & chain B

r 168 & a NZ & chain A
r 168 & a NZ & chain B

r 208 & a NZ & chain A
r 208 & a NZ & chain B

r 212 & a OH & chain A
r 212 & a OH & chain B

r 246 & a OE2 & chain A
r 246 & a OE2 & chain B

1 | 3
name 17 23_A_16_A

2 | 4
name 18 23_B_16_B

5 | 13
name 19 24_A_212_A

6 | 14
name 20 24_B_212_B

5 | 11
name 21 24_A_208_A

6 | 12
name 22 24_B_208_B

7 | 9
name 23 157_A_168_A

8 | 10
name 24 157_B_168_B

9 | 15
name 25 168_A_246_A

10 | 16
name 26 168_B_246_B

q
EOF
