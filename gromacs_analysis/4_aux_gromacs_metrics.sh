#!/bin/bash

# script to extract additional metrics from GROMACS postprocessed trajectories
# usage: sh all_analysis.sh traj.xtc md_topol.tpr distance_index.ndx new_index.ndx res_chi.txt distance_grps.txt

traj=$1
top=$2
dist_index=$3
index=$4
chi_res=$5
dist=$6

# need to load latest version or chi command will fail
module load gromacs/2020.3

# calculate ALL chi angles
gmx chi -s $top -f $traj -maxchi 2 -all

mkdir chi
mkdir chi_aux
mv *.xvg chi.log chi

# pull out specified chi angles
while read p; do
	echo "$p"
	# renumber residues to match b/w chains
	resnum="${p:7}"
	if [[ $resnum -gt 254 ]]
	then
		echo "Chain B"
		resnum=$(($resnum-254))
		new_name="${p::7}${resnum}_chB.xvg"
	else
		echo "Chain A"
		new_name="${p::7}${resnum}_chA.xvg"
	fi
	cp "chi/${p}.xvg" "chi_aux/${new_name}"
done < $chi_res

mkdir dist_aux

while read grp; do
        echo "$grp"
        gmx distance -s $top -f $traj -n $dist_index -oall -oallstat -select "$grp"
        mv "dist.xvg" "dist_aux/${grp}_dist.xvg"
        mv "diststat.xvg" "dist_aux/${grp}_diststat.xvg"
done < $dist

# clean up temporary files
rm "#"*".ndx."*
rm -rf chi
