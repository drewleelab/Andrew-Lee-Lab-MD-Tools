#!/bin/bash

# script to extract chi angles, distances, and rms values for residues of interest and compile into nice format
# usage: sh all_analysis.sh traj.xtc md_topol.tpr distance_index.ndx new_index.ndx res_chi.txt distance_grps.txt res_rms.txt

traj=$1
top=$2
dist_index=$3
index=$4
chi_res=$5
dist=$6
rms_res=$7

# need to load latest version or chi command will fail
module load gromacs/2020.3

# calculate ALL chi angles
gmx chi -s $top -f $traj -maxchi 2 -all

mkdir chi
mkdir chi_selected
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
	cp "chi/${p}.xvg" "chi_selected/${new_name}"
done < $chi_res

mkdir dist

while read grp; do
        echo "$grp"
        gmx distance -s $top -f $traj -n $dist_index -oall -oallstat -select "$grp"
        mv "dist.xvg" "dist/${grp}_dist.xvg"
        mv "diststat.xvg" "dist/${grp}_diststat.xvg"
done < $dist

# calculate rmsf per-residue
echo -e "1 \n" | gmx rmsf -s $top -f $traj -o "rmsf-per-res.xvg" -res

# initialize index file if not existing, keep only protein group
echo -e "keep 1 \n q \n" | gmx make_ndx -f $top -o $index
mkdir rmsd

while read p; do
        echo $p
        # make new group with single res in it
        echo -e "ri ${p} \n q \n" | gmx make_ndx -f $top -n $index -o $index

        if [[ $p -gt 254 ]]
        then
                echo "Chain B"
                resnum=$(($p-254))
                new_name="${resnum}_chB_rmsd.xvg"
        else
                echo "Chain A"
                new_name="${resnum}_chA_rmsd.xvg"
        fi

        # calculate rms for single res
        echo -e "0 \n 1 \n" | gmx rms -s $top -f $traj -n $index -o "rmsd/$new_name"
        # remove single res group
        echo -e "del 1 \n q \n" | gmx make_ndx -f $top -n $index -o $index

done < $rms_res

# calculate RMSD relative to 1csm/2csm xtal structures on a per-subunit basis

ref="/proj/kpoplab/users/henry/CM/ref_structures/rebuilt"

# using C-alphas only for backbone fit
echo -e "3 \n 3 \n" | gmx rms -s md_0_1.tpr -f $traj -o rmsd-vs-self.xvg -fit rot+trans
echo -e "14 \n 14 \n" | gmx rms -s "${ref}/1csm-prepped-FINAL.pdb" -f $traj -fit rot+trans -o "rmsd-vs-1csm.xvg" -n "${ref}/index_1csm.ndx"
echo -e "14 \n 14 \n" | gmx rms -s "${ref}/2csm-prepped-FINAL.pdb" -f $traj -fit rot+trans -o "rmsd-vs-2csm.xvg" -n "${ref}/index_2csm.ndx"


# clean up temporary files
rm "#"*".ndx."*

