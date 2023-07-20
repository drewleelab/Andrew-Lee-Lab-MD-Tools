#!/bin/bash
#SBATCH --job-name=analysis
#SBATCH -p debug_queue
#SBATCH -N 1
#SBATCH --ntasks-per-node=44
#SBATCH --time=1:00:00 # time (D-HH:MM)

module load openmpi_4.0.1/gcc_11.2.0
module load gromacs/2020.3

# NOTE: need to submit this from the directory w/the traj data in it

traj="WT-Apo-T"
scripts="/proj/kpoplab/users/henry/CM/analysis_scripts"

echo "===================="
echo "Initial GROMACS Postprocessing"
# do initial gromacs postprocessing
sh "${scripts}/1_gromacs_postproc.sh" $traj

# make index file for custom atom-pair distance calculations
sh "${scripts}/2_custom_ndx.sh" "${traj}.xtc" md_0_1.tpr dist.ndx

echo "==================="
echo "Calculating GROMACS Metrics"
# calculate gromacs-based metrics: distances, chi angles, RMSF, RMSD
sh "${scripts}/3_gromacs_metrics.sh" "${traj}.xtc" md_0_1.tpr dist.ndx tmp.ndx "${scripts}/res_chi.txt" "${scripts}/distance_grps.txt" "${scripts}/res_rms.txt"

echo "==================="
echo "Dynamic Cross-Correlation"
# calculate dynamic cross-correlation
/nas/longleaf/home/dieckhau/scripts/MD-TASK/venv/bin/python \
	/nas/longleaf/home/dieckhau/scripts/MD-TASK/calc_correlation.py \
	--topology snap.pdb --trajectory "${traj}.xtc" --step 10 --prefix $traj --title $traj

# calculate helical/subunit descriptors
ref_Txtal="/proj/kpoplab/users/henry/CM/ref_structures/xtal/2csm_chainsAB.pdb"

echo "==================="
echo "Subunit Angles"
/nas/longleaf/home/dieckhau/miniconda3/envs/mdanalysis/bin/python \
	/nas/longleaf/home/dieckhau/scripts/trajectory_analysis/subunit_angles.py \
	-f "${traj}.pdb" -ref $ref_Txtal -o "${traj}_subunits.csv"

ref_Txtal="/proj/kpoplab/users/henry/CM/ref_structures/rebuilt/2csm-prepped-FINAL.pdb"
echo "==================="
echo "Contact Map"
# calculate contact map
/nas/longleaf/home/dieckhau/miniconda3/envs/mdanalysis/bin/python \
	/nas/longleaf/home/dieckhau/scripts/trajectory_analysis/contact_map.py \
	-f "${traj}.pdb" -t 10. -o "${traj}_contact_map.npy" -r $ref_Txtal

echo "================="
echo "Compiling Analysis Data"

mkdir analysis
mv dist chi chi_selected rmsd rmsf* *.xvg analysis
mv tmp.ndx dist.ndx snap.pdb analysis
mv *.png *.txt *.npy *.csv analysis


echo "=================="
echo "Analysis Complete"
