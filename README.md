# CM-Allostery MD Analysis
This repository includes all of the scripts used to run, analyze,
and visualize molecular dynamics trajectories of yeast chorismate mutase (CM)
with different mutations and ligand parameters.

## Installation
To clone this repository, use ```git clone {link}``` or download manually.
All necessary dependencies are provided in the ```environment.yaml``` file.
It is recommended to use python package manager ```mamba``` for most convenient installation.

## Contents

Input PDBs used as starting points for MD runs are provided in ```pdbs/```.
Example runscripts for GROMACS are provided in ```gromacs_run/```.
Initial analysis using a combination of GROMACS and shell scripting are provided in ```gromacs_analysis/```.
Custom python scripts for further analysis with pymol etc. are provided in the base directory.
Visualization utilities are provided in ```visualization/```.