import os
from plot_utils import MultiTrajectoryPlotter

# script for plotting multiple MD trajectories together

loc = '/Users/henrydieckhaus/OneDrive - University of North Carolina at Chapel Hill/popov-lab/CM-allostery/trajectories/'

run_list = [
    os.path.join(loc, 'WT-Tyr-T'),
    os.path.join(loc, 'T226I-Trp-R'),

    os.path.join(loc, 'T226I-Trp-R'),
    os.path.join(loc, 'WT-Tyr-T'),

    os.path.join(loc, 'WT-Trp-T'),
]

mtp = MultiTrajectoryPlotter(loc, run_list, add_xtal=True)

# can easily plot raw metrics on top of each other
mtp.plot_subunit_metrics(project=False)
mtp.plot_local_metrics(project=False)

mtp.plot_local_metrics(project=True, dist=False, use_trajs=[0, 1, 2, 3, 4, 5], metric='chi_selected', chi=True, savefig='CMTrp_metrics.svg')

# use_trajs indicates which trajs to use for PCA fitting - if not specified, all trajs will be used
# NOTE: if each traj is in triplicate, you must use ALL 3 indices to select it (0, 1, 2, etc.)
mtp.plot_local_metrics(project=True, dist=False, use_trajs=[0, 1, 2, 3, 4, 5], metric='chi_selected', chi=True, savefig='fig2A_top.svg')





