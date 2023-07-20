from plot_utils import TrajectoryPlotter

 # How to plot an individual trajectory (triplicate runs)

 # specify where the trajectories are and which one you want to use
loc = '/Users/henry/OneDrive - University of North Carolina at Chapel Hill/popov-lab/CM-allostery/trajectories/'

run = 'WT-Trp-T'
reps = ['t1/', 't2/', 't3/']

# initialize a TrajectoryPlotter object from plot_utils and use various functions to plot metrics
tp = TrajectoryPlotter(loc, run, reps)

tp.plot_local_metrics(metric='chi_selected', chi=True, savefig='CMTrp_T_chi.svg')
tp.plot_local_metrics(metric='dist_aux', chi=True, savefig='atmp.png')
tp.plot_rmsd()
tp.plot_rmsf()
tp.plot_subunit_metrics(add_xtal=True)

tp.plot_contact_map()
ref = 'data/1csm_contact.npy'
tp.plot_contact_map(ref=ref)

