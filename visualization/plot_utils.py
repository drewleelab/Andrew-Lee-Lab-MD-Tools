import matplotlib.pyplot as plt
import os
import seaborn as sns
import matplotlib.gridspec as gridspec
import sys
sys.path.append('../')
from data_utils import *
from copy import deepcopy

import umap
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import colorcet as cc
from scipy.spatial.distance import cdist, pdist, squareform


class TrajectoryPlotter:
    """Shared Plotter class used for plotting a SINGLE MD Trajectory analysis data"""
    def __init__(self, loc, run, reps):

        self.loc = loc
        self.run = run
        self.reps = reps

    def plot_local_metrics(self, metric, chi, savefig=''):
        """Plot local (rotamer) metrics for a Trajectory"""
        traj_list, label_list = [], []
        for traj in self.reps:
            df = xvg_to_df(os.path.join(self.loc, self.run, traj + metric), chi=chi)
            traj_list.append(df)
            label_list.append(self.run + '(%s)' % (traj.rstrip('/')))

        total_df = pd.concat(traj_list)
        plot_grid(traj_list[0].columns, traj_list, label_list, save=savefig)

    def plot_rmsd(self, savefig=''):
        """Plot RMSDs from MD run in triplicate"""
        fig, ax = plt.subplots(1, 3, figsize=(12, 4))
        for n, rmsd in enumerate(['vs-self', 'vs-1csm', 'vs-2csm']):
            for traj in self.reps:
                arr = read_xvg(os.path.join(self.loc, self.run, traj + 'rmsd-' + rmsd + '.xvg'))[0]
                ax[n].scatter(arr[:, 0], arr[:, 1] * 10)
                ax[n].set_xlabel('Time (ps)')
                ax[n].set_ylabel('RMSD (A)')
                ax[n].set_title(rmsd)
                ax[n].set_ylim((0, 7.5))

        plt.suptitle('RMSD')
        plt.tight_layout()
        if savefig != '':
            plt.savefig(savefig)

        plt.show()

    def plot_rmsf(self, savefig=''):
        """Plot RMSF for whole protein + key resgions"""
        traj_list, label_list = [], []
        for traj in self.reps:
            df = read_xvg(os.path.join(self.loc, self.run, traj + 'rmsf-per-res.xvg'))[0]
            traj_list.append(df)
            label_list.append(self.run + '(%s)' % (traj.rstrip('/')))

        for traj, label in zip(traj_list, label_list):
            plt.bar(traj[:, 0], traj[:, 1], label=label, alpha=0.5)

        plt.ylim((0, 0.8))
        plt.legend()
        if savefig != '':
            plt.savefig(savefig)

        plt.show()

        # show EBR region
        for traj, label in zip(traj_list, label_list):
            plt.bar(traj[43:107, 0], traj[43:107, 1], label=label, alpha=0.5)
        plt.ylim((0, 0.8))
        plt.legend()

        if savefig != '':
            plt.savefig('EBR.png')
        plt.show()

        # show 11-12 loop region
        for traj, label in zip(traj_list, label_list):
            plt.bar(traj[212:226, 0], traj[212:226, 1], label=label, alpha=0.5)
        plt.ylim((0, 0.8))
        plt.legend()

        if savefig != '':
            plt.savefig('11-12-loop.png')
        plt.show()

    def plot_subunit_metrics(self, add_xtal=False, savefig=''):
        """Plotting subunit metrics (helical angles, etc.)"""
        traj_list, label_list = [], []
        run_name = self.run.split('/')[-1]
        for traj in self.reps:
            df = pd.read_csv(os.path.join(self.loc, self.run,  traj + run_name + '_subunits.csv'))
            traj_list.append(df)
            label_list.append(run_name + ' (%s)' % (traj.rstrip('/')))

        savenames = ['', '', '', '']
        if savefig != '':
            savenames = [savefig + 'subunit_metrics.png', savefig + 'helix_angles_1.png',
                         savefig + 'helix_angles_2.png', savefig + 'helix_angles_3.png']

        plot_grid(traj_list[0].columns[0:2], traj_list, label_list, save=savenames[0], xtal=True)
        plot_grid(traj_list[0].columns[0:15], traj_list, label_list, save=savenames[1], xtal=True)
        plot_grid(traj_list[0].columns[15:30], traj_list, label_list, save=savenames[2], xtal=True)
        plot_grid(traj_list[0].columns[30:], traj_list, label_list, save=savenames[3], xtal=True)

    def plot_contact_map(self, savefig='', ref=None):
        """Plot avg contact frequency map"""
        traj_list = []
        run_name = self.run.split('/')[-1]
        for traj in self.reps:
            arr = np.load(os.path.join(self.loc, self.run, traj + run_name + '_contact_map.npy'))
            traj_list.append(arr)

        avg_map = np.abs(np.mean(traj_list, axis=0))
        plt.imshow(avg_map, cmap='RdBu', origin='lower', vmin=-1, vmax=1)
        plt.colorbar()
        plt.title('Contact Frequency Map')

        if savefig != '':
            plt.savefig(savefig)

        plt.show()

        # if provided, use ref matrix to calculate difference map
        if ref is not None:
            ref = np.load(ref)
            diff = avg_map - ref
            plt.imshow(diff, cmap='RdBu', origin='lower', vmin=-1, vmax=1)
            plt.colorbar()
            plt.title('Contact Difference Map (vs ref)')
            if savefig != '':
                plt.savefig('diff_' + savefig)

            plt.show()


class MultiTrajectoryPlotter:
    """Class for plotting MULTIPLE trajectories in triplicate"""
    def __init__(self, loc, run_list, add_xtal=False):
        self.loc = loc
        self.run_list = run_list
        self.add_xtal = add_xtal

    def plot_local_metrics(self, project=False, savefig='', dist=False, method='pca', use_trajs=None,
                           metric='chi_selected', chi=True):
        """Plot or project local metrics with PCA/UMAP methods"""
        traj_list = []
        label_list = []
        for run in self.run_list:
            run_name = run.split('/')[-1]
            for traj in ['t1/', 't2/', 't3/']:
                df = xvg_to_df(os.path.join(self.loc, run, traj + metric), chi=chi)
                traj_list.append(df)
                label_list.append(run_name + ' (%s)' % (traj.rstrip('/')))

        if self.add_xtal:
            df_1csm, df_2csm = load_xtal_as_df(impute=True, subunit=False)
            traj_list = traj_list + [df_1csm, df_2csm]
            label_list = label_list + ['1csm xtal (R)', '2csm xtal (T)']

        if project:
            project_data(traj_list, label_list, method=method, use_trajs=use_trajs,
                         save=savefig, legend=False, dist=dist)
        else:
            plot_grid(traj_list[0].columns, traj_list, label_list, save=savefig)

    def plot_subunit_metrics(self, project=False, savefig='', dist=False, method='pca', use_trajs=None):
        """Plot or project subunit metrics with PCA/UMAP methods"""
        traj_list = []
        label_list = []
        for run in self.run_list:
            run_name = run.split('/')[-1]
            for traj in ['t1/', 't2/', 't3/']:
                df = pd.read_csv(os.path.join(self.loc, run,  traj + run_name + '_subunits.csv'))
                traj_list.append(df)
                label_list.append(run_name + ' (%s)' % (traj.rstrip('/')))

        if self.add_xtal:
            df_1csm, df_2csm = load_xtal_as_df(impute=True, subunit=True)
            traj_list = traj_list + [df_1csm, df_2csm]
            label_list = label_list + ['1csm xtal (R)', '2csm xtal (T)']

        if project:
            project_data(traj_list, label_list, method=method, use_trajs=use_trajs,
                         save=savefig, legend=False, dist=dist)
        else:
            savenames = ['', '', '', '']
            if savefig != '':
                savenames = [savefig + 'subunit_metrics.png', savefig + 'helix_angles_1.png',
                             savefig + 'helix_angles_2.png', savefig + 'helix_angles_3.png']
            plot_grid(traj_list[0].columns[0:2], traj_list[:-2], label_list[:-2], save=savenames[0], xtal=True)
            plot_grid(traj_list[0].columns[0:15], traj_list[:-2], label_list[:-2], save=savenames[1], xtal=True)
            plot_grid(traj_list[0].columns[15:30], traj_list[:-2], label_list[:-2], save=savenames[2], xtal=True)
            plot_grid(traj_list[0].columns[30:], traj_list[:-2], label_list[:-2], save=savenames[3], xtal=True)


def grouped_cluster_distances(df, data, savefig=''):
    """Calculates distances b/w groups of 3 MD traj run clusters"""
    trajs = df['TRAJ'].unique()
    centroids = []
    spreads = []
    for n in range(0, len(trajs), 3):
        points = np.concatenate([data[df['TRAJ'] == n, :],
                                 data[df['TRAJ'] == (n + 1), :],
                                 data[df['TRAJ'] == (n + 2), :]], axis=0)
        # calculate centroid of each cluster
        centroid = np.mean(points, axis=0)
        centroids.append(centroid)
        print('CENTROID:', centroid)
        cdistance = np.mean(cdist(points, [centroid]))
        spreads.append(cdistance)

    print('CENTROIDS:', centroids)
    print('SPREADS', spreads)
    dist = pdist(centroids)
    print('PAIRWISE DISTANCES:\n', squareform(dist))
    run_list = ['WT-Apo-R', 'WT-Apo-T', 'T226I-Apo-R', 'T226I-Apo-T', 'WT-Trp-R', 'WT-Trp-T',
                'T226I-Trp-R', 'T226I-Trp-T', 'WT-Tyr-R', 'WT-Tyr-T', 'T226I-Tyr-R', 'T226I-Tyr-T']

    norm_dist = squareform(dist)
    norm_dist = (norm_dist + np.min(norm_dist)) / np.max(norm_dist)
    plt.imshow(norm_dist, cmap='Blues')
    plt.xticks(np.arange(0, 12), run_list, rotation=90)
    plt.yticks(np.arange(0, 12), run_list)

    for i in range(len(run_list)):
        for j in range(len(run_list)):
            text = plt.text(j, i, round(norm_dist[i, j], 2), fontdict={'fontsize': 8},
                           ha="center", va="center", color="black")

    plt.title('PCA Centroid\nNormalized Distances')
    plt.colorbar()
    plt.tight_layout()

    if savefig != '':
        plt.savefig('savefig')
    plt.show()

    quit()
    return


def project_data(traj_list, label_list, method='pca', use_trajs=None, save='', legend=False, dist=False):
    """Project multidimensional trajectory data into 2 dimensional plot with PCA or UMAP."""

    SMALL_SIZE = 18
    MEDIUM_SIZE = 24
    BIGGER_SIZE = 32

    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    dftot = pd.concat(traj_list).reset_index(drop=True)
    nan_cols = [i for i in dftot.columns if dftot[i].isnull().any()]
    tally = 0
    dftot.loc[:, 'TRAJ'] = 0

    for n, t in enumerate(traj_list):
        dftot.loc[tally:, 'TRAJ'] = n
        tally += t.shape[0]

    sc = StandardScaler()
    if method == 'umap':
        reducer = umap.UMAP()
    else:
        reducer = PCA(n_components=2)

    if use_trajs is not None:
        # using this to keep constant PCA basis for displaying different runs
        sel_trajs = [traj_list[i] for i in use_trajs]

        # remove selected labels from label list in reverse order for legend making
        kept_labels = deepcopy(label_list)
        for index in sorted(use_trajs, reverse=True):
            del kept_labels[index]
            del traj_list[index]

        dftot = pd.concat(traj_list).reset_index(drop=True)
        nan_cols = [i for i in dftot.columns if dftot[i].isnull().any()]
        tally = 0
        dftot.loc[:, 'TRAJ'] = 0

        for n, t in enumerate(traj_list):
            dftot.loc[tally:, 'TRAJ'] = n
            tally += t.shape[0]

        xtals = pd.concat(sel_trajs).reset_index(drop=True)
        nan_cols = nan_cols + [i for i in xtals.columns if xtals[i].isnull().any()]
        xtals = xtals[[x for x in xtals.columns if x not in nan_cols]]
        sc = sc.fit(xtals)
        reducer = reducer.fit(sc.transform(xtals))
        sc = sc.transform(dftot[[c for c in dftot.columns[:-1] if c not in nan_cols]])
        reducer = reducer.transform(sc)
    else:
        # otherwise, just fit PCA on the whole thing
        kept_labels = label_list
        dftot = pd.concat(traj_list).reset_index(drop=True)
        nan_cols = [i for i in dftot.columns if dftot[i].isnull().any()]
        tally = 0
        dftot.loc[:, 'TRAJ'] = 0

        for n, t in enumerate(traj_list):
            dftot.loc[tally:, 'TRAJ'] = n
            tally += t.shape[0]

        sc = sc.fit_transform(dftot[[c for c in dftot.columns[:-1] if c not in nan_cols]])
        reducer = reducer.fit_transform(sc)

    l = []
    # default SNS color palette only has 10 colors - handle larger dataset with glasbey palette
    palette = sns.color_palette(cc.glasbey, n_colors=11)

    if dist:
        grouped_cluster_distances(dftot, reducer, save)

    # for CM-Trp (T) with xtals
    colors = [0, 0, 0, 3, 3, 3, 5, 5, 5, 1, 2]
    # colors = [0, 0, 0, 3, 3, 3, 1, 2]

    # for CM-Tyr (T) and T226I-Trp (R) with xtals
    # colors = [3, 3, 3, 4, 4, 4, 1, 2]
    for n in dftot['TRAJ'].unique():
        # special handling for xtal data points to make more visible
        if reducer[dftot['TRAJ'] == n, 0].size < 2:
            alpha = 1
            s = 100
            edgecolor = 'black'
        else:
            alpha = 0.5
            s = 30
            edgecolor = None
        l.append(plt.scatter(reducer[dftot['TRAJ'] == n, 0], reducer[dftot['TRAJ'] == n, 1], s=s,
                             c=palette[colors[n]], edgecolor=edgecolor,
                             alpha=alpha, label=kept_labels[n]))

    plt.title('Rotamer Angles', fontsize=32, fontweight='bold')
    # plt.title('Helical Angles', fontsize=32, fontweight='bold')

    plt.xlabel('PC1', fontsize=24, fontweight='bold')
    plt.ylabel('PC2', fontsize=24, fontweight='bold')
    # if legend:
    #     plt.legend(l, kept_labels, bbox_to_anchor=(1, 1))
    legend = False
    if legend:
        new_l = [l[0]] + [l[3]] + [l[6]] + [l[9]] + [l[10]]
        new_labels = ['T226I-Trp', 'CM-Tyr', 'CM-Trp', 'R (1CSM)', 'T (2CSM)']
        plt.legend(new_l, new_labels, bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.ylim(-6, 6)
    plt.xlim(-6, 6)
    plt.gca().set_aspect('equal')

    plt.gcf().set_size_inches(10, 10)
    if save != '':
        plt.savefig(save,  format='svg', dpi=1200)
    plt.show()
    return


def plot_grid(columns, traj_list, label_list, save='', xtal=True):
    """Wrapper function to plot many trajectory metrics (chi angles, distances, etc.) on same grid figure"""
    num_cols = len(columns)
    fig = plt.figure(figsize=(24, 2))
    gs = gridspec.GridSpec(1, num_cols)

    i = 0
    for target in columns:
        g0 = plot_multi_traj(traj_list, label_list, target, legend=False, xtal=xtal)
        mg0 = SeabornFig2Grid(g0, fig, gs[i])
        i += 1

    return


def plot_multi_traj(traj_list, label_list, col, skip=False, xlim=(0, 1e7), legend=False, xtal=False, alpha=0.5):
    """Jointplot of stats from N trajectories overlaid with histogram axes."""

    SMALL_SIZE = 18
    MEDIUM_SIZE = 24
    BIGGER_SIZE = 32

    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    xlim = xlim
    if 'chi' in col:
        ylims = [0, 360]
        corr_factor = 1
    elif 'H' in col or 'angle' in col:
        ylims = None
        corr_factor = 1
    elif 'subunit' in col.lower() and 'rmsd' in col.lower():
        ylims = [0, 5]
        corr_factor = 1
    else:
        ylims = [0, 20]
        corr_factor = 10
    if skip:
        n = 2
    else:
        n = 0

    if True:
        palette = sns.color_palette(cc.glasbey, n_colors=11)
    else:
        palette = sns.color_palette(n_colors=len(traj_list))

    dftot = pd.DataFrame()
    col = col.split('_')[0]
    print('TARGET:', col)
    for traj, lab in zip(traj_list, label_list):
        traj = traj.filter(regex=col)
        dftot = dftot.merge(traj, how='right', left_index=True, right_index=True, suffixes=['', lab])


    # first column doesn't get suffix for some reason - need to correct it manually
    dftot = dftot.rename(columns={col + '_chA': col + '_chA' + label_list[0],
                                  col + '_chB': col + '_chB' + label_list[0]})
    print(dftot.columns)
    # filter out all other data
    dftot = dftot * corr_factor

    print([col + lab for lab in label_list])
    dftot_A = pd.melt(dftot.reset_index(), value_vars=[col + '_chA' + lab for lab in label_list], var_name='Run', value_name=col,
                    id_vars='index')
    dftot_B = pd.melt(dftot.reset_index(), value_vars=[col + '_chB' + lab for lab in label_list], var_name='Run', value_name=col,
                    id_vars='index')
    dftot = pd.concat([dftot_A, dftot_B], axis=0)
    print(dftot.shape)

    dftot = dftot.rename(columns={'index': 'Time (ps)'})
    dftot['Time (ps)'] = dftot['Time (ps)'] / 1000.


    j = sns.jointplot(data=dftot, x='Time (ps)', y=col, hue='Run', s=20,
                      ylim=ylims, joint_kws=dict(alpha=alpha),
                      palette=[palette[5]])

    # j.ax_joint.set(xlabel='Time (ns)', ylabel='Chi', fontsize=MEDIUM_SIZE, fontweight='bold')
    # JointGrid has a convenience function
    j.set_axis_labels('Time (ns)', 'Chi', fontsize=24, fontweight='bold')

    legend = True
    if not legend:
        j.ax_joint.legend_.remove()

    if xtal:
        if 'chi' in col:
            tmp1, tmp2 = chi_angles[col + '_chA']["1csm (R)"], chi_angles[col + '_chA']["2csm (T)"]
            # do chi correction on-the-fly
            if tmp1 < 0:  tmp1 = 360 + tmp1
            if tmp2 < 0:  tmp2 = 360 + tmp2
        elif 'dist' in col:
            tmp1, tmp2 = distances[col.rstrip('_dist')]["1csm (R)"], distances[col.rstrip('_dist')]["2csm (T)"]
        else:  # subunit angle data
            tmp1, tmp2 = subunits_2csm[col]["1csm (R)"], subunits_2csm[col]["2csm (T)"]
        j.ax_joint.axhline(tmp1, ls='--', c=palette[1])
        j.ax_joint.axhline(tmp2, ls='--', c=palette[2])

    sns.move_legend(j.ax_joint, "upper left", title='', frameon=True, bbox_to_anchor=(1.1, 0.5))
    plt.savefig('legend' + '.svg', format='svg', dpi=1200, bbox_inches='tight')

    # j.ax_joint.ticklabel_format(style='sci', scilimits=(0, 0), axis='x')
    # j.ax_marg_x.remove()
    # j.fig.suptitle(col, fontsize=18, fontweight='bold')
    # j.fig.tight_layout()
    # j.fig.subplots_adjust(top=0.99)  # Reduce plot to make room
    # j.fig.set_figwidth(4)
    # j.fig.set_figheight(6)
    # plt.savefig(col + '.svg', format='svg', dpi=1200, bbox_inches='tight')
    return j


class SeabornFig2Grid:
    """Class that lets you plot figure-level seaborn plots onto a grid like matplotlib axes objects"""
    def __init__(self, seaborngrid, fig, subplot_spec):
        self.fig = fig
        self.sg = seaborngrid
        self.subplot = subplot_spec
        if isinstance(self.sg, sns.axisgrid.FacetGrid) or \
                isinstance(self.sg, sns.axisgrid.PairGrid):
            self._movegrid()
        elif isinstance(self.sg, sns.axisgrid.JointGrid):
            self._movejointgrid()
        self._finalize()

    def _movegrid(self):
        """ Move PairGrid or Facetgrid """
        self._resize()
        n = self.sg.axes.shape[0]
        m = self.sg.axes.shape[1]
        self.subgrid = gridspec.GridSpecFromSubplotSpec(n, m, subplot_spec=self.subplot)
        for i in range(n):
            for j in range(m):
                self._moveaxes(self.sg.axes[i, j], self.subgrid[i, j])

    def _movejointgrid(self):
        """ Move Jointgrid """
        h = self.sg.ax_joint.get_position().height
        h2 = self.sg.ax_marg_x.get_position().height
        r = int(np.round(h / h2))
        self._resize()
        self.subgrid = gridspec.GridSpecFromSubplotSpec(r + 1, r + 1, subplot_spec=self.subplot)

        self._moveaxes(self.sg.ax_joint, self.subgrid[1:, :-1])
        #         self._moveaxes(self.sg.ax_marg_x, self.subgrid[0, :-1])
        self._moveaxes(self.sg.ax_marg_y, self.subgrid[1:, -1])

    def _moveaxes(self, ax, gs):
        # https://stackoverflow.com/a/46906599/4124317
        ax.remove()
        ax.figure = self.fig
        self.fig.axes.append(ax)
        self.fig.add_axes(ax)
        ax._subplotspec = gs
        ax.set_position(gs.get_position(self.fig))
        ax.set_subplotspec(gs)

    def _finalize(self):
        plt.close(self.sg.fig)
        self.fig.canvas.mpl_connect("resize_event", self._resize)
        self.fig.canvas.draw()

    def _resize(self, evt=None):
        self.sg.fig.set_size_inches(self.fig.get_size_inches())
