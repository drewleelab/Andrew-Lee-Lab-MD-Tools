import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm
from Bio.PDB import *
from Bio.PDB.Polypeptide import is_aa

import pymol
from pymol import cmd
from psico.orientation import angle_between_helices, angle_between_domains
from data_utils import helix_ids


def valid(r):
    """Check if residue is valid (not solvent or hetero res)"""
    return is_aa(r, standard=True) and 'H_' not in r.full_id[-1][0]


def get_subunit_angle_rmsd(state):
    """Get subunit info for a single state (frame)"""
    # do quick re-alignment of chain A with current state
    cmd.align('ref and c. A', '%s and c. A' % state)
    cmd.select('ref_chB', 'c. B and ref')
    cmd.select('sample_chB', '%s and c. B' % state)

    # angle between domain A and B
    angle = angle_between_domains('ref_chB', 'sample_chB')
    # rmsd of domain B once domain A is aligned
    data = cmd.align('ref_chB', 'sample_chB', cycles=0, transform=0)
    return angle, data


def subunit_rotation_pymol(ref, sample, n=1):
    """Calculate subunit rotation and RMSD"""

    cmd.load(ref, "ref")
    cmd.load(sample, "sample")
    states = cmd.count_states('sample')
    print('States in sample:', states)
    # reducing states by factor of N for speed-up
    red_states = states // n
    print('Reducing states by factor of %s to %s' % (str(n), str(red_states)))
    # we have 45 total variables (28+28+15 helix angles + subunit angle + subunit RMSD
    all_data = np.zeros((red_states, 45))

    # only need to do first alignment once w/all states
    initial = cmd.align('ref and c. A', ' sample and c. A')

    # iterate over states in object
    for state in tqdm(np.arange(red_states) + 1):
        if states == 1:
            current_sample = 'sample'
            true_state = 1
        else:
            true_state = state * n
            current_sample = 'sample_' + str(true_state).zfill(4)
        cmd.split_states('sample', first=true_state, last=true_state)

        # get subunit angle and RMSD
        angle, rmsd_data = get_subunit_angle_rmsd(current_sample)
        all_data[state - 1, 0] = angle
        all_data[state - 1, 1] = rmsd_data[0]

        # calculate helix angles and add to data matrix
        all_data, labels = calc_helix_angles(all_data, current_sample, state)
        cmd.delete(current_sample)

    cmd.delete('all')
    labels.insert(0, 'subunitB_angle')
    labels.insert(1, 'subunitB_RMSD')
    return all_data, labels


def helix_angle_robust(s1, s2):
    """If helix criteria aren't met, PyMOL will fail - use fallback cafit instead"""
    try:
        theta = angle_between_helices(s1, s2)
    except pymol.CmdException:  # if helix criteria aren't met, use fallback
        theta = angle_between_helices(s1, s2, method='cafit_orientation')
    return theta


def calc_helix_angles(all_data, current_sample, state):
    """Calculate angles between selected helix pairs"""
    # define each helix w/pymol selection algebra
    for helix in helix_ids.keys():
        idx1, idx2 = helix_ids[helix][0], helix_ids[helix][1]
        cmd.select(helix + '_A', "c. A and i. %s-%s and %s" % (str(idx1), str(idx2), current_sample))
        cmd.select(helix + '_B', "c. B and i. %s-%s and %s" % (str(idx1), str(idx2), current_sample))

    labels = []
    col = 2
    for n in range(1, 15):
        # calculate angles b/w adjacent helices & opposite subunit pair
        for combo in [('A', 'A', 1), ('B', 'B', 1), ('A', 'B', 0)]:
            s1, s2 = 'H%s_%s' % (str(n), combo[0]), 'H%s_%s' % (str(n + combo[2]), combo[1])
            theta = helix_angle_robust(s1, s2)
            # print(s1, s2, theta, col)
            all_data[state - 1, col] = theta
            labels.append(s1 + '_' + s2)
            col += 1

    # get last one manually, since loop misses 15-15 angle
    s1, s2 = 'H%s_%s' % (str(15), 'A'), 'H%s_%s' % (str(15), 'B')
    theta = helix_angle_robust(s1, s2)
    all_data[state - 1, col] = theta
    labels.append(s1 + '_' + s2)

    return all_data, labels


def main():
    """Script to calculate inter-helical angles and submit RMSD/angles via PyMOL."""
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='static or dynamic pdb to plot', type=str, default='')
    parser.add_argument('-ref', help='reference pdb for subunit comparison', type=str, default='')
    parser.add_argument('-o', help='output prefix for saved results', type=str, default='')
    parser.add_argument('-n', help='use every n-th frame for analysis (recommend using when >600 frames)',
                        type=int, default=1)
    args = parser.parse_args()

    data, labels = subunit_rotation_pymol(args.ref, args.f, args.n)

    if len(args.o) > 0:
        df = pd.DataFrame(columns=labels, data=data)
        df.to_csv(args.o, index=False)

main()
