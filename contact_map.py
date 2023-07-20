from Bio.PDB import *
from Bio.PDB.Polypeptide import is_aa
from scipy.spatial import distance_matrix
import numpy as np
import matplotlib.pyplot as plt
import argparse
from tqdm import tqdm


def valid(r):
    """Check if residue is valid (not solvent or hetero res)"""
    return is_aa(r, standard=True) and 'H_' not in r.full_id[-1][0]


def contact_map(model, thresh=10., method='Ca'):
    """Calculate contact map for a single model (frame)"""
    valid_chains = [c for c in model.get_chains()][0:2]  # filter out extra Trp/Tyr ligands
    valid_res = []
    for vc in valid_chains:
        for r in vc.get_residues():
            if valid(r):
                valid_res.append(r)

    num_res = len(valid_res)

    if method == 'Ca':
        atoms = np.zeros((num_res, 3))
        for n, vr in enumerate(valid_res):
            atoms[n, :] = vr['CA'].get_coord()
        return distance_matrix(atoms, atoms) <= thresh

    # use nearest-atom distance instead of just Ca distances
    else:
        total_mat = np.zeros((num_res, num_res))

        # load all atoms and mark which ones belong to which residue
        for n1, vr1 in enumerate(valid_res):
            for n2, vr2 in enumerate(valid_res):
                if n1 == n2:  # handle diagonals trivially
                    total_mat[n1, n2] = 0.
                else:
                    atoms1 = [a.get_coord() for a in vr1.get_atoms()]
                    atoms2 = [a.get_coord() for a in vr2.get_atoms()]
                    mat = distance_matrix(atoms1, atoms2)
                    total_mat[n1, n2] = np.min(mat)

        return total_mat <= thresh


def prob_contact_map(models, thresh=10., method='Ca'):
    """Calculate probability-of-contact map for a trajectory"""
    valid_chains = [c for c in models[0].get_chains()][0:2]  # filter out extra Trp/Tyr ligands
    valid_res = []
    for vc in valid_chains:
        for r in vc.get_residues():
            if valid(r):
                valid_res.append(r)
    num_res = len(valid_res)
    net_map = np.zeros((num_res, num_res), dtype=np.int32)
    i = 1

    for i, m in tqdm(enumerate(models)):
        net_map += contact_map(m, thresh, method)
        plt.show()

    return np.abs(net_map / float(i))


def animated_contact_map(models, thresh=10., method='Ca'):
    """Animate and save contact map video looping over many frames."""
    num_res = len([r for r in models[0].get_residues() if valid(r)])
    import matplotlib.animation as animation
    frames = []
    fig = plt.figure()
    for i, m in tqdm(enumerate(models)):
        frames.append([plt.imshow(contact_map(m, thresh, method), cmap='Greys', animated=True)])

    ani = animation.ArtistAnimation(fig, frames, interval=50, blit=True, repeat_delay=True)
    ani.save('animation.mp4')
    plt.show()
    return


def main():
    """Script for plotting contact maps from multi-frame PDB files."""

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='input .pdb file (any # of frames)', type=str)
    parser.add_argument('-t', help='threshold for distance matrix calculation', default=10., type=float)
    parser.add_argument('-o', help='output file to save distance matrix (optional)', default='', type=str)
    parser.add_argument('-m', help='method of contact mapping (Ca, nearest)', default='Ca', type=str)
    parser.add_argument('-r', help='reference pdb (static) for difference map calculation', default='', type=str)
    args = parser.parse_args()

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("raw", args.f)
    models = [m for m in structure.get_models()]

    # if only one model exists, it's a static structure
    if len(models) == 1:
        print('Static Structure Detected')
        mat = contact_map(models[0], args.t, args.m)
        cmap = 'Greys'
        if len(args.o) > 0:
            np.save(args.o, mat)

    # if many models exist, it's a trajectory
    else:
        print('Trajectory Detected')
        mat = prob_contact_map(models, args.t, args.m)
        cmap = 'plasma'

    plt.imshow(mat, cmap=cmap, origin='lower')
    if cmap == 'plasma':
        plt.colorbar()
    plt.title('Protein Contact Map at ' + str(args.t) + 'A')
    plt.show()

    # calculate reference structure and make difference map
    if len(args.r) > 0:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("raw", args.r)
        models = [m for m in structure.get_models()]
        # if only one model exists, it's a static structure
        if len(models) == 1:
            print('Static Structure Detected')
            ref = contact_map(models[0], args.t, args.m)

        # if many models exist, it's a trajectory
        else:
            print('Trajectory Detected')
            ref = prob_contact_map(models, args.t, args.m)

        diff = np.abs(mat.astype(np.int32) - ref.astype(np.int32))
        # save data if requested
        if len(args.o) > 0:
            np.save(args.o, mat)
            np.save('diff_' + args.o, diff)

        plt.imshow(diff, cmap=cmap, origin='lower')
        plt.colorbar()
        plt.title('Difference Map at ' + str(args.t) + 'A')
        plt.show()

main()
