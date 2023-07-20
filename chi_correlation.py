import matplotlib.pyplot as plt
import argparse
from scipy.stats import pearsonr, spearmanr
import numpy as np
from data_utils import xvg_to_df


def main(args):
    """Script to calculate chi correlations between a directory of individual chi xvg files"""
    used = xvg_to_df(args.xvg_dir, chi=True)

    flag = np.zeros((len(used.columns), len(used.columns)))
    for n, c in enumerate(used.columns):
        for n2, c2 in enumerate(used.columns):
            if args.corr_method == 'pearson':
                s, p = pearsonr(used[c], used[c2])
            else:
                s, p = spearmanr(used[c], used[c2])
            s = np.abs(s)
            flag[n, n2] = s

    plt.imshow(flag, cmap='Greys', vmin=0, vmax=1)
    plt.title(args.title)
    plt.xticks(np.arange(0, len(used.columns)), used.columns, rotation=90)
    plt.yticks(np.arange(0, len(used.columns)), used.columns)
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(f"{args.title}_chi_correlation.png")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--xvg_dir', help='directory of xvg files (one per chi angle) for analysis')
    parser.add_argument('--title', help='Name of run for display purposes')
    parser.add_argument('--corr_method', help='which correlation to use', default='pearson')
    args = parser.parse_args()
    main(args)