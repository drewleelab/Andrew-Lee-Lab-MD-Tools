import numpy as np
import os
import pandas as pd
from tqdm import tqdm


def read_xvg(fname):
    """Read columns of data from file fname. Returns numpy array of data plus labels for axes."""
    with open(fname, 'r') as f:
        for i, line in enumerate(f):
            if not line.startswith(('#', '@')):
                skip = i
                break
            elif "title" in line:
                title = line.split('"')[-2]
            elif "xaxis" in line:
                xlabel = line.split('"')[-2]
            elif "yaxis" in line:
                ylabel = line.split('"')[-2]

    return np.genfromtxt(fname, skip_header=skip).astype(np.float32), title, xlabel, ylabel


def xvg_to_df(d, chi=False):
    """Load all xvg files in a folder into one df with a shared time index."""
    df = pd.DataFrame()
    print('Reading %s  files from folder' % str(len(os.listdir(d))))
    for f in tqdm(os.listdir(d)):
        if f.endswith('.xvg') and 'histo' not in f and 'stat' not in f and 'order' not in f and 'rmsf' not in f:
            data, title, xlabel, ylabel = read_xvg(os.path.join(d, f))
            assert not np.any(np.isnan(data))
            # add data to df as new column
            if chi:  # rescale data to [0, 360] instead of [-180, 180] for better interpretability
                data[:, 1] = correct_180(data[:, 1])
            df[f.split('.')[0]] = data[:, 1]
            df.index = data[:, 0] # keep time column as index
    return df


def correct_180(tmp):
    """Adjusts scale of chi angles from [-180, 180] to [0, 360] without distortion."""
    return [360 + t if t < 0 else t for t in tmp]


def correct_traj(data):
    """Corrects timescale of a trajectory ?"""
    ts = data[1, 0] - data[0, 0]
    data[:, 0] = np.arange(data[0, 0], data[0, 0] + (ts * data.shape[0]), ts)
    return data


def load_xtal_as_df(impute=True, subunit=False):
    """Load xtal data from 1csm & 2csm for plotting/visualization."""
    df_1csm, df_2csm = pd.DataFrame(), pd.DataFrame()

    if subunit:
        for sub in subunits_2csm.keys():
            df_1csm.loc[0, sub], df_2csm.loc[0, sub] = subunits_2csm[sub]["1csm (R)"], subunits_2csm[sub]["2csm (T)"]

    else:
        for d in distances.keys():
            df_1csm.loc[0, d + '_dist'], df_2csm.loc[0, d + '_dist'] = distances[d]["1csm (R)"] / 10., distances[d][
                "2csm (T)"] / 10.

            if df_2csm.loc[0, d + '_dist'] == 0. and impute:  # 2csm is missing a few values - may need to impute them
                df_2csm.loc[0, d + '_dist'] = df_1csm.loc[0, d + '_dist']

        for c in chi_angles.keys():
            tmp1, tmp2 = chi_angles[c]["1csm (R)"], chi_angles[c]["2csm (T)"]
            if tmp1 < 0: tmp1 = 360 + tmp1
            if tmp2 < 0: tmp2 = 360 + tmp2
            df_1csm.loc[0, c], df_2csm.loc[0, c] = tmp1, tmp2
    return df_1csm, df_2csm

# Key distances measured in xtal structures for reference
distances = {
            "234_A_157_A": {"1csm (R)": 5.2, "2csm (T)": 3.9},
            "234_B_157_B": {"1csm (R)": 5.2, "2csm (T)": 3.9},
            "23_A_157_A": {"1csm (R)": 7.2, "2csm (T)": 2.8},
            "23_B_157_B": {"1csm (R)": 7.2, "2csm (T)": 2.8},
            "23_A_234_A": {"1csm (R)": 4.8, "2csm (T)": 2.7},
            "23_B_234_B": {"1csm (R)": 4.8, "2csm (T)": 2.7},
            "23_A_208_A": {"1csm (R)": 2.8, "2csm (T)": 11.2},
            "23_B_208_B": {"1csm (R)": 2.8, "2csm (T)": 11.2},
            "23_A_204_A": {"1csm (R)": 2.7, "2csm (T)": 5.4},
            "23_B_204_B": {"1csm (R)": 2.7, "2csm (T)": 5.4},
            "215_A_204_B": {"1csm (R)": 6.5, "2csm (T)": 0.},
            "215_A_208_B": {"1csm (R)": 5.4, "2csm (T)": 0.},
            "23_A_16_A": {"1csm (R)": 13.8, "2csm (T)": 12.1},
            "23_B_16_B": {"1csm (R)": 13.8, "2csm (T)": 12.1},
            "24_A_212_A": {"1csm (R)": 3.0, "2csm (T)": 11.7},
            "24_B_212_B": {"1csm (R)": 3.0, "2csm (T)": 11.7},
            "24_A_208_A": {"1csm (R)": 3.1, "2csm (T)": 4.0},
            "24_B_208_B": {"1csm (R)": 3.1, "2csm (T)": 4.0},
            "157_A_168_A": {"1csm (R)": 9.6, "2csm (T)": 7.8},
            "157_B_168_B": {"1csm (R)": 9.6, "2csm (T)": 7.8},
            "168_A_246_A": {"1csm (R)": 3.0, "2csm (T)": 4.1},
            "168_B_246_B": {"1csm (R)": 3.0, "2csm (T)": 4.1},

}

# Key chi angles measured in xtal structures for reference
chi_angles = {
            "chi1GLU23_chA": {"1csm (R)": 80, "2csm (T)": -70},
            "chi1GLU23_chB": {"1csm (R)": 80, "2csm (T)": -70},
            "chi1TYR212_chA": {"1csm (R)": -58, "2csm (T)": -179},
            "chi1TYR212_chB": {"1csm (R)": -58, "2csm (T)": -179},
            "chi1TYR234_chA": {"1csm (R)": -108, "2csm (T)": 180},
            "chi1TYR234_chB": {"1csm (R)": -108, "2csm (T)": 180},
            "chi1GLU246_chA": {"1csm (R)": -68, "2csm (T)": -68},
            "chi1GLU246_chB": {"1csm (R)": -68, "2csm (T)": -68},
            "chi2ILE26_chA": {"1csm (R)": 174, "2csm (T)": 169},
            "chi2ILE26_chB": {"1csm (R)": 174, "2csm (T)": 169},
            "chi2ILE27_chA": {"1csm (R)": 174, "2csm (T)": -63},
            "chi2ILE27_chB": {"1csm (R)": 174, "2csm (T)": -63},
            "chi2ILE225_chA": {"1csm (R)": -54, "2csm (T)": -60},
            "chi2ILE225_chB": {"1csm (R)": -54, "2csm (T)": -60},
            "chi2LEU230_chA": {"1csm (R)": 62, "2csm (T)": 54},
            "chi2LEU230_chB": {"1csm (R)": 62, "2csm (T)": 54},
            "chi2ILE237_chA": {"1csm (R)": 174, "2csm (T)": 164},
            "chi2ILE237_chB": {"1csm (R)": 174, "2csm (T)": 164},
}

# Residue numbers corresponding to each helix in CM
helix_ids = {'H1': (5, 10),
             'H2': (11, 34),
             'H3': (39, 44),
             'H4': (58, 74),
             'H5': (75, 79),
             'H6': (107, 111),
             'H7': (113, 125),
             'H8': (125, 130),
             'H9': (139, 160),
             'H10': (160, 172),
             'H11': (172, 182),
             'H12': (184, 193),
             'H13': (194, 212),
             'H14': (226, 237),
             'H15': (237, 252)
}


# Subunit and helix angle metrics measured in xtal structures for reference
subunits_2csm = {
    "subunitB_angle": {"1csm (R)": 11.28, "2csm (T)": 0},
    "subunitB_RMSD": {"1csm (R)": 2.13, "2csm (T)": 0},
    "H1_A_H2_A": {"1csm (R)": 63.56, "2csm (T)": 63.78167837787629},
    "H1_B_H2_B": {"1csm (R)": 62.72335152, "2csm (T)": 63.781650258168845},
    "H1_A_H1_B": {"1csm (R)": 99.58908673, "2csm (T)": 93.15788556878083},
    "H2_A_H3_A": {"1csm (R)": 121.7723003, "2csm (T)": 123.29190243626056},
    "H2_B_H3_B": {"1csm (R)": 123.1879803, "2csm (T)": 123.29190210771257},
    "H2_A_H2_B": {"1csm (R)": 164.1018727, "2csm (T)": 157.4785626420974},
    "H3_A_H4_A": {"1csm (R)": 55.93330621, "2csm (T)": 63.68536545068476},
    "H3_B_H4_B": {"1csm (R)": 55.93327067, "2csm (T)": 63.6853942655176},
    "H3_A_H3_B": {"1csm (R)": 82.61502121, "2csm (T)": 93.35726569602409},
    "H4_A_H5_A": {"1csm (R)": 82.47464817, "2csm (T)": 87.11411074765424},
    "H4_B_H5_B": {"1csm (R)": 84.2918594, "2csm (T)": 87.11410170515713},
    "H4_A_H4_B": {"1csm (R)": 161.2878125, "2csm (T)": 173.3799494688368},
    "H5_A_H6_A": {"1csm (R)": 108.8003401, "2csm (T)": 114.46953891366991},
    "H5_B_H6_B": {"1csm (R)": 108.3305626, "2csm (T)": 114.46949445469369},
    "H5_A_H5_B": {"1csm (R)": 120.0601129, "2csm (T)": 129.68176441613448},
    "H6_A_H7_A": {"1csm (R)": 23.43002957, "2csm (T)": 20.183336638399457},
    "H6_B_H7_B": {"1csm (R)": 20.66846198, "2csm (T)": 20.183380629014408},
    "H6_A_H6_B": {"1csm (R)": 152.641816, "2csm (T)": 160.63805960501367},
    "H7_A_H8_A": {"1csm (R)": 41.25473603, "2csm (T)": 39.45377701962559},
    "H7_B_H8_B": {"1csm (R)": 44.40248051, "2csm (T)": 39.453849090606724},
    "H7_A_H7_B": {"1csm (R)": 118.399003, "2csm (T)": 125.28794963152114},
    "H8_A_H9_A": {"1csm (R)": 115.5724493, "2csm (T)": 117.42882166094844},
    "H8_B_H9_B": {"1csm (R)": 116.1658066, "2csm (T)": 117.42883521462876},
    "H8_A_H8_B": {"1csm (R)": 97.05706094, "2csm (T)": 99.74969609028666},
    "H9_A_H10_A": {"1csm (R)": 8.244653367, "2csm (T)": 11.255154017584998},
    "H9_B_H10_B": {"1csm (R)": 7.222508475, "2csm (T)": 11.255126622646534},
    "H9_A_H9_B": {"1csm (R)": 150.1815758, "2csm (T)": 145.44674801993207},
    "H10_A_H11_A": {"1csm (R)": 114.1855535, "2csm (T)": 113.4917688},
    "H10_B_H11_B": {"1csm (R)": 113.7078078, "2csm (T)": 113.4918225},
    "H10_A_H10_B": {"1csm (R)": 138.2296311, "2csm (T)": 125.7293342},
    "H11_A_H12_A": {"1csm (R)": 150.6368696, "2csm (T)": 150.9066053},
    "H11_B_H12_B": {"1csm (R)": 151.1705648, "2csm (T)": 150.9065670482528},
    "H11_A_H11_B": {"1csm (R)": 121.5313391, "2csm (T)": 135.93519084921104},
    "H12_A_H13_A": {"1csm (R)": 54.47062746, "2csm (T)": 58.17825611225685},
    "H12_B_H13_B": {"1csm (R)": 56.76869263, "2csm (T)": 58.17825213817734},
    "H12_A_H12_B": {"1csm (R)": 110.7791648, "2csm (T)": 127.08168943666061},
    "H13_A_H14_A": {"1csm (R)": 145.695972, "2csm (T)": 134.31493474545806},
    "H13_B_H14_B": {"1csm (R)": 145.8246801, "2csm (T)": 134.3148684986954},
    "H13_A_H13_B": {"1csm (R)": 178.6218329, "2csm (T)": 158.7612766512734},
    "H14_A_H15_A": {"1csm (R)": 12.78882613, "2csm (T)": 28.68194175446889},
    "H14_B_H15_B": {"1csm (R)": 14.24330981, "2csm (T)": 28.681938447819167},
    "H14_A_H14_B": {"1csm (R)": 112.1384253, "2csm (T)": 110.06720671322432},
    "H15_A_H15_B": {"1csm (R)": 119.8289384, "2csm (T)": 136.97756183668844}
}

