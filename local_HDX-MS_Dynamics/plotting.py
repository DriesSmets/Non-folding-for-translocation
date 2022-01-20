from functions import *
from pyhdx.fileIO import csv_to_protein, csv_to_dataframe
from pyhdx.alignment import align_dataframes
from pyhdx.support import pprint_df_to_file
import proplot as pplt
import scipy

from pathlib import Path
import pandas as pd

current_dir = Path(__file__).parent
r1, r2 = 0.1, 0.05
width = 100/25.4

names = ['PPiA_WT', 'PPiB_WT']

all_alignments = {
    'ppia_nma':      '------------------------AKGDPHVLLTTSAGNIELELDKQKAPVSVQNFVDYVNSGFYNNTTFHRVIPGFMIQGGGFTEQMQQKKPNPPIKNEADNGLRNTRGTIAMARTADKDSATSQFFINVADNAFLDHG---QRDFGYAVFGKVVKGMDVADKISQVPTHDVGPYQNVPSKPVVILSATVLP',
    'ppib_nma':      '-----------------------------MVTFHTNHGDIVIKTFDDKAPETVKNFLDYCREGFYNNTIFHRVINGFMIQGGGFEPGMKQKATKEPIKNEANNGLKNTRGTLAMARTQAPHSATAQFFINVVDNDFLNFSGESLQGWGYCVFAEVVDGMDEVDKIKGVATGRSGMHQDVPKEDVIIESVTVSE',
    'ppia_hdx':      'MFKSTLAAMAAVFALSALSPAAMAAKGDPHVLLTTSAGNIELELDKQKAPVSVQNFVDYVNSGFYNNTTFHRVIPGFMIQGGGFTEQMQQKKPNPPIKNEADNGLRNTRGTIAMARTADKDSATSQFFINVADNAFLDHG---QRDFGYAVFGKVVKGMDVADKISQVPTHDVGPYQNVPSKPVVILSAKVLP',
    'ppib_hdx':      '-----------------------------MVTFHTNHGDIVIKTFDDKAPETVKNFLDYCREGFYNNTIFHRVINGFMIQGGGFEPGMKQKATKEPIKNEANNGLKNTRGTLAMARTQAPHSATAQFFINVVDNDFLNFSGESLQGWGYCVFAEVVDGMDVVDKIKGVATGRSGMHQDVPKEDVIIESVTVSE'

}


gibbs_dfs = csv_to_protein(current_dir / 'fit' / f"fit_output_r1_{r1}_r2_{r2}.csv", column_depth=2)
nma_dfs = {name: csv_to_dataframe(current_dir / 'normal_modes' / f'{name}_NMA.txt').set_index('r_number') for name in names}

dfs_dict = {
    'ppia_nma': nma_dfs['PPiA_WT'],
    'ppib_nma': nma_dfs['PPiB_WT'],
    'ppia_hdx': gibbs_dfs['PPiA_WT'],
    'ppib_hdx': gibbs_dfs['PPiB_WT']
}

for key, df in dfs_dict.items():
    pprint_df_to_file(df, f'{key}.txt')

aligned_dataframes = align_dataframes(dfs_dict, all_alignments)
aligned_dataframes.index -= 24

bools1 = aligned_dataframes.columns.get_level_values(1) == 'r_number'
bools2 = aligned_dataframes.columns.get_level_values(1) == 'sequence'
bools = np.logical_or(bools1, bools2)
selected_df = aligned_dataframes.iloc[:, bools]
pprint_df_to_file(selected_df, 'All_data_aligned.txt')


names = ['PpiA', 'PpiB']
fig, axes = pplt.subplots(nrows=2, aspect=2.5, width=width,sharex=False)
for name, ax in zip(names, axes):
    print(name)
    df = pd.concat([aligned_dataframes[name.lower() + '_hdx', 'deltaG'], aligned_dataframes[name.lower() + '_nma', 'displacement']],
                   axis=1, keys=['deltaG', 'displacement'])
#
    na_removed = df.dropna(how='any')
    x = na_removed['displacement']
    y = na_removed['deltaG']
    rho, p = scipy.stats.pearsonr(x, y)
    print('rho, p:', rho, p)
#
    ax.plot(aligned_dataframes.index, aligned_dataframes[name.lower() + '_nma', 'displacement'],
              color='magenta', label='Displacement', zorder=-10)
    ax1 = ax.twinx()
    single_deltaG_scatter(ax1, aligned_dataframes[name.lower() + '_hdx'])
    dG = aligned_dataframes[name.lower() + '_hdx', 'deltaG']
    print(dG.mean(), dG.min(), dG.max(), dG.sum())

    ax1.set_ylim(40, 5)

    single_deltaG_inverted_colorbar(ax1, ax1, [10, 25, 40])
    title = name
    ax.format(title=title, ylabel='NMA Displacement', xlabel='Alignment Index')
    ax.yaxis.label.set_color('magenta')
    ax.tick_params(colors='magenta', axis='y', which='both')

output = 'save'

if output == 'show':
    plt.show()
elif output == 'save':

    fname = 'PpiA_PpiB_aligned'
    plt.savefig(f'{fname}.png', transparent=False)
    plt.savefig(f'{fname}.pdf', transparent=False)
    #plt.savefig(f'{fname}.eps', transparent=False)