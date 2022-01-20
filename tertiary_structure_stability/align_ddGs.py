#%%
import numpy as np
import pandas as pd
from pathlib import Path
from tertiary_structure_stability.hal.functions import align_dataframes
import proplot as pplt
import matplotlib.pyplot as plt
from pyhdx.support import apply_cmap, series_to_pymol

#%%
cwd = Path().resolve() / 'tertiary_structure_stability'

#%%
alignments = {
'ppia':      'AKGDPHVLLTTSAGNIELELDKQKAPVSVQNFVDYVNSGFYNNTTFHRVIPGFMIQGGGFTEQMQQKKPNPPIKNEADNGLRNTRGTIAMARTADKDSATSQFFINVADNAFLDHG---QRDFGYAVFGKVVKGMDVADKISQVPTHDVGPYQNVPSKPVVILSATVLP',
'ppib':      '-----MVTFHTNHGDIVIKTFDDKAPETVKNFLDYCREGFYNNTIFHRVINGFMIQGGGFEPGMKQKATKEPIKNEANNGLKNTRGTLAMARTQAPHSATAQFFINVVDNDFLNFSGESLQGWGYCVFAEVVDGMDEVDKIKGVATGRSGMHQDVPKEDVIIESVTVSE',
}

#%%
dA = pd.read_csv(cwd / "1v9t_SSM_ddg_NAN.csv", header=None, index_col=0).T
dA['sequence'] = list(alignments['ppia'].replace('-', ''))
dB = pd.read_csv(cwd / "1lop_SSM_ddg_NAN.csv", header=None, index_col=0).T
dB['sequence'] = list(alignments['ppib'].replace('-', ''))

# align dataframes
aligned_df = align_dataframes({'ppia': dA, 'ppib': dB}, alignments)

# Check resuling alignment by looking at the sequences
aligned_sequences = aligned_df.xs('sequence', level=1, axis=1)


#%%
aligned_df[('ppia', 'sequence')] = aligned_df[('ppia', 'sequence')].fillna('')
aligned_df[('ppib', 'sequence')] = aligned_df[('ppib', 'sequence')].fillna('')

#%%
def make_index(sub_df):
    lst = [f'{a}\n{int(n)}' if n % 10 == 0 else a for a, n in zip(sub_df['sequence'], sub_df['r_number'].astype(float))]
    return pd.Index(lst)

a_index = make_index(aligned_df['ppia'])
b_index = make_index(aligned_df['ppib'])

#%%

pplt.rc.update({'tick.labelsize': 5})
ppia_data = aligned_df['ppia'].drop(['r_number', 'sequence'], axis=1).iloc[:, ::-1].set_index(a_index)
ppib_data = aligned_df['ppib'].drop(['r_number', 'sequence'], axis=1).iloc[:, ::-1].set_index(b_index)

fig, axes = pplt.subplots(nrows=3, aspect=6, width=200/25.4, share=0)

collectionA = axes[0].pcolor(ppia_data.T, robust=True, edgefix=False)
collectionB = axes[1].pcolor(ppib_data.T, cmap=collectionA.cmap, norm=collectionA.norm, edgefix=False)

axes[0].format(xlabel='PpiA sequence')
axes[0].set_facecolor('#8c8c8c')
axes[1].format(xlabel='PpiB sequence')
axes[1].set_facecolor('#8c8c8c')


axes[0].colorbar(collectionA, label='ddG')
axes[1].colorbar(collectionB, label='ddG')

combined_sequence = [f'{sa}\n{sb}' for sa, sb in zip(aligned_df[('ppia', 'sequence')], aligned_df[('ppib', 'sequence')])]

idx = pd.Index(combined_sequence)
ddGs = aligned_df['ppia'].drop(['sequence', 'r_number'], axis=1) - aligned_df['ppib'].drop(['sequence', 'r_number'], axis=1)
ddGs = ddGs.set_index(idx).iloc[:, ::-1]
cmap = pplt.Colormap('PiYG')
cmap.set_bad('#8c8c8c')
levels = np.arange(-12, 14, 2)
norm = pplt.Norm('linear', vmin=-12, vmax=12)
norm = pplt.DiscreteNorm(levels=levels, norm=norm, clip=True)

collection = axes[2].pcolor(ddGs.T, cmap=cmap, norm=norm)
axes[2].colorbar(collection, label='dddG', width=0.1)
axes[2].format(title='PpiA - PpiB')
axes[2].set_facecolor('#8c8c8c')


plt.savefig('ppia_ppib_tertiary_structure.png')

plt.show()


#%%

fig, axes = pplt.subplots(ncols=3, share=False)
axes[0].hist(ppia_data.to_numpy().flatten(), bins='fd')
axes[0].format(xlabel='ddG', ylabel='Occurence', title='PpiA')
axes[1].hist(ppib_data.to_numpy().flatten(), bins='fd')
axes[1].format(xlabel='ddG', ylabel='Occurence', title='PpiB')
axes[2].hist(ddGs.to_numpy().flatten(), bins='fd')
axes[2].format(xlabel='dddG', ylabel='Occurence', title='PpiA - PpiB')

plt.savefig('histograms.png')
plt.show()


#%%

ppia_core_indices = [9, 18, 20, 33, 36, 88, 135, 162]
ppib_core_indices = [4, 13, 15, 28, 31, 83, 133, 160]

idx = aligned_df[('ppia', 'r_number')][ppia_core_indices].index
sub_df = aligned_df.loc[idx]

zipper = (sub_df[('ppia', 'sequence')], sub_df[('ppib', 'sequence')], ppia_core_indices, ppib_core_indices)
combined_sequence = [f'{sa}\n{sb}\n{ia}\n{ib}' for sa, sb, ia, ib in zip(*zipper)]

idx = pd.Index(combined_sequence)
ddGs_sub = sub_df['ppia'].drop(['sequence', 'r_number'], axis=1) - sub_df['ppib'].drop(['sequence', 'r_number'], axis=1)
ddGs_sub = ddGs_sub.set_index(idx).iloc[:, ::-1]
cmap = pplt.Colormap('PiYG')
cmap.set_bad('#8c8c8c')
levels = np.arange(-12, 14, 2)
norm = pplt.Norm('linear', vmin=-12, vmax=12)
norm = pplt.DiscreteNorm(levels=levels, norm=norm, clip=True)

fig, axes = pplt.subplots(nrows=1, aspect=1, width=100/25.4, share=0)

collection = axes[0].pcolor(ddGs_sub.T, cmap=cmap, norm=norm, diverging=True)
axes[0].colorbar(collection, label='dddG', width=0.1)
axes[0].format(title='PpiA - PpiB')
axes[0].set_facecolor('#8c8c8c')

plt.savefig('core_dddG.png')
plt.show()

#%%


#%%

ddGs = aligned_df['ppia'].drop(['sequence', 'r_number'], axis=1) - aligned_df['ppib'].drop(['sequence', 'r_number'], axis=1)
ddG_arr = ddGs.to_numpy()

ddG_values = ddG_arr[~np.isnan(ddG_arr)]


#%%
from scipy.optimize import fsolve

def get_frac(cutoff):
    frac = (np.abs(ddG_values) < cutoff[0]).sum() / ddG_values.size
    return frac - 0.95

thd = fsolve(get_frac, np.array([12]))[0]
thd
#%%

projected = np.nanmean(np.clip(ddG_arr, -thd, thd), axis=1)

projected_df = pd.DataFrame(
    {'r_number': aligned_df[(('ppia', 'r_number'))],
     'sequence': aligned_df[(('ppia', 'sequence'))],
     'projected_dddG': projected
     }
)



#%%
projected_df
#%%
std = projected_df['projected_dddG'].std()
std

#%%
dropped = projected_df.dropna()
dropped


#%%
dropped['projected_dddG'].abs() > 2*std
#%%

outliers = dropped[dropped['projected_dddG'].abs() > 2*std]
outliers

#%%
pplt.rc.update({'tick.labelsize': 12})  # with rc context ...

fig, ax = pplt.subplots(aspect=3, width=160/25.4)
ax.plot(dropped['r_number'].to_numpy().astype(int), dropped['projected_dddG'].to_numpy())
ax.scatter(outliers['r_number'].to_numpy().astype(int), outliers['projected_dddG'], color='r')
ax.axhline(std, color='k', alpha=0.5, linestyle='--')
ax.axhline(2*std, color='k', alpha=0.5, linestyle='--')
ax.axhline(-std, color='k', alpha=0.5, linestyle='--')
ax.axhline(-2*std, color='k', alpha=0.5, linestyle='--')

for i in outliers.index:
    x = outliers.loc[i]['r_number']
    y = outliers.loc[i]['projected_dddG']
    s = outliers.loc[i]['sequence']
    ax.annotate(f'{x}{s}', (x, y))
ax.format(xlabel='Residue Number')
plt.savefig('ppia_projected.png')
plt.show()

#%%
# plot projected dddG values on ppia / ppib residues
for prot in ['ppia', 'ppib']:
    projected_df = pd.DataFrame(
        {'r_number': aligned_df[((prot, 'r_number'))],
         'sequence': aligned_df[((prot, 'sequence'))],
         'projected_dddG': projected
         }
    )

    std = projected_df['projected_dddG'].std()
    dropped = projected_df.dropna()

    outliers = dropped[dropped['projected_dddG'].abs() > 2*std]
    pplt.rc.update({'tick.labelsize': 12})  # with rc context ...

    fig, ax = pplt.subplots(aspect=3, width=160/25.4)
    ax.plot(dropped['r_number'].to_numpy().astype(int), dropped['projected_dddG'].to_numpy())
    ax.scatter(outliers['r_number'].to_numpy().astype(int), outliers['projected_dddG'], color='r')
    ax.axhline(std, color='k', alpha=0.5, linestyle='--')
    ax.axhline(2*std, color='k', alpha=0.5, linestyle='--')
    ax.axhline(-std, color='k', alpha=0.5, linestyle='--')
    ax.axhline(-2*std, color='k', alpha=0.5, linestyle='--')

    for i in outliers.index:
        x = outliers.loc[i]['r_number']
        y = outliers.loc[i]['projected_dddG']
        s = outliers.loc[i]['sequence']
        ax.annotate(f'{x}{s}', (x, y))
    ax.format(xlabel='Residue Number')
    plt.savefig(f'{prot}_projected.png')
    plt.show()

    cmap = pplt.Colormap('PiYG')
    cmap.set_bad('#8c8c8c')
    levels = np.arange(-12, 14, 2)
    norm = pplt.Norm('linear', vmin=-12, vmax=12)
    norm = pplt.DiscreteNorm(levels=levels, norm=norm, clip=True)

    value_series = dropped.set_index('r_number')['projected_dddG']
    colors = apply_cmap(value_series, cmap, norm)
    pymol_str = series_to_pymol(colors)

    Path(f'{prot}_colors.pml').write_text(pymol_str)




