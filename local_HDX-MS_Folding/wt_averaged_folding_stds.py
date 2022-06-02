#%%

# This scripts extracts RFUs folding PPiX folding local HDX-MS experiments.
# Needs PyHDX branchs stds: https://github.com/Jhsmit/PyHDX/pull/290

#%%

from pathlib import Path
from pyhdx.fileIO import read_dynamx, dataframe_to_file
from pyhdx.models import PeptideMasterTable, HDXMeasurement
from pyhdx.plot import *
import proplot as pplt

import matplotlib.pyplot as plt

cwd = Path(__file__).parent
output_dir = cwd / 'rfus'
output_dir.mkdir(exist_ok=True, parents=True)

ND_state = 'FD'
exp_state = 'Folding'


sp_length_A = 24
sp_length_B = 30

files = list((cwd / 'DynamX csv files').iterdir())

f = files[0]

#%%

#
for f in files:
    print(f.stem)
    data = read_dynamx(f)

    if f.stem.startswith('pro'):
        FD_state = 'Folding'
        control_FD = (FD_state, 0.)
    else:
        FD_state = 'Native'
        s_data = data.query(f'state == "{FD_state}"')
        control_FD = (FD_state, s_data['exposure'].unique()[-1])

    s_data = data.query(f'state == "{ND_state}"')
    control_ND = (ND_state, s_data['exposure'].unique()[-1])

    pmt = PeptideMasterTable(data, d_percentage=90.)
    pmt.set_control(control_FD, control_ND)
    state_data = pmt.get_state(exp_state)
    state_data = state_data.query('exposure != 0')

    n_term = 1
    if f.stem.lower().startswith('proppia'):
        state_data['start'] -= sp_length_A
        state_data['end'] -= sp_length_A
        n_term -= sp_length_A
    if f.stem.lower().startswith('proppib'):
        state_data['start'] -= sp_length_B
        state_data['end'] -= sp_length_B
        n_term -= sp_length_B

    hdxm = HDXMeasurement(state_data, n_term=n_term)


    txt_output = output_dir / 'txt'
    txt_output.mkdir(exist_ok=True)
    dataframe_to_file(txt_output / f"{f.stem}_rfu.csv", hdxm.rfu_residues)
    dataframe_to_file(txt_output / f"{f.stem}_rfu.txt", hdxm.rfu_residues, fmt='pprint')

    dataframe_to_file(txt_output / f"{f.stem}_rfu_sd.csv", hdxm.rfu_residues_sd)
    dataframe_to_file(txt_output / f"{f.stem}_rfu_sd.txt", hdxm.rfu_residues_sd, fmt='pprint')


    fig, ax = pplt.subplots(aspect=1.5, axwidth='120mm')
    norm = pplt.Norm('log', hdxm.timepoints.min(), hdxm.timepoints.max())
    cmap = pplt.Colormap('viridis')

    for t_idx, t in enumerate(hdxm.timepoints):
        sds = hdxm.rfu_residues_sd.iloc[:, t_idx].to_numpy()
        rfus = hdxm.rfu_residues.iloc[:, t_idx]
        bardata = np.stack([rfus - sds, rfus + sds])
        color = cmap(norm(t))
        ax.scatter(
            rfus.index,
            rfus,
            bardata=bardata,
            color=color,
            barzorder=0,
            s=7,
            barcolor='gray5')


    ax.format(ylim=(0, 1), ylabel='RFU', xlabel=r_xlabel, title=f.stem)

    locator = pplt.Locator(norm(hdxm.timepoints))
    cbar_ax = ax.colorbar(cmap, width=CBAR_KWARGS['width']*1.5, ticks=locator)
    formatter = pplt.Formatter("simple", precision=1)
    cbar_ax.ax.set_yticklabels([formatter(t) for t in hdxm.timepoints])
    cbar_ax.set_label("Folding time (s)", labelpad=-0)

    fig_output = output_dir / 'fig'
    fig_output.mkdir(exist_ok=True)

    fig.savefig(fig_output / f"{f.stem}_rfu_sd.png")
