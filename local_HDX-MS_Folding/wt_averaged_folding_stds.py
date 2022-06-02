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

#%%

# format is <f.stem>: (<label>, <100% folded timepoint in minutes>)
lut = {
    'ppiAWT_25CFolding_February2021_new': ("PpiA 25C", 30),
    "ppiAWT_4Cfolding_February2021_new": ("PpiA 4C", 960),
    "ppiBWT_4Cfold_Febr2021_new": ("PpiB 4C", 30),
    "ppiB_WT_folding_25C_Febr2021_new": ("PpiB 25C", 30),
    "proppiAWT_25CFolding_April2021_2": ("pro PpiA 25C", 1440),
    "proppiBcorrect_25C_Aug2021_2": ("pro PpiB 25C", 1440)
}


#%%

def process_rfu(name, data, control_FD, control_ND, suffix):
    pmt = PeptideMasterTable(data, d_percentage=95.)
    pmt.set_control(control_FD, control_ND)
    state_data = pmt.get_state(exp_state)
    state_data = state_data.query('exposure > 0')

    n_term = 1
    hdxm = HDXMeasurement(state_data, n_term=n_term)

    rfu = hdxm.rfu_residues
    rfu_sd = hdxm.rfu_residues_sd

    if name.startswith('pro PpiA'):
        shift = sp_length_A
    elif name.startswith('pro PpiB'):
        shift = sp_length_B
    else:
        shift = 0

    rfu.index -= shift
    rfu_sd.index -= shift
    txt_output = output_dir / 'txt'
    txt_output.mkdir(exist_ok=True)
    dataframe_to_file(txt_output / f"{name.replace(' ', '_')}_rfu_{suffix}.csv", rfu)
    dataframe_to_file(txt_output / f"{name.replace(' ', '_')}_rfu_{suffix}.txt", rfu, fmt='pprint')

    dataframe_to_file(txt_output / f"{name.replace(' ', '_')}_rfu_sd_{suffix}.csv", rfu_sd)
    dataframe_to_file(txt_output / f"{name.replace(' ', '_')}_rfu_sd_{suffix}.txt", rfu_sd, fmt='pprint')

    fig, ax = pplt.subplots(aspect=1.5, axwidth='120mm')
    norm = pplt.Norm('log', hdxm.timepoints[1], hdxm.timepoints[-1])
    cmap = pplt.Colormap('viridis')

    for t_idx, t in enumerate(hdxm.timepoints[1:]):
        sds = rfu_sd.iloc[:, t_idx].to_numpy()
        rfus = rfu.iloc[:, t_idx]
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

    ax.format(ylim=(0, 1), ylabel='Fraction folded', xlabel=r_xlabel, title=name)

    locator = pplt.Locator(norm(hdxm.timepoints))
    cbar_ax = ax.colorbar(cmap, width=CBAR_KWARGS['width']*1.5, ticks=locator)
    formatter = pplt.Formatter("simple", precision=1)
    cbar_ax.ax.set_yticklabels([formatter(t) for t in hdxm.timepoints])
    cbar_ax.set_label("Folding time (s)", labelpad=-0)

    fig_output = output_dir / 'fig'
    fig_output.mkdir(exist_ok=True)

    fig.savefig(fig_output / f"{name.replace(' ', '_')}_rfu_sd_{suffix}.png")
    fig.savefig(fig_output / f"{name.replace(' ', '_')}_rfu_sd_{suffix}.pdf")
    pplt.close(fig)

#%%
for f in files:
    data_ = read_dynamx(f)

    name_, timepoint = lut[f.stem]

    s_data = data_.query(f'state == "{ND_state}"')
    control_ND_ = (ND_state, s_data['exposure'].unique()[-1])  # = FD

    # Native as 100% folded
    FD_state = 'Native'
    s_data = data_.query(f'state == "{FD_state}"')
    print('native', name_, s_data['exposure'].unique())
    control_FD_ = (FD_state, s_data['exposure'].unique()[-1])
    print(name_, control_FD_, control_ND_)
    process_rfu(name_, data_, control_FD_, control_ND_, 'native')

    # last timepoint as 100% folded
    FD_state = 'Folding'
    s_data = data_.query(f'state == "{FD_state}"')
    t_bool = s_data['exposure'].unique().round() == timepoint*60
    FD_exposure = float(s_data['exposure'].unique()[t_bool])

    print(name_, timepoint*60, FD_exposure)
    control_FD_ = (FD_state, FD_exposure)
    data_ = data_.query(f'exposure <= {FD_exposure}')
    process_rfu(name_, data_, control_FD_, control_ND_, 'timepoint')
