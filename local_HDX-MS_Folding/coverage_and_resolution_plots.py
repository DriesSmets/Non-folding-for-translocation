from pathlib import Path
from pyhdx.fileIO import read_dynamx
from pyhdx.models import PeptideMasterTable, HDXMeasurement
from pyhdx.plot import *
import proplot as pplt

import matplotlib.pyplot as plt

cwd = Path(__file__).parent


ND_state = 'FD'
exp_state = 'Folding'


sp_length_A = 24
sp_length_B = 30

files = list((cwd / 'DynamX csv files').iterdir())
for f in files:
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

    cmap = pplt.Colormap(['#545454'])

    fig, axes = pplt.subplots(nrows=3, aspect=3, ref=1, hratios=[6, 1, 1],
                              width='160mm')
    cbar = peptide_coverage(axes[0], hdxm[0].data, cmap=cmap, cbar=False,
                            )
    axes[0].format(title=f'Peptides')

    red = redundancy(axes[1], hdxm)

    levels = [0., 1., 2., 5., 10., 20.][::-1]
    norm = pplt.Norm('segmented', levels=levels)
    norm = pplt.DiscreteNorm(levels=levels, norm=norm)

    colors = ['008832', '72D100', 'FFFF04', 'FFB917', 'FF8923'][::-1]
    colors = ['#' + c for c in colors]
    over='#FE2B2E'
    cmap_res = pplt.Colormap(colors, discrete=True, N=len(colors), listmode='discrete')
    cmap_res.set_under(over)

    res = resolution(axes[2], hdxm, norm=norm, cmap=cmap_res)

    kwargs = dict(width=CBAR_KWARGS['width'], extendsize=2.5, extendrect=True)

    fig.colorbar(red, label='Redundancy (peptides)', **kwargs)
    fig.colorbar(res, label='Local Resolution (residues)', **kwargs)

    axes[2].format(xlabel='Residue Number')

    plt.savefig(cwd / 'figures_coverage' / f'{f.stem}.png')
    plt.close(fig)
