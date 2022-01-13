from pathlib import Path
from pyhdx.fileIO import read_dynamx
from pyhdx.models import PeptideMasterTable, HDXMeasurement
from pyhdx.plot import *
import proplot as pplt

import matplotlib.pyplot as plt

cwd = Path(__file__).parent


folding_colors = ['#1f4e78', '#5a9bd5', '#9dc3e6', '#deeaf6'][::-1]
cmap = pplt.Colormap(folding_colors, listmode='discrete')

data = read_dynamx(cwd / 'data' / 'ppiAWT_4Cfolding_February2021.csv')
FD_state = 'Native folded'
s_data = data.query(f'state == "{FD_state}"')
control_FD = (FD_state, s_data['exposure'].unique()[-1])

ND_state = 'FD'
s_data = data.query(f'state == "{ND_state}"')
control_ND = (ND_state, s_data['exposure'].unique()[-1])


pmt = PeptideMasterTable(data, d_percentage=90.)
pmt.set_control(control_FD, control_ND)
state_data = pmt.get_state('folding_4C_10secLabelling')
state_data = state_data.query('exposure != 0')

hdxm = HDXMeasurement(state_data)

i = 13
fig, axes = pplt.subplots(nrows=3, aspect=3, ref=1, hratios=[6,1,1],
                          width='160mm')
cbar = peptide_coverage(axes[0], hdxm[i].data, cmap=cmap,
                        cbar_kwargs={'ticks': [0, 0.25, 0.5, 0.75, 1]})
axes[0].format(title=f'Peptide folded fraction t={hdxm.timepoints[i]:.0f} s')
cbar.set_label('Folded Fraction')

red = redundancy(axes[1], hdxm)
res = resolution(axes[2], hdxm)

kwargs = dict(width=CBAR_KWARGS['width'], extendsize=2.5, extendrect=True)
fig.colorbar(res, label='Local Resolution (residues)', **kwargs)
fig.colorbar(red, label='Redundancy (peptides)', **kwargs)


plt.savefig('ppia_4c_resolution.png')