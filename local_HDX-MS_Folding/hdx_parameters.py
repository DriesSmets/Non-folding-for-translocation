#%%

# This scripts extracts RFUs folding PPiX folding local HDX-MS experiments.
# Needs PyHDX >= ef89abdf9407914fc6d79c04041d171421e42487

#%%

from pathlib import Path
from pyhdx.fileIO import read_dynamx, dataframe_to_file
from pyhdx.models import PeptideMasterTable, HDXMeasurement
from pyhdx.plot import *
import proplot as pplt

import matplotlib.pyplot as plt

cwd = Path(__file__).parent
output_dir = cwd / 'stats'
output_dir.mkdir(exist_ok=True, parents=True)

FD_state = 'FD'
exp_state = 'Folding'

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
for f in files:
    data = read_dynamx(f)

    name_, timepoint = lut[f.stem]

    s_data = data.query(f'state == "{FD_state}"')
    control_FD = (FD_state, s_data['exposure'].unique()[-1])  # = FD

    pmt = PeptideMasterTable(data, d_percentage=90.)
    pmt.set_control(control_FD)
    state_data = pmt.get_state(exp_state)
    state_data = state_data.query('exposure > 0')

    n_term = 1
    hdxm = HDXMeasurement(state_data, n_term=n_term)
    print(hdxm)

    s = str(hdxm)
    s += f"Repeatability {hdxm.data['uptake sd'].mean():.2} Da \n"

    be_frac = 1 - (hdxm[1].data['fd_uptake'] / (0.9*hdxm[1].data['maxuptake']))
    iqr = np.diff(be_frac.quantile([0.25, 0.75])).item()
    be_mean = be_frac.mean()

    s += f"Back exchange: {100*be_mean:.1f} % (IQR: {100*iqr:.1f}) \n"

    (output_dir / f"{name_}_stats.txt").write_text(s)


