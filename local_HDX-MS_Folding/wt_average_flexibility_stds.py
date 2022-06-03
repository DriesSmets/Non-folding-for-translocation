#%%

# This scripts extracts RFUs folding PPiX flexibility local HDX-MS experiments.
# Needs PyHDX branchs stds: https://github.com/Jhsmit/PyHDX/pull/290


from pathlib import Path
from pyhdx.fileIO import read_dynamx, dataframe_to_file
from pyhdx.models import PeptideMasterTable, HDXMeasurement
from pyhdx.plot import *
import proplot as pplt

import matplotlib.pyplot as plt

cwd = Path(__file__).parent
output_dir = cwd / 'flexibility'
output_dir.mkdir(exist_ok=True, parents=True)

FD_state = 'FD'
ND_state = 'Folding'


sp_length_A = 24
sp_length_B = 30

files = list((cwd / 'DynamX csv files').iterdir())

#%%

def process_rfu(name, data, control_FD, control_ND, suffix, exp_state):
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

for f in files:
    data_ = read_dynamx(f)

    # timepoint information is not used here
    name_, timepoint = lut[f.stem]

    s_data = data_.query(f'state == "{FD_state}"')
    control_FD_ = (FD_state, s_data['exposure'].unique()[-1])

    s_data = data_.query(f'state == "{ND_state}"')
    control_ND_ = (ND_state, s_data['exposure'].unique()[0])

    process_rfu(name_, data_, control_FD_, control_ND_, 'folding_flexibility', 'Folding')
    if name_.startswith("Ppi"):
        process_rfu(name_, data_, control_FD_, control_ND_, 'native_flexibility', 'Native')
