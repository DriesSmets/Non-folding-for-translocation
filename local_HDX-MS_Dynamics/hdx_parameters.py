#%%

# This scripts extracts RFUs folding PPiX folding local HDX-MS experiments.
# Needs PyHDX >= ef89abdf9407914fc6d79c04041d171421e42487

#%%

from pathlib import Path
from pyhdx.fileIO import read_dynamx, dataframe_to_file
from pyhdx.models import PeptideMasterTable, HDXMeasurement
from pyhdx.batch_processing import StateParser
from pyhdx.plot import *
import proplot as pplt
import yaml

import matplotlib.pyplot as plt

cwd = Path(__file__).parent

data_dir = (cwd / 'data')
state_spec = yaml.safe_load((cwd / 'states_pyhdx_v040.yaml').read_text())

output_dir = cwd / 'stats'
output_dir.mkdir(exist_ok=True, parents=True)

state_spec

#%%

parser = StateParser(state_spec, data_src=data_dir)
hdxm_set = parser.load_hdxmset()



#%%
for hdxm in hdxm_set:
    print(hdxm)

    s = str(hdxm)
    s += f"Repeatability {hdxm.data['uptake sd'].mean():.2} Da \n"

    be_frac = 1 - (hdxm[1].data['fd_uptake'] / (0.9*hdxm[1].data['maxuptake']))
    iqr = np.diff(be_frac.quantile([0.25, 0.75])).item()
    be_mean = be_frac.mean()

    s += f"Back exchange: {100*be_mean:.1f} % (IQR: {100*iqr:.1f}) \n"

    (output_dir / f"{hdxm.name}_stats.txt").write_text(s)


