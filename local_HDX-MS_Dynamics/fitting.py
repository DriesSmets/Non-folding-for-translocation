from pyhdx.batch_processing import load_from_yaml
from pathlib import Path
from pyhdx.fitting import fit_gibbs_global_batch_aligned
import yaml
import time
from datetime import datetime
from pyhdx.local_cluster import default_client
from pyhdx.fileIO import csv_to_protein
from pyhdx.models import HDXMeasurementSet
from pyhdx.support import pprint_df_to_file
from pyhdx import VERSION_STRING

data_dir = Path('data')
yaml_stream = Path('states_dict.yaml').read_text()
data_dict = yaml.safe_load(yaml_stream)

current_dir = Path(__file__).parent
output_dir = current_dir / 'fit'
output_dir.mkdir(exist_ok=True)

client = default_client()
names = ['PPiA_WT', 'PPiB_WT']
#names = ['PPiB_WT']

alignment = {
    'PPiA_WT': 'MFKSTLAAMAAVFALSALSPAAMAAKGDPHVLLTTSAGNIELELDKQKAPVSVQNFVDYVNSGFYNNTTFHRVIPGFMIQGGGFTEQMQQKKPNPPIKNEADNGLRNTRGTIAMARTADKDSATSQFFINVADNAFLDHG---QRDFGYAVFGKVVKGMDVADKISQVPTHDVGPYQNVPSKPVVILSAKVLP',
    'PPiB_WT': '-----------------------------MVTFHTNHGDIVIKTFDDKAPETVKNFLDYCREGFYNNTIFHRVINGFMIQGGGFEPGMKQKATKEPIKNEANNGLKNTRGTLAMARTQAPHSATAQFFINVVDNDFLNFSGESLQGWGYCVFAEVVDGMDEVDKIKGVATGRSGMHQDVPKEDVIIESVTVSE'
}

hdxm_list = [load_from_yaml(dic, data_dir=data_dir, name=name) for name, dic in data_dict.items()]
rates_list = [csv_to_protein(current_dir / 'guesses' / f'{name}_rates_guess.txt')['rate'] for name in data_dict.keys()]
hdx_set = HDXMeasurementSet(hdxm_list)

hdx_set.add_alignment(list(alignment.values()))

pprint_df_to_file(hdx_set.aligned_dataframes, 'aligned_ppia_ppib.txt')
gibbs_guess = hdx_set.guess_deltaG(rates_list)

log_file = output_dir / f"fitting_log_r2.txt"
now = datetime.now()
date = f'# {now.strftime("%Y/%m/%d %H:%M:%S")} ({int(now.timestamp())})'

lines = [VERSION_STRING, date]

r1 = 0.1
r2 = 0.05

t0 = time.time()
result = fit_gibbs_global_batch_aligned(hdx_set, gibbs_guess, epochs=100000, r1=r1, r2=r2)
t1 = time.time()

block = '--------------------------'
regularizers = f'Regualizer 1: {r1}  Regualizer 2: {r2}'
loss = f'Total_loss {result.total_loss:.2f}, mse_loss {result.mse_loss:.2f}, reg_loss {result.reg_loss:.2f}' \
       f'({result.regularization_percentage:.2f}%)'
time_elapsed = f"Time elapsed: {(t1 - t0):.2f} s"
epochs = f"Number of epochs: {len(result.metadata['total_loss'])}"

result.output.to_csv(output_dir / f"fit_output_r1_{r1}_r2_{r2}.csv")#, na_rep='NaN')
result.output.to_file(output_dir / f"fit_output_r1_{r1}_r2_{r2}.txt", fmt='pprint', na_rep='NaN')

lines += ['', block, regularizers, loss, epochs, time_elapsed, block]

log_file.write_text('\n'.join(lines))

