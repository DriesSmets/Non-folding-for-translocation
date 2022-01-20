from pyhdx.batch_processing import load_from_yaml
from pathlib import Path
from pyhdx.fitting import fit_rates_weighted_average
import yaml
from pyhdx.local_cluster import default_client

data_dir = Path('data')
yaml_stream = Path('states_dict.yaml').read_text()
data_dict = yaml.safe_load(yaml_stream)

current_dir = Path(__file__).parent
#PPIA HAS THE SIGNAL PEPTIDE
#print(data_dict)
client = default_client()
names = ['PPiA_WT', 'PPiB_WT']
#names = ['PPiB_WT']

output_dir = current_dir / 'guesses'
output_dir.mkdir(exist_ok=True)


for name in names:
    print(name)
    dic = data_dict[name]
    hdxm = load_from_yaml(dic, data_dir=data_dir)

    hdxm.coverage.protein.to_file(f'{name}_sequence_info.txt', fmt='pprint')

    fr = fit_rates_weighted_average(hdxm, client=client)
    fr.output.to_file(output_dir / f'{name}_rates_guess.txt')



