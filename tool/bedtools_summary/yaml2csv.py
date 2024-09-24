import sys, os, re
#sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
#sys.path.insert(0, '/mnt/storage/dept/pedonc/src/lab_pipelines')
sys.path.insert(0, '/home/yw900/lab_pipelines')
from utility import YamlReader
from collections import defaultdict
import pandas as pd
import argparse
parser = argparse.ArgumentParser(description='Bedtool Summary Yaml to CSV tool')
parser.add_argument('-d', '--dataset', type=str)
parser.add_argument('-o', '--output', type=str)
args = parser.parse_args()

def yml2csv(dataset, out_dir):
    sm = YamlReader('metaconfig/sample_meta.yaml')
    sm_yml = defaultdict()
    out_dict = defaultdict()
    try:
        for i in sm.read_yaml():
            sm_yml = i[dataset]
    except KeyError:
        print(f'{dataset} isn\'t in sample_meta.yaml!')
        return 0
    out_dict['Index'] = ['sample_cpct_name', 's_n', 'input_sample', 'is_input',
                         'short_name', 'group_id', 'group_position', 'group_leader', 'color',
                         'is_group_leader', 'seq_paired_sorted', 'prop_paired_sorted',
                         'seq_paired_rmdup', 'prop_paired_rmdup', 'prop_paired_sorted_sp', 
                         'prop_paired_rmdup_sp']
    for k,v in sm_yml.items():
        out_dict[k] = [vi for ki,vi in v.items()]

    df = pd.DataFrame(out_dict)
    df.to_csv(out_dir, index=False)

def main():
    dataset = args.dataset
    out_dir = args.output
    yml2csv(dataset, out_dir)
    
    

if __name__ == "__main__":
    main()

