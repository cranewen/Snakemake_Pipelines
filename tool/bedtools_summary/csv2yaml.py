from collections import defaultdict
import pandas as pd
import fileinput
import shutil
import sys, os
#sys.path.insert(0, '/mnt/storage/dept/pedonc/src/lab_pipelines')
sys.path.insert(0, '/home/yw900/lab_pipelines')
import yaml
from utility import YamlReader, YamlDict
import argparse
parser = argparse.ArgumentParser(description='Bedtool Summary Assistance Tool')
parser.add_argument('-d', '--dataset', type=str)
parser.add_argument('-c', '--csv', type=str)
parser.add_argument('-r', '--replace', type=str)
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()
from icecream import ic

def write_sample_meta_from_csv(csv_path, dataset_name):
    sample_meta_check = YamlReader('metaconfig/sample_meta.yaml')
    for smc in sample_meta_check.read_yaml():
        try:
            if smc.get(dataset_name): 
                print(f'{dataset_name} exists!')
                return 0
        except KeyError:
            print(f'Writing {dataset_name} sample meta!')

    df = pd.read_csv(csv_path, keep_default_na=False)
    ic(df) 
    idx = df['Index']
    samples = list(df.columns)[1:]
    n = len(idx)

    with open('metaconfig/sample_meta.yaml', 'a') as f:
        f.write(f'{dataset_name}:\n')
        for s in samples:
            f.write(f'  {s}:\n')
            for i in range(n):
                f.write(f'    {idx[i]}:\n') 
                if df[s][i] == '':
                    f.write(f'      null\n') 
                elif i == 0 or i == 2 or i == 4 or i == 7 or i == 8: 
                    f.write(f'      \"{df[s][i]}\"\n') 
                else:
                    f.write(f'      {df[s][i]}\n') 
        

# before copying the sample_meta.yaml, comparing with the backup yaml first
# if sample in backup file but not in sample_meta.yaml, stop and check
def check_backup():
    sample_meta_backup = YamlDict('metaconfig/sample_meta_backup.yaml')
    sample_meta = YamlDict('metaconfig/sample_meta.yaml')
    sample_meta_backup_keys = sample_meta_backup.read_yaml().keys()
    sample_meta_keys = sample_meta.read_yaml().keys()
    intersect = set(sample_meta_backup_keys) & set(sample_meta_keys)
    if sorted(list(intersect)) == sorted(list(sample_meta_backup_keys)):
        return True
    print(f'{set(sample_meta_backup_keys) - intersect} from sample_meta_backup.yaml ARE NOT in sample_meta.yaml')
    print(f'Please check the sample_meta and sample_meta_backup before overwriting sample_meta.yaml!')
    return False
    
def write_sample_meta_from_csv_replace(dataset, csv_path):
    if check_backup():
        shutil.copy2('metaconfig/sample_meta.yaml', 'metaconfig/sample_meta_backup.yaml')
        start_line = 0
        edit_dict = defaultdict()
        with open('metaconfig/sample_meta.yaml') as f:
            for n,line in enumerate(f, 1):
                if f'{dataset}:' in line:
                    start_line = n
                    print(f'replace_test line number : {n}')
        with open('metaconfig/sample_meta.yaml', 'r') as f:
            data = f.readlines()


        df = pd.read_csv(csv_path, keep_default_na=False)
        ic(df) 
        idx = df['Index']
        samples = list(df.columns)[1:]
        n = len(idx)

        l_n = start_line # line number for overwrite
        for s in samples:
            data[l_n] = f'  {s}:\n'
            l_n += 1
            for i in range(n):
                data[l_n] = f'    {idx[i]}:\n'
                l_n += 1
                if df[s][i] == '':
                    data[l_n] = f'      null\n'
                    l_n += 1
                elif i == 0 or i == 2 or i == 4 or i == 7 or i == 8: 
                    data[l_n] = f'      \"{df[s][i]}\"\n'
                    l_n += 1
                else:
                    data[l_n] = f'      {df[s][i]}\n'
                    l_n += 1
        
        #ic(data[13260]) # ChIPseq_231003 test
        #data = data[start_line:l_n]
        print(l_n)
        with open('metaconfig/sample_meta.yaml', 'w') as f:
            f.writelines(data)

def main():
    dataset = args.dataset
    csv_path = args.csv_path
    if args.replace:
        write_sample_meta_from_csv_replace(dataset, csv_path)
    else:
        write_sample_meta_from_csv(csv_path, dataset)
        
        
    
    

if __name__ == "__main__":
    main()

#write_sample_meta_from_csv(csv_path = args.csv, dataset_name = args.dataset)
#test()
#text_replace('ChIPseq_231003', 'xlsx_test/chipseq_231003_template.csv')
#write_sample_meta_from_csv_replace('modify_test', 'xlsx_test/modify_test_csv.csv')
#print(check_backup())
