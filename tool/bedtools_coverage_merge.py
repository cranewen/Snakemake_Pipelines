import pandas as pd
import glob, sys, os, yaml
from collections import defaultdict
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import argparse

parser = argparse.ArgumentParser(description = "Inputing a path that includes all the bedtools_coverage outputs")
parser.add_argument('-f', '--filepath', type = str, required = True)
args = parser.parse_args()

class BedtoolsMerger:
    # output csv file will be in the same directory
    def __init__(self, filepath):
        self.filepath = filepath
        

    def merge_files(self):
        file_list = glob.glob(self.filepath + '/*.txt')
        df_dict = defaultdict()
        # using the first file in the list and get fixed columns(Interval, Chr, Start, End)
        df1 = pd.read_csv(file_list[0], sep='\t', header=None)
        df_dict['Interval'] = df1.iloc[:, 3].values
        df_dict['Chr'] = df1.iloc[:, 0].values
        df_dict['Start'] = df1.iloc[:, 1].values
        df_dict['End'] = df1.iloc[:, 2].values
        for f in file_list:
            df = pd.read_csv(f, sep='\t', header=None)
            col_name = f.split('/')[-1].split('rmdup')[0][:-6]
            df_dict[col_name] = df.iloc[:, 6].values

        fixed_column_name_list = ['Interval', 'Chr', 'Start', 'End']
        file_name_list = [x.split('/')[-1].split('rmdup')[0][:-6] for x in file_list]
        df_column_name_list = fixed_column_name_list + file_name_list
        df_merged = pd.DataFrame(df_dict)
        df_merged.to_csv(self.filepath + '/tss_all_samples.csv')


def main():
    bm = BedtoolsMerger(args.filepath)
    bm.merge_files()


if __name__ == '__main__':
    main()
