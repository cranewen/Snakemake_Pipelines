import numpy as np
import pandas as pd
import argparse
from scipy import stats
import sys, os, re
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

parser = argparse.ArgumentParser(description = "Inputing a file's path and name")
parser.add_argument('-f', '--filepath', type = str, required = True)
args = parser.parse_args()
np.seterr(divide='ignore', invalid='ignore')


def zscore(filepath):
    df = pd.read_csv(filepath)
    ensembl_ids = df.iloc[:, 0]
    gene_names = df.iloc[:, -2]
    biotype = df.iloc[:, -1]
    #print(protein_coding)
    column_names = df.columns[1:-2].values
    print(column_names)
    df_np_arr = df.iloc[:, 1:-2].to_numpy()
    #print(df_np_arr)
    z_scores = stats.zscore(df_np_arr, axis=1)
    new_df = pd.DataFrame(z_scores, columns = column_names)
    new_df.insert(0, 'gene', ensembl_ids)
    new_df['gene'] = gene_names
    new_df['biotype'] = biotype
    zscore_path = filepath.split('rlog_expression')[0] + 'rlog_expression_zscore/'
    zscore_filename = filepath.split('/')[-1].split('.')[0] + '_zscore.csv'
    new_df.to_csv(zscore_path + zscore_filename, index = False)


zscore(args.filepath)

