import pandas as pd
import glob, sys, os, yaml
from collections import defaultdict
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from metaparser.deseq2parser import DESeq2Parser
from metaparser.sampleparser import SampleParser

class FileConvertor:
    def __init__(self, dataset_name):
        self.filepath = ''
        self.cell_line = []
        self.treatment = []
        self.count_matrix_column_file = ''
        self.count_matrix_column_df = ''
        self.file_list_by_cellline = []
        self.colData_dict = defaultdict()
        self.dataset_name = dataset_name

    def concat_count_files(self, filepath, ref_genome):
        file_list = glob.glob(filepath + '/*.txt')
        file_list = [f for f in file_list if 'consolidated' not in f]
        # extract all the column names from txt file names
        # here mm10 needs to be replaced by ref_genome
        file_columns_name_list = [x.split('/')[-1].split(ref_genome)[0][:-1] for x in file_list]

        index_column = ''
        columns = []
        count = 0
        for f in file_list:
            if count == 0:
                index_column = pd.read_csv(f, sep='\t', header=None).iloc[:, 0]
                columns.append(index_column)
            columns.append(pd.read_csv(f, sep='\t', header=None).iloc[:, 1])
            count += 1

        file_list_len = len(file_list)
        consolidated_dict = defaultdict()
        consolidated_dict[''] = columns[0]
        for i in range(file_list_len):
            consolidated_dict[file_columns_name_list[i]] = columns[i+1]
        
        consolidated_df = pd.DataFrame(consolidated_dict)
        consolidated_df.drop(consolidated_df.tail(5).index, inplace = True)
        if ('sorted' == filepath.split('/')[-1]):
            consolidated_df.to_csv(filepath + '/consolidated_counts_sorted.csv', index = False)
        if ('sorted_rmdup' == filepath.split('/')[-1]):
            consolidated_df.to_csv(filepath + '/consolidated_counts_sorted_rmdup.csv', index = False)

    def create_col_data_file(self, filepath):
        deseq = DESeq2Parser(self.dataset_name)
        deseq.get_col_data()
        deseq.create_col_data_df()
        df = deseq.col_data_df
        if ('sorted' == filepath.split('/')[-1]):
            df.to_csv(filepath + '/col_data_sorted.csv', index = False)
        if ('sorted_rmdup' == filepath.split('/')[-1]):
            df.to_csv(filepath + '/col_data_sorted_rmdup.csv', index = False)

    def create_deseq2_meta_yaml(self):
        samples = SampleParser()
        samples.get_sample_list(self.dataset_name)
        samples = samples.sample_list
        data = pd.read_csv("path/metaconfig/deseq2_sample_info.csv")
        print(data['Unique_sample_name'][0:6])
        

def main():
    f = FileConvertor('SLAMseq_test')
    f.create_deseq2_meta_yaml()

if __name__ == "__main__":
    main()
