import sys, os, re
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from metaparser.yamlhandle import YamlHandle
from metaparser.metaconf import DESeq2Meta
from collections import defaultdict
import pandas as pd

class DESeq2Parser:
    def __init__(self, dataset_name):
        self.dataset_name = dataset_name
        self.col_data = ''
        self.col_data_df = ''
        
    def get_col_data(self):
        yaml_handle = YamlHandle(DESeq2Meta.DESEQ2_YAML_PATH.value).read_yaml()
        for i in yaml_handle:
            self.col_data = i[self.dataset_name]

    def create_col_data_df(self):
        col_data_list = []
        column_names = ['']
        for k,v in self.col_data.items():
            temp_list = [k]
            column_names = ['']
            for v_k, v_v in v.items():
                temp_list.append(v_v)
                column_names.append(v_k)
            col_data_list.append(temp_list)
                    
        self.col_data_df = pd.DataFrame(col_data_list, columns = column_names)
        
    # This function is used with DESeq2StrategyParser to select the samples
    # that compatible with the cellline's strategy
    # The strategy argument is formated as "mm160x|mm160z"
    # It returns a string like "sample1|sample2|sample3..."
    def select_samples(self, strategy):
        sample_strategy = ''
        samples = strategy.split(",")
        strategy = [s.split("=") for s in samples]
        for k,v in self.col_data.items():
            temp = ''
            for s in strategy:
                if (v[s[0]] in s[1]) == False:
                    temp = False
                    break    
                else: 
                    temp = True
            if (temp == True):
                sample_strategy += k + '|'
        sample_strategy = sample_strategy[:-1]
        return(sample_strategy)
            

def main():
    deseq = DESeq2Parser('RNAseq_test')
    deseq.get_col_data()
    deseq.create_col_data_df()
    deseq.select_samples("cellline=A673,timepoint=4d")

if __name__ == "__main__":
    main()
