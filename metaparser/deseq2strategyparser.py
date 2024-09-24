import sys, os, re
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from metaparser.yamlhandle import YamlHandle
from metaparser.deseq2parser import DESeq2Parser
from metaparser.metaconf import DESeq2StrategyMeta
from collections import defaultdict

class DESeq2StrategyParser:
    def __init__(self, dataset_name):
        self.dataset_name = dataset_name
        # a list of all the names of the csv files
        self.deseq2_result = ''
        self.deseq2_strategy = defaultdict()
        # get meta data for repective dataset
        self.deseq2_meta = ''

    def get_deseq2_strategy(self):
        yaml_handle = YamlHandle(DESeq2StrategyMeta.DESEQ2_STRATEGY_META_YAML_PATH.value).read_yaml()       
        for i in yaml_handle:
            self.deseq2_strategy = i[self.dataset_name]
        self.deseq2_result = [k for k in self.deseq2_strategy.keys()]
        temp = self.deseq2_strategy
        temp_meta = DESeq2Parser(self.dataset_name)
        temp_meta.get_col_data()
        for k,v in temp.items():
            self.deseq2_strategy[k] = '\"' + temp_meta.select_samples(v[0]) + '\" ' +'\"' + '\" \"'.join(v[1:]) + '\"'
            

def main():
    test = DESeq2StrategyParser('RNAseq_test')
    test.get_deseq2_strategy()
    print(test.deseq2_strategy)

if __name__ == "__main__":
    main()
            
