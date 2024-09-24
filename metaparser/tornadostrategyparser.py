import sys, os, re
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from metaparser.yamlhandle import YamlHandle
from metaparser.tornadoparser import TornadoParser
from metaparser.genomeannotationsparser import GenomeAnnotationsParser
from metaparser.metaconf import TornadoStrategyMeta
from collections import defaultdict

class TornadoStrategyParser:
    def __init__(self, dataset_name):
        self.dataset_name = dataset_name
        # Taking the bam output files with corresponded samples as the input, using "sorted_readgps_rmdup.bam"
        self.tornado_strategy = ''
        # tornado_result is a list of the final output in snakemake, they are the keys in tornado_strategy.yaml used for ngsplot
        self.tornado_result = ''
        self.tornado_sample_bam = ''
        self.tornado_label = ''
        self.tornado_meta_dict = defaultdict()
        # a dictionary for writing the tornado config file
        self.tornado_config_dict = defaultdict()
        # variables for computeMatrix section
        # original dict from yaml
        self.computeMatrix_dict = defaultdict()
        # after parsing, it's the final dict be used in Snakefile
        self.computeMatrix_results_dict = defaultdict()
        # used for ngsplot to make -R -L arguments flexible 
        self.tornado_ngsplot_dict = defaultdict()

    def get_tornado_strategy(self):
        yaml_handle = YamlHandle(TornadoStrategyMeta.TORNADO_STRATEGY_META_YAML_PATH.value).read_yaml()
        for i in yaml_handle:
            self.tornado_strategy = i[self.dataset_name]
        self.tornado_result = [k for k in self.tornado_strategy.keys()]
        tp = TornadoParser(self.dataset_name)
        tp.get_tornado_sample_dict()
        self.tornado_meta_dict = tp.tornado_sample_s_dict
        self.tornado_ngsplot_dict = self.tornado_strategy
        for k,v in self.tornado_ngsplot_dict.items():
            for v_k,v_v in v.items():
                if v_k == 'point':
                    v[v_k] = v_v.split('|')[0]
                if v_k == 'type':
                    v[v_k] = v_v.lower()

    def get_tornado_config(self):
        for k,v in self.tornado_strategy.items():
            self.tornado_config_dict[k] = [[self.tornado_meta_dict[s][0], self.tornado_meta_dict[s][1] + " " + v['type']] for s in v['sample'].split('|')]
    # taking a string like "S1|S2" and split S1 and S2 to a full sample name array
    def get_samples(self, sample_str):
        samples = sample_str.split('|')
        tp = TornadoParser(self.dataset_name)
        tp.get_tornado_sample_dict()
        samples_arr = [tp.tornado_sample_s_dict[s][0] for s in samples]
        return(samples_arr)
    
    # turning '1000|4000|1000' to --upstream 1000 --regionBodyLength 4000 --downstream 1000
    def region_param(self, param):
        params = param.split('|')
        if len(params) == 3:
            return('--upstream ' + params[0] + ' --regionBodyLength ' + params[1] + ' --downstream ' + params[2]) 
        else:
            return('error')
        
    # turning '1000|1000' to --upstream 1000 --downstream 1000
    def point_param(self, param):
        params = param.split('|')
        if len(params) == 2:
            return('--upstream ' + params[0] + ' --downstream ' + params[1]) 
        else:
            return('error')
    
    def get_labels(self, sample_str):
        samples = sample_str.split('|')
        tp = TornadoParser(self.dataset_name)
        tp.get_tornado_sample_dict()
        labels_list = [tp.tornado_sample_s_dict[s][1] for s in samples]
        labels_list = ['_'.join(l.split(' ')) for l in labels_list]
        labels_str = ' '.join(labels_list)
        return(labels_str)
 

    def get_computeMatrix_results(self, ref_genome):
        yaml_handle = YamlHandle(TornadoStrategyMeta.TORNADO_STRATEGY_META_YAML_PATH.value).read_yaml()
        for i in yaml_handle:
            self.computeMatrix_dict = i[self.dataset_name]
        gap = GenomeAnnotationsParser(ref_genome)
        for k,v in self.computeMatrix_dict.items():
            if v['type'].lower() == 'genebody':
                self.computeMatrix_results_dict[k] = {'sample': self.get_samples(v['sample']), 'interval_ref': gap.get_gene_spans(v['interval_ref']),\
                                                      'blacklist': gap.get_blacklist('blacklist'), 'region': self.region_param(v['region']),\
                                                      'point': self.point_param(v['point']), 'label': self.get_labels(v['sample'])}
    
            if v['type'].lower() == 'tss':
                self.computeMatrix_results_dict[k] = {'sample': self.get_samples(v['sample']), 'interval_ref': gap.get_tss_spans(v['interval_ref']),\
                                                      'blacklist': gap.get_blacklist('blacklist'), 'region': self.region_param(v['region']),\
                                                      'point': self.point_param(v['point']), 'label': self.get_labels(v['sample'])}
        
        

def main():
    tornado_test = TornadoStrategyParser('ChIPseq_211026')
    tornado_test.get_tornado_strategy()
    tornado_test.get_computeMatrix_results('mm10')
    print(tornado_test.computeMatrix_results_dict)
    print(tornado_test.tornado_result)
    print(tornado_test.tornado_strategy)
    print(tornado_test.tornado_ngsplot_dict)


if __name__ == "__main__":
    main()
            
        
