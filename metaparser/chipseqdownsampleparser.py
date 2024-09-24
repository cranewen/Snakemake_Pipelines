import sys, os, re
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from metaparser.yamlhandle import YamlHandle
from metaparser.metaconf import ChipseqDownsampleMeta
from collections import defaultdict


class ChipseqDownsampleParser:
    def __init__(self, dataset_name):
        self.dataset_name = dataset_name
        self.chipseq_downsample_meta = defaultdict()
        self.chipseq_downsample_method_meta = defaultdict()
        # making a new dictionary with S# as the keys

    def get_chipseq_downsample_meta(self):
        yaml_handle = YamlHandle(ChipseqDownsampleMeta.CHIPSEQ_DOWNSAMPLE_META_YAML_PATH.value).read_yaml()
        for i in yaml_handle:
            self.chipseq_downsample_meta = i[self.dataset_name]
            
        return self.chipseq_downsample_meta

    def get_meta_by_method(self, ds_method):
        for k,v in self.chipseq_downsample_meta.items():
            self.chipseq_downsample_method_meta[k] = defaultdict()
            if v['down_million_' + ds_method] != 'None':
                self.chipseq_downsample_method_meta[k]['name'] = 'ds_' + ds_method + '_' + v['down_million_' + ds_method]
                self.chipseq_downsample_method_meta[k]['frac'] = v['frac_' + ds_method]
            else:
                self.chipseq_downsample_method_meta[k]['name'] = 'ds_' + ds_method + '_' + str(round(v['rmdup_aligned_reads'] / 1000000, 1)) + 'M'
                self.chipseq_downsample_method_meta[k]['frac'] = 0.99
            
        return self.chipseq_downsample_method_meta
        



def main():
    print('some tests')

if __name__ == "__main__":
    main()
