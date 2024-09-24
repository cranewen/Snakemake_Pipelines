import sys, os, re
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from metaparser.yamlhandle import YamlHandle
from metaparser.metaconf import TornadoStrategyMeta, TornadoMeta
from collections import defaultdict


class TornadoParser:
    def __init__(self, dataset_name):
        self.dataset_name = dataset_name
        self.tornado_sample_dict = defaultdict()
        # making a new dictionary with S# as the keys
        self.tornado_sample_s_dict = defaultdict()

    def get_tornado_sample_dict(self):
        yaml_handle = YamlHandle(TornadoMeta.TORNADO_META_YAML_PATH.value).read_yaml()
        for i in yaml_handle:
            self.tornado_sample_dict = i[self.dataset_name]
        for k,v in self.tornado_sample_dict.items():
            self.tornado_sample_s_dict[k.split('_')[-1]] = [k, v['label']]
            



def main():
    tp = TornadoParser('ChIPseq_210322_HU_test')
    tp.get_tornado_sample_dict()
    print(tp.tornado_sample_s_dict)

if __name__ == "__main__":
    main()
