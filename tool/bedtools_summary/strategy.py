# read strategy from yaml
import sys, os
sys.path.insert(0, '/home/yw900/lab_pipelines')
import yaml
from metaparser.yamlhandle import YamlHandle
from metaparser.metaconf import BedtoolSummaryConfig
from pprint import pprint

class Strategy:
    def __init__(self, dataset: str):
        self.bedtool_summary_config = YamlHandle(BedtoolSummaryConfig.BEDTOOL_SUMMARY_CONFIG_YAML_PATH.value).read_yaml()
        self.strategy = None
        try:
            for b in self.bedtool_summary_config:
                self.strategy = b[dataset]
        except KeyError:
            print('check config file if the dataset is there')

    def __repr__(self):
        print(f'===== Strategy object members =====')
        pprint(self.__dict__)
        return f'===== Strategy object members ====='

def main():
    s = Strategy('ChIPseq_221118_JR')
    print(s)

if __name__ == '__main__':
    main()
    

