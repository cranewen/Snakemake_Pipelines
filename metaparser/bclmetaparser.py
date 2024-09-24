import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import yaml
import datetime
import os, errno
from metaparser.yamlhandle import YamlHandle
from metaparser.metaconf import BclMeta, DirConfig
from collections import defaultdict


class BclMetaParser:
    def __init__(self):
        self.dataset_path = '' # project bcl2fastq demultiplex path
        self.demultiplex_path = '' # original demultiplex data path
        self.bcl2fastq_path = '' # a subdirectory of demultiplex path, which stores dataset output
        self.selected_dataset = '' # dataset data in bcl_meta.yaml
        
        self.raw_bcl_path = ''
        # not the full path, it's just a folder name like "220622_NB552101_0365_AHK2LMBGXL"
        self.raw_bcl_path_name = ''
        self.sample_sheet_path = ''
        self.bcl_output_dir = ''
        self.bcl_meta = YamlHandle(BclMeta.BCL_YAML_PATH.value).read_yaml()
        self.bcl_config = YamlHandle(DirConfig.DIR_CONFIG_PATH.value).read_yaml()

       
    def set_dataset_path(self, dataset_name):
        for c in self.bcl_config:
            self.demultiplex_path = c[DirConfig.DATASET_BCL_DEMULTIPLEX.value]
            self.bcl2fastq_path = self.demultiplex_path + c[DirConfig.DATASET_BCL2FASTQ.value]
            self.bcl_output_dir = self.bcl2fastq_path + dataset_name + '/'

        for m in self.bcl_meta:
            self.selected_dataset = m[dataset_name]

        self.raw_bcl_path = self.demultiplex_path + self.selected_dataset[BclMeta.RAW_BCL_PATH.value]
        self.sample_sheet_path = self.demultiplex_path + self.selected_dataset[BclMeta.SAMPLE_SHEET_NAME.value] + '/'
        self.raw_bcl_path_name = self.selected_dataset[BclMeta.RAW_BCL_PATH.value]

    def get_all_info(self):
        all_info = defaultdict()
        count = 1
        for b in self.bcl_meta:
            for k, v in b.items():
                all_info[count] = [k, v['run_date']]
                count += 1

        return(all_info)


def main():
    meta = BclMetaParser()
    meta.set_dataset_path('test_dataset')
    print(meta.selected_dataset)
    temp = meta.get_all_info()
    print(temp)

if __name__ == "__main__":
    main()
