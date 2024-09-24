import glob, sys, os, yaml
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from metaparser.dirbuilder import DirBuilder
from metaparser.tornadostrategyparser import TornadoStrategyParser
from metaparser.genomeannotationsparser import GenomeAnnotationsParser
#from metaparser.yamlhandle import YamlHandle
#from metaparser.metaconf import *
from collections import defaultdict
import subprocess
import argparse

parser = argparse.ArgumentParser(description='2 arguments are needed. dataset_name & ref_genome')
parser.add_argument('-d', '--dataset', type=str, required=True)
parser.add_argument('-f', '--ref_genome', type=str)
# we currently only have "protcoding_genes" as the input
parser.add_argument('-p', '--gene_list', type=str)
args = parser.parse_args()

class TornadoWriter:
    def __init__(self, dataset_name, ref_name):
        self.dataset_name = dataset_name
        self.ref_name = ref_name
        self.tornado_bam_sorted_path = ''
        self.tornado_bam_rmdup_path = ''
        self.tornado_config_dict = ''
        #self.tornado_config_path = ''
        self.tornado_dataset_path = ''
        self.tornado_sorted_dataset_path = ''
        self.tornado_rmdup_dataset_path = ''
        self.tornado_sorted_config_path = ''
        self.tornado_rmdup_config_path = ''
        self.gene_lists_txt_path = ''
        
    def get_tornado_config_dict(self):
        tsp = TornadoStrategyParser(self.dataset_name)
        tsp.get_tornado_strategy()
        tsp.get_tornado_config()
        self.tornado_config_dict = tsp.tornado_config_dict

    def get_tornado_path(self):
        dirs = DirBuilder(self.dataset_name)
        dirs.build_tornado_dirs(self.ref_name)
        #self.tornado_config_path = dirs.tornado_config_dir
        self.tornado_dataset_path = dirs.tornado_dataset_dir
        self.tornado_sorted_dataset_path = dirs.tornado_sorted_dataset_dir
        self.tornado_rmdup_dataset_path = dirs.tornado_rmdup_dataset_dir
        self.tornado_bam_sorted_path = dirs.tornado_bam_sorted_dir
        self.tornado_bam_rmdup_path = dirs.tornado_bam_rmdup_dir
        self.tornado_sorted_config_path = dirs.tornado_sorted_config_dir
        self.tornado_rmdup_config_path = dirs.tornado_rmdup_config_dir
    
    # gene_list_type, e.g. 'protcoding_genes')
    def get_gene_lists_path(self, gene_list_type):
        gsp = GenomeAnnotationsParser(self.ref_name)
        self.gene_lists_txt_path = gsp.get_gene_list(gene_list_type)

    def write_config_txt_file(self):
        #os.makedirs(self.tornado_config_path, exist_ok = True)
        os.makedirs(self.tornado_sorted_config_path, exist_ok = True)
        os.makedirs(self.tornado_rmdup_config_path, exist_ok = True)
        for k,v in self.tornado_config_dict.items():
            with open(self.tornado_sorted_config_path + k + '_config.txt', 'w') as config_file:
                for s in v:
                    config_file.write(self.tornado_bam_sorted_path + s[0] + '_' + self.ref_name + '_sorted_readgps.bam' + ' ' + self.gene_lists_txt_path + ' \"' + s[1] + '\"\n')

            with open(self.tornado_rmdup_config_path + k + '_config.txt', 'w') as config_file:
                for s in v:
                    config_file.write(self.tornado_bam_rmdup_path + s[0] + '_' + self.ref_name + '_sorted_readgps_rmdup.bam' + ' ' + self.gene_lists_txt_path + ' \"' + s[1] + '\"\n')
             


def main():
    tw = TornadoWriter(args.dataset, args.ref_genome)
    tw.get_tornado_config_dict()
    tw.get_tornado_path()
    #tw.get_gene_lists_path('protcoding_genes')
    tw.get_gene_lists_path(args.gene_list)
    
    #tw.write_config_txt_file()

    

if __name__ == "__main__":
    main()
