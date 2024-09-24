import sys, os, re
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from metaparser.yamlhandle import YamlHandle
from metaparser.dirbuilder import DirBuilder
from metaparser.metaconf import CRISPRStrategyMeta, CRISPRMeta
from collections import defaultdict
from pathlib import Path
import shutil

class CRISPRStrategyParser():
    def __init__(self, dataset_name, ref_genome):
        self.dataset_name = dataset_name
        self.ref_genome = ref_genome
        self.crispr_meta_dict = defaultdict()
        self.strategy_dict = defaultdict() 
        self.matrix_names = []
        self.matrix_dict = defaultdict()
        self.sample_dict = defaultdict()
        self.hEpi_sgRNA_library =  '/mnt/storage/dept/pedonc/Reference/CRISPR_libraries/hEpi_sgRNA_library_trimmed_1.txt'
        self.species = {'hg38': 'homo_sapiens', 'hg19': 'homo_sapiens'}

    def get_crispr_strategy(self):
        yaml_handle = YamlHandle(CRISPRStrategyMeta.CRISPR_STRATEGY_META_YAML_PATH.value).read_yaml()
        for i in yaml_handle:
            self.strategy_dict = i[self.dataset_name]

    def get_crispr_meta(self):
        yaml_handle = YamlHandle(CRISPRMeta.CRISPR_YAML_PATH.value).read_yaml()
        for i in yaml_handle:
            self.crispr_meta_dict = i[self.dataset_name]

    def get_samples_from_attribute(self, attr):
        # a dictionary for catagorizing samples. 
        # e.g. DMSO and VTP, so DMSO and VTP will be the keys and samples name will be the values
        for k,v in self.crispr_meta_dict.items():
            try:
                self.sample_dict[v[attr]].append(k)
            except KeyError:
                self.sample_dict[v[attr]] = [k]
                
        
    def get_samples_from_attribute_single_strategy(self, attr):
        sample_dict = defaultdict()
        for k,v in self.crispr_meta_dict.items():
            try:
                sample_dict[v[attr]].append(k)
            except KeyError:
                sample_dict[v[attr]] = [k]
        return sample_dict
            

    # writing all the matrix files and config and Snakefiles 
    def crispr_prep(self):
        dirs = DirBuilder(self.dataset_name)
        dirs.build_sgrnaseq_dirs(self.ref_genome)
        mvispr_output_dir = dirs.mvispr_output_dir
        attribute = ''
        sample_dict = ''
        samples_for_config_yaml = []
        for k,v in self.strategy_dict.items():
            attribute = v.split('=')[0]
            sample_dict = self.get_samples_from_attribute_single_strategy(attribute)
            sample_group = v.split('=')[1].split(':')
            file_name = self.dataset_name + '_design_matrix_' + k + '.txt'
            # test file path
            matrix_path = 'crispr_input_test/' + file_name
            with open(matrix_path, 'w') as f:
                f.write('Sample    baseline    ')
                f.write(sample_group[0] + '\t' + sample_group[1] + '\n')
                f.write(sample_dict['none'][0] + '\t1\t' + str(0) + '\t0\n')
                for sa in sample_group:
                    for s in sample_dict[sa]:
                        samples_for_config_yaml.append(s)
                        if sa == sample_group[0]:
                            f.write(s + '\t1\t' + str(1) + '\t0\n')
                        else:
                            f.write(s + '\t1\t' + str(0) + '\t1\n')
                            
            
            config_filepath = mvispr_output_dir + k
            Path(config_filepath).mkdir(parents=True, exist_ok=True)
            with open(config_filepath + '/config.yaml', 'w+') as f:
                print('library: ' + self.hEpi_sgRNA_library + '\n')
                print('species: ' + self.species[self.ref_genome] + '\n')
                print('assembly: ' + self.ref_genome + '\n')
                print('targets:\n    genes: true\n')
                print('sgrnas:\n    update-efficiency: false\n    trim-5: AUTO\n    len: AUTO\n    annotate-sgrna: false\n    annotate-sgrna-efficiency: false\n')
                print('samples:\n')
                for s in samples_for_config_yaml:
                    print('    ' + s + ':\n')
                    print('         - ' + self.crispr_meta_dict[s]['path'] + '\n')
                print('countpair: false\n')
                print('threads: 4\n')
                print('correct_cnv: false\n')
                print('cnv_norm: /dev/null\n')
                print('cnv_cell_line: HL60\n')
                print('experiments:\n    \"mle\":\n        designmatrix: ' + matrix_path + '\n')
                
            samples_for_config_yaml = []
            

        

def main():
    crisprtest = CRISPRStrategyParser('sgRNAseq_test', 'hg38')
    crisprtest.get_crispr_strategy()
    crisprtest.get_crispr_meta()
    #crisprtest.get_samples_from_attribute('treatment')
    crisprtest.crispr_prep()
    strategy = crisprtest.strategy_dict
    meta = crisprtest.crispr_meta_dict



if __name__ == "__main__":
    main()
