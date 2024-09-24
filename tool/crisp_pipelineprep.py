import yaml
import os, errno
from collections import defaultdict
from pathlib import Path
import shutil

class ConfigParser:
    def __init__(self, dataset, species = 'homo_sapiens', assembly = 'hg38'):
        self.dataset = dataset
        self.strategy_dict = defaultdict()
        # a list for Snakefile
        self.strategies = []
        self.meta_dict = defaultdict()
        self.attributes = ''
        self.samples = ''
        self.d0_sample = ''
        self.species = species
        self.assembly = assembly
        #self.hEpi_sgRNA_library = '/mnt/storage/dept/pedonc/Reference/CRISPR_libraries/hEpi_sgRNA_library.txt'
        self.hEpi_sgRNA_library = '/mnt/storage/dept/pedonc/Reference/CRISPR_libraries/hEpi_sgRNA_library_trimmed_1.txt'
        self.output_dir = ''

    def read_yaml(self, yaml_path):
        try:
            return(yaml.safe_load_all(open(yaml_path)))
        except FileNotFoundError:
            print(str(yaml_path) + ' can not be found')
    
    # return a sample list, strategy is the key in self.strategy_dict; e.g. VTP_DMSO 
    def get_samples_from_strategy(self, strategy):
        '''
        strategies = self.read_yaml('strategy.yaml')
        for s in strategies:
            self.strategy_dict = s[self.dataset]
        '''
        samples = defaultdict()
        # add baseline 1 to total column length
        ncol_matrix = len(self.strategy_dict[strategy][0].split(':')) + 1
        meta = self.read_yaml('meta.yaml')
    
        d0_sample = defaultdict()
        for m in meta:
            self.meta_dict = m[self.dataset]
        for k,v in self.meta_dict.items():
            temp_list = [0] * ncol_matrix 
            if v['attribute'] == 'D0':
                temp_list[0] = 1
                temp_list.append(v['path'])
                samples[k] = temp_list
                
            try:
                # get the row index of 1 in the matrix
                index1 = self.strategy_dict[strategy][0].split(':').index(v['attribute'])
                temp_list[0] = 1
                temp_list[index1 + 1] = 1
                temp_list.append(v['path'])
                samples[k] = temp_list
            except ValueError:
                #print('Attribute doesn\'t exist!')
                continue
               
        return(samples)
        
        

    def write_all_files(self):
        matrix_dict = defaultdict()
        strategy = self.read_yaml('strategy.yaml')
        for s in strategy:
            self.strategy_dict = s[self.dataset]
        for k,v in self.strategy_dict.items():
            self.strategies.append(k)
            matrix_path = 'Experiment_path/mageck-vispr/' + self.dataset + '/input/'
            matrix_file = matrix_path + self.dataset + '_design_matrix_' + k + '.txt'
            config_path = 'Experiment_path/mageck-vispr/' + self.dataset + '/mvispr_' + self.dataset + '/'
            config_strategy_path = config_path + k
            Path(matrix_path).mkdir(parents=True, exist_ok=True)
            Path(config_path).mkdir(parents=True, exist_ok=True)
            Path(config_strategy_path).mkdir(parents=True, exist_ok=True)
            shutil.copy('Snakefile', config_strategy_path + '/Snakefile')
            samples = self.get_samples_from_strategy(k)
            with open(matrix_file, 'w+') as f:
                head_line = self.strategy_dict[k][0].split(':')
                f.write('Sample\t baseline\t' + '\t'.join(head_line) + '\n')
                for k,v in samples.items():
                    f.write(k + '\t')
                    f.write('\t'.join(str(x) for x in v[:-1]))
                    f.write('\n')

            with open(config_strategy_path + '/config.yaml', 'w+') as f:
                f.write('library: ' + self.hEpi_sgRNA_library + '\n')
                f.write('species: ' + self.species + '\n')
                f.write('assembly: ' + self.assembly + '\n')
                f.write('targets:\n    genes: true\n')
                f.write('sgrnas:\n    update-efficiency: false\n    trim-5: AUTO\n    len: AUTO\n    annotate-sgrna: false\n    annotate-sgrna-efficiency: false\n')
                f.write('samples:\n')
                for k,v in samples.items():
                    f.write('    ' + k + ':\n')
                    f.write('         - ' + v[-1] + '\n')
                f.write('countpair: false\n')
                f.write('threads: 4\n')
                f.write('correct_cnv: false\n')
                f.write('cnv_norm: /dev/null\n')
                f.write('cnv_cell_line: HL60\n')
                f.write('experiments:\n    \"mle\":\n        designmatrix: ' + matrix_file + '\n')
                

def main():
    cp = ConfigParser('sgRNAseq_test')
    cp.write_all_files()

if __name__ == '__main__':
    main()


