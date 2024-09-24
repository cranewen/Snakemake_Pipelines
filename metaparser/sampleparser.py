import sys, os, re
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from metaparser.yamlhandle import YamlHandle
from metaparser.metaconf import SampleMeta
from metaparser.metaconf import ChipseqInputSampleMeta
from fileio.fileutil import FileUtil
from collections import defaultdict



class SampleParser:
    def __init__(self):
        self.lane_list = ['L001_R1_001', 'L001_R2_001', 'L002_R1_001', 'L002_R2_001', 
        'L003_R1_001', 'L003_R2_001', 'L004_R1_001', 'L004_R2_001']
        self.sample_list = ''
        self.fastq_list = ''
        self.merged_fastq_list = ''
        self.fastqc_html_list = ''
        self.fastqc_zip_list = ''
        self.sample_no_input_list = ''
        self.input_list = ''
        self.input_sample_list = '' # with S#
        self.fastq_sample_list = '' # files with S# and from a directory independent to yamls
        self.fastq_sample_no_input_list = '' # only sample files with S#
        self.sample_input_dict = defaultdict() # a dictionary of sample:input_sample with S#

    # Get only sample names for dataset with no prediction of "_S" number
    def get_sample_name_list(self, dataset_name):
        yaml_handle = YamlHandle(SampleMeta.SAMPLE_YAML_PATH.value).read_yaml()
        for s in yaml_handle:
            self.sample_list = s[dataset_name]

        self.input_list = [s for s in self.sample_list if 'input' in s]

    # Get the fastq list from current directory (with S#)
    # Also create a input sample list with S#
    def get_sample_and_input_list(self, filepath):
        file_util = FileUtil(filepath)
        file_util.get_fastq_sample_list(self.sample_list)
        self.fastq_sample_list = file_util.fastq_sample_name_list # with S#
        self.input_sample_list = [f for f in self.fastq_sample_list if 'input' in f] # with S#
        self.fastq_sample_no_input_list = [f for f in self.fastq_sample_list if 'input' not in f] # with S#


    # Get sample names for dataset and append "_S" number sequentially
    def get_sample_list(self, dataset_name):
        yaml_handle = YamlHandle(SampleMeta.SAMPLE_YAML_PATH.value).read_yaml()
        for s in yaml_handle:
            self.sample_list = s[dataset_name]
        
        sample_len = len(self.sample_list)
        for i in range(sample_len):
            self.sample_list[i] = self.sample_list[i] + '_S' + str(i + 1)

    def get_sample_input_dict(self, dataset_name):
        input_yaml_handle = YamlHandle(ChipseqInputSampleMeta.CHIPSEQ_INPUT_SAMPLE_YAML_PATH.value).read_yaml()
        for i in input_yaml_handle:
            self.sample_input_dict = i[dataset_name]
        self.sample_no_input_list = [s for s in self.sample_list if 'input' not in s]


    def get_sample_dict(self, dataset_name):
        input_yaml_handle = YamlHandle(ChipseqInputSampleMeta.CHIPSEQ_INPUT_SAMPLE_YAML_PATH.value).read_yaml()
        dataset_input_dict = defaultdict()
        for i in input_yaml_handle:
            dataset_input_dict = i[dataset_name]

        for s_n in self.fastq_sample_no_input_list:
            s = '_'.join(s_n.split('_')[:-1])
            if s in dataset_input_dict:
                s_k = [v for v in self.input_sample_list if dataset_input_dict[s] in v]
                self.sample_input_dict[s_n] = s_k[0]


    def get_fastq_name_list(self, dataset_name):
        yaml_handle = YamlHandle(SampleMeta.SAMPLE_YAML_PATH.value).read_yaml()
        fastq_list = ''
        for s in yaml_handle:
            fastq_list = s[dataset_name]
        
        # A list for [S1, S2, S3, ...]
        s_n = ['S'+ str(i + 1) for i in range(len(fastq_list))]
        fastq_list = [f + '_' + s for f, s in zip(fastq_list, s_n)]
        fastq_list = [f + '_' + l + '.fastq.gz' for f in fastq_list for l in self.lane_list]
        self.fastq_list = fastq_list

    def get_merged_fastq_name_list(self, dataset_name):
        yaml_handle = YamlHandle(SampleMeta.SAMPLE_YAML_PATH.value).read_yaml()
        fastq_list = ''
        r_n = ['_R1_merged.fastq.gz', '_R2_merged.fastq.gz']
        for s in yaml_handle:
            fastq_list = s[dataset_name]

        # A list for [S1, S2, S3, ...]
        s_n = ['S'+ str(i + 1) for i in range(len(fastq_list))]
        fastq_list = [f + '_' + s for f, s in zip(fastq_list, s_n)]
        fastq_list = [f + r for f in fastq_list for r in r_n]
        self.merged_fastq_list = fastq_list

    def get_fastqc_name_list(self, dataset_name):
        yaml_handle = YamlHandle(SampleMeta.SAMPLE_YAML_PATH.value).read_yaml()
        fastq_list = ''
        fastqc_html_list = ''
        fastqc_zip_list = ''
        r_n = ['_R1', '_R2']
        for s in yaml_handle:
            fastq_list = s[dataset_name]

        # A list for [S1, S2, S3, ...]
        s_n = ['S'+ str(i + 1) for i in range(len(fastq_list))]
        fastq_list = [f + '_' + s for f, s in zip(fastq_list, s_n)]
        fastq_list = [f + r for f in fastq_list for r in r_n]
        fastqc_html_list = [f + '_merged_fastqc.html' for f in fastq_list]
        self.fastqc_html_list = fastqc_html_list
        fastqc_zip_list = [f + '_merged_fastqc.zip' for f in fastq_list]
        self.fastqc_zip_list = fastqc_zip_list
    

# Tests
def main():
    sample = SampleParser()
    sample.get_sample_list('ChIPseq_test')
    print(sample.sample_list)


if __name__ == "__main__":
    main()





