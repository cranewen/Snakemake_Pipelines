import sys, os, glob
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from metaparser.metaconf import ChipseqInputSampleIndexMeta
from metaparser.sampleparser import SampleParser
from metaparser.yamlhandle import YamlHandle
from fileio.fileutil import FileUtil
from collections import defaultdict

class SampleIndexWriter:
    def __init__(self, dataset_name):
        self.dataset_name = dataset_name
        self.realtime_sample_list = '' # from IO 
        self.sample_list = '' # from original sample.yaml
        self.sample_index_dict = defaultdict() # a dictionary from chipseq_input_sample_index.yaml
        self.realtime_sample_dict = defaultdict() # a dictionary stores the final data for MACS2

    # if the bcl_fastq pipeline's output files have been moved or deleted, we are looking up the file system for real samples
    def get_realtime_sample(self, filepath):
        sample_parser = SampleParser()
        sample_parser.get_sample_list(self.dataset_name)
        self.sample_list = sample_parser.sample_list

        file_util = FileUtil(filepath)
        file_util.get_fastq_sample_list(self.sample_list)
        self.realtime_sample_list = file_util.fastq_sample_name_list

        index_yaml_handle = YamlHandle(ChipseqInputSampleIndexMeta.CHIPSEQ_INPUT_SAMPLE_INDEX_YAML_PATH.value).read_yaml()
        for i in index_yaml_handle:
            self.sample_index_dict = i[self.dataset_name]

    # the bcl_fastq pipeline's output files hasn't been touched, we use as it is
    def get_original_sample(self):
        sample_parser = SampleParser()
        sample_parser.get_sample_list(self.dataset_name)
        self.sample_list = sample_parser.sample_list
        self.realtime_sample_list = self.sample_list

        index_yaml_handle = YamlHandle(ChipseqInputSampleIndexMeta.CHIPSEQ_INPUT_SAMPLE_INDEX_YAML_PATH.value).read_yaml()
        for i in index_yaml_handle:
            self.sample_index_dict = i[self.dataset_name]


    def write_chipseq_input_sample_yaml(self):
        for k, v in self.sample_index_dict.items():
            for r_s in v:
                for s in self.sample_list:
                    if r_s in s:
                        self.realtime_sample_dict[s] = k
        chipseq_input_sample_yaml = '/mnt/storage/dept/pedonc/src/snakemake_pipelines/metaconfig/chipseq_input_sample.yaml'
        
        with open(chipseq_input_sample_yaml, 'a') as f:
            if os.stat(chipseq_input_sample_yaml).st_size != 0:
                f.write('\n')
            f.write(self.dataset_name + '\n')
            for k, v in self.realtime_sample_dict.items():
                f.write(' ' + k + ':\n')
                f.write('  \'' + v + '\'\n')



#Tests
def main():
    test_obj = SampleIndexWriter('ChIPseq_200807_test')
    test_obj.get_original_sample()
    test_obj.write_chipseq_input_sample_yaml()
    print(test_obj.realtime_sample_dict)



if __name__ == "__main__":
    main()
