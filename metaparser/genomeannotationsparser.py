import sys, os, re
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from metaparser.yamlhandle import YamlHandle
from metaparser.metaconf import GenomeAnnotationsMeta, DirConfig
from collections import defaultdict

class GenomeAnnotationsParser:
    def __init__(self, ref_genome):
        self.ref_genome = ref_genome
        # a dictionary with "hg38, mm10, etc" as the keys, stores the information from gene_spans_meta.yaml
        self.genome_dict = ''
        # gene_spans and tss_spans are sharing the same directory, which is the gene_spans_root
        self.gene_lists_root_path = ''
        self.gene_spans_root_path = ''
        self.gene_spans_file_path = ''
        self.tss_spans_file_path = ''
        self.blacklist_file_path = ''
        self.enhancers_file_path = ''
        config_path_handle = YamlHandle(DirConfig.DIR_CONFIG_PATH.value).read_yaml()
        # define gene_spans root file path, later condat the string together to get the final file path
        for i in config_path_handle:
            self.gene_spans_root_path = i[DirConfig.GENE_SPANS_PATH.value]
            self.gene_lists_root_path = i[DirConfig.GENE_LISTS_PATH.value]
        

    def get_gene_spans(self, param):
        yaml_handle = YamlHandle(GenomeAnnotationsMeta.YAML_PATH.value).read_yaml()

        for i in yaml_handle:
            self.genome_dict = i[self.ref_genome]

        for k,v in self.genome_dict[GenomeAnnotationsMeta.GENE_SPANS.value].items():
            if k == param:
                self.gene_spans_file_path = self.gene_spans_root_path + v[GenomeAnnotationsMeta.DIR.value]
        return(self.gene_spans_file_path)
        
    def get_tss_spans(self, param):
        yaml_handle = YamlHandle(GenomeAnnotationsMeta.YAML_PATH.value).read_yaml()

        for i in yaml_handle:
            self.genome_dict = i[self.ref_genome]

        for k,v in self.genome_dict[GenomeAnnotationsMeta.TSS_SPANS.value].items():
            if k == param:
                self.tss_spans_file_path = self.gene_spans_root_path + v[GenomeAnnotationsMeta.DIR.value]
        return(self.tss_spans_file_path)
 


    def get_blacklist(self, param):
        yaml_handle = YamlHandle(GenomeAnnotationsMeta.YAML_PATH.value).read_yaml()

        for i in yaml_handle:
            self.genome_dict = i[self.ref_genome]

        for k,v in self.genome_dict[GenomeAnnotationsMeta.BLACKLIST.value].items():
            if k == param:
                self.blacklist_file_path = self.gene_spans_root_path + v[GenomeAnnotationsMeta.DIR.value]

        return(self.blacklist_file_path)
        
    def get_gene_list(self, param):
        yaml_handle = YamlHandle(GenomeAnnotationsMeta.YAML_PATH.value).read_yaml()

        for i in yaml_handle:
            self.genome_dict = i[self.ref_genome]
    
        for k,v in self.genome_dict[GenomeAnnotationsMeta.GENE_LISTS.value].items():
            if k == param:
                return(self.gene_lists_root_path + v[GenomeAnnotationsMeta.DIR.value])

    # Enhancer params are GRCh38, MOLM13, MV411
    def get_enhancers(self, param):
        if param == '':
            return('')
        yaml_handle = YamlHandle(GenomeAnnotationsMeta.YAML_PATH.value).read_yaml()

        for i in yaml_handle:
            self.genome_dict = i[self.ref_genome]
    
        # Enhancers root directory is the same as gene_spans
        for k,v in self.genome_dict[GenomeAnnotationsMeta.ENHANCERS.value].items():
            if k == param:
                return(self.gene_spans_root_path + v[GenomeAnnotationsMeta.DIR.value])


def main():
    gsp = GenomeAnnotationsParser('hg38')
    print(gsp.get_tss_spans('') == '')



if __name__ == "__main__":
    main()
