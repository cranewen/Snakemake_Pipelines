import glob, sys, os, yaml
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from metaparser.dirbuilder import DirBuilder
from metaparser.genomeannotationsparser import GenomeAnnotationsParser
from .permutationtool import PermutationTool
#from metaparser.yamlhandle import YamlHandle
#from metaparser.metaconf import *
from collections import defaultdict
import subprocess
import argparse

class IntersectWriter:
    # gene_spans_tss_ref are the gene_spans files, ref_genome is "hg38, etc"
    def __init__(self, dataset_name, ref_genome):
        self.dataset_name = dataset_name
        self.ref_genome = ref_genome
        self.sample_list = '' 

        self.enhancers_list = ''
        self.gene_spans_list = ''
        self.tss_spans_list = ''
        self.blacklist_list = ''

        # gap is short for genome annotations parser
        self.gap = GenomeAnnotationsParser(ref_genome)
        self.perm_tool = PermutationTool(dataset_name, ref_genome)
        
    def __repr__(self):
        return "Display-- dataset:%s ref_genome:%s enhancers:%s " \
                % (self.dataset_name, self.ref_genome, self.enhancers_list)



    # self.enhancers_list only contains the "raw" inputed params split into a list, it doesn't include the entire file path
    def set_enhancers(self, enhancers):
        enhancers_list = enhancers.split('|')
        self.enhancers_list = enhancers_list

    def get_enhancers(self):
        self.enhancers_list = [self.gap.get_enhancers(p) for p in self.enhancers_list]
        return self.enhancers_list

    # set the gene_spans list by the keys
    def set_gene_spans(self, gene_spans):
        gene_spans_list = gene_spans.split('|')
        self.gene_spans_list = gene_spans_list

    # get the gene_spans list of whole paths 
    def get_gene_spans(self):
        self.gene_spans_list = [self.gap.get_gene_spans(g) for g in self.gene_spans_list]
        return self.gene_spans_list


    # set the tss_spans list by the keys
    def set_tss_spans(self, tss_spans):
        tss_spans_list = tss_spans.split('|')
        self.tss_spans_list = tss_spans_list

    # get the tss_spans list of whole paths 
    def get_tss_spans(self):
        self.tss_spans_list = [self.gap.get_tss_spans(g) for g in self.tss_spans_list]
        return self.tss_spans_list

    # blacklist only has one file, so we pass the default parameter as 'blacklist'
    def set_blacklist(self, blacklist):
        blacklist = blacklist.split('|')
        self.blacklist_list = blacklist

    def get_blacklist(self):
        self.blacklist = [self.gap.get_blacklist(b) for b in self.blacklist_list]
        return self.blacklist

        

    # a is the list of samples or samples+enhancers, b is the list of refs
    def write_intersect_sh(self, filepath, intersect_tool, samples, enhancers, gene_spans, tss_spans, blacklist, chrNameLength, input_dir, output_dir):
        # samples list only includes the samples, not the reference files
        os.makedirs(filepath, exist_ok = True)
        #samples_list = []
        samples_list = [input_dir + s for s in samples]
        samples_perm_list = self.perm_tool.perm2_samples(samples_list)

        enhancers_list = []
        gene_spans_list = []
        tss_spans_list = []
        blacklist_list = [] 

        self.set_enhancers(enhancers)
        self.set_gene_spans(gene_spans)
        self.set_tss_spans(tss_spans)
        self.set_blacklist(blacklist)

        # setting all the reference to sample pairs
        if self.enhancers_list[0] != '':
            enhancers_list = self.get_enhancers()
        enhancers_samples_list = self.perm_tool.ref_compare(samples_list, enhancers_list)

        if self.gene_spans_list[0] != '':
            gene_spans_list = self.get_gene_spans()
        gene_spans_samples_list = self.perm_tool.ref_compare(samples_list, gene_spans_list)

        if self.tss_spans_list[0] != '':
            tss_spans_list = self.get_tss_spans()
        tss_spans_samples_list = self.perm_tool.ref_compare(samples_list, tss_spans_list)

        if self.blacklist_list[0] != '':
            blacklist_list = self.get_blacklist()
        blacklist_samples_list = self.perm_tool.ref_compare(samples_list, blacklist_list)
        
           
        
        with open(filepath + '/bedtool_intersect.sh', 'w') as f:
            f.write('#!/bin/bash\n')

            for s in samples_perm_list: 
                f.write(intersect_tool + ' intersect -g ' + chrNameLength + ' -u -a ' + s[0] + '_' + self.ref_genome + '_peaks.narrowPeak' + ' -b ' + s[1] + '_' + self.ref_genome + '_peaks.narrowPeak' + ' > ' + output_dir + '/' + s[0].split('/')[-1] + '_vs_' + s[1].split('/')[-1] + '_' + self.ref_genome + '_nofloor.narrowPeak_overlaps.txt' + '\n') 
       
            for es in enhancers_samples_list:
                enhancer_name = es[1].split('/')[-1].split('.')[0]
                f.write(intersect_tool + ' intersect -g ' + chrNameLength + ' -u -a ' + es[0] + '_' + self.ref_genome + '_peaks.narrowPeak' + ' -b ' + es[1] + ' > ' + output_dir + '/' + es[0].split('/')[-1] + '_' + self.ref_genome + '_' + enhancer_name + '_nofloor.narrowPeak_gene_overlaps.txt' + '\n') 
                
            for gs in gene_spans_samples_list:
                f.write(intersect_tool + ' intersect -g ' + chrNameLength + ' -u -a ' + gs[0] + '_' + self.ref_genome + '_peaks.narrowPeak' + ' -b ' + gs[1] + ' > ' + output_dir + '/' + gs[0].split('/')[-1] + '_' + self.ref_genome + '_nofloor.narrowPeak_gene_overlaps.txt' + '\n') 
                f.write(intersect_tool + ' intersect -g ' + chrNameLength + ' -wo -a ' + gs[0] + '_' + self.ref_genome + '_peaks.narrowPeak' + ' -b ' + gs[1] + ' > ' + output_dir + '/' + gs[0].split('/')[-1] + '_' + self.ref_genome + '_nofloor.narrowPeak_gene_identity.txt' + '\n') 

            for ts in tss_spans_samples_list:
                f.write(intersect_tool + ' intersect -g ' + chrNameLength + ' -u -a ' + ts[0] + '_' + self.ref_genome + '_peaks.narrowPeak' + ' -b ' + ts[1] + ' > ' + output_dir + '/' + ts[0].split('/')[-1] + '_' + self.ref_genome + '_nofloor.narrowPeak_gene_TSS_overlaps.txt' + '\n') 
                f.write(intersect_tool + ' intersect -g ' + chrNameLength + ' -wo -a ' + ts[0] + '_' + self.ref_genome + '_peaks.narrowPeak' + ' -b ' + ts[1] + ' > ' + output_dir + '/' + ts[0].split('/')[-1] + '_' + self.ref_genome + '_nofloor.narrowPeak_gene_TSS_identity.txt' + '\n') 

            for bs in blacklist_samples_list:
                blacklist_name = bs[1].split('/')[-1].split('.')[0]
                f.write(intersect_tool + ' intersect -g ' + chrNameLength + ' -u -a ' + bs[0] + '_' + self.ref_genome + '_peaks.narrowPeak' + ' -b ' + bs[1] + ' > ' + output_dir + '/' + bs[0].split('/')[-1] + '_' + blacklist_name + '_nofloor.narrowPeak_gene_overlaps.txt' + '\n') 


        
def main():
    iw = IntersectWriter('ChIPseq_211106', 'mm10')
    print(iw.sample_list)
    


if __name__ == "__main__":
    main()
