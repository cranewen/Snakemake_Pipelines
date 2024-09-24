# a parser for flagstats from ChIPseq 
import sys, os
sys.path.insert(0, '/home/yw900/lab_pipelines')
from collections import defaultdict
from typing import Dict
import pandas as pd
from pprint import pprint
from db.dbconnection import DBconnector
from db.cpctdb import DatasetSeqRunSample, DatasetSeqRun, Sample, SampleGroup, SampleAttribute, Protocol
from metaparser.sampleparser import SampleParser
from icecream import ic

class Stats:
    def __init__(self, dataset, tss_all_samples_path: str, flyreads: bool):
        self.dataset = dataset
        self.tss_all_samples_path = tss_all_samples_path
        self.stat_dict = defaultdict()
        self.flyreads = flyreads # passing from Dataset class

    def __repr__(self):
        print(f'===== Stats object members =====')
        pprint(self.__dict__)
        return f'===== Stats object members ====='
        
    # flagstat_path['sorted'], flagstat_path['rmdup'], etc.
    def parse_flagstat(self, flagstat_path: Dict):
        #sample_names = self.raw_counts_df.columns.values[5:]
        s = SampleParser()
        s.get_sample_list(self.dataset)
        sample_names = s.sample_list
        print(sample_names)

        self.stat_dict['sorted'] = defaultdict()
        self.stat_dict['rmdup'] = defaultdict()
        # for dm6 there is only 1 element in the array, because we only need the spike-in properly paired data
        # so it will be {sample_name:['properly paired']}
        self.stat_dict['sorted_dm6'] = defaultdict()
        self.stat_dict['rmdup_dm6'] = defaultdict()

    def parse_zscore(self):
        self.raw_counts_df = pd.read_csv(self.tss_all_samples_path, sep=',')

class FlagStats:
    def __init__(self, dataset, flagstat_path, ref, ref_sp=None):
        self.dataset = dataset
        self.flagstat_path = flagstat_path
        self.stat_dict = defaultdict()
        self.ref = ref
        self.ref_sp = ref_sp # passing from Dataset class

    def parse_flagstat(self):
        #sample_names = self.raw_counts_df.columns.values[5:]
        s = SampleParser()
        s.get_sample_list(self.dataset)
        sample_names = s.sample_list
        print(sample_names)

        '''
        self.stat_dict['sorted'] = defaultdict()
        self.stat_dict['rmdup'] = defaultdict()
        # for dm6 there is only 1 element in the array, because we only need the spike-in properly paired data
        # so it will be {sample_name:['properly paired']}
        self.stat_dict['sorted_spikein'] = defaultdict()
        self.stat_dict['rmdup_spikein'] = defaultdict()
        '''

        for s in sample_names:
            self.stat_dict[s] = defaultdict()
            if self.ref_sp:
                try:
                    with open(f'{self.flagstat_path}BAM_sorted/{self.dataset}_{self.ref_sp}/samtools_flagstat/{s}_{self.ref_sp}_sorted_readgps_samtools_flagstat.txt') as f:
                        for line in f:
                            if 'properly paired' in line:
                                self.stat_dict[s]['sorted_spikein'] = int(line.split(' ')[0])
                    with open(f'{self.flagstat_path}BAM_rmdup/{self.dataset}_{self.ref_sp}/samtools_flagstat/{s}_{self.ref_sp}_sorted_readgps_rmdup_samtools_flagstat.txt') as f:
                        for line in f:
                            if 'properly paired' in line:
                                self.stat_dict[s]['rmdup_spikein'] = int(line.split(' ')[0])
                except FileNotFoundError:
                    return None
            try:
                with open(f'{self.flagstat_path}BAM_sorted/{self.dataset}_{self.ref}/samtools_flagstat/{s}_{self.ref}_sorted_readgps_samtools_flagstat.txt') as f:
                    for line in f:
                        if 'properly paired' in line:
                            self.stat_dict[s]['sorted'] = int(line.split(' ')[0])
                with open(f'{self.flagstat_path}BAM_rmdup/{self.dataset}_{self.ref}/samtools_flagstat/{s}_{self.ref}_sorted_readgps_rmdup_samtools_flagstat.txt') as f:
                    for line in f:
                        if 'properly paired' in line:
                            self.stat_dict[s]['rmdup'] = int(line.split(' ')[0])
            except FileNotFoundError:
                return None
                                        
        ic(self.stat_dict)
        
        


def main():
    #s = Stats('ChIPseq_221118_JR', '/mnt/storage/dept/pedonc/CPCT/projects/Partner_lab_projects/bedtools_coverage/ChIPseq_221118_JR_hg38/sorted_rmdup/tss_all_samples.csv', True)
    #s = Stats('ChIPseq_221118_JR','/mnt/storage/dept/pedonc/CPCT/projects/Partner_lab_projects/BAM_sorted/ChIPseq_221118_JR_hg38/samtools_flagstat/A549_ChIP_input_10ug_dex_10uM_2h_REP1_221118_TGACCA_S4_hg38_sorted_readgps_samtools_flagstat.txt', True)
    #s = FlagStats('RNAseq_2311x2', '/mnt/storage/dept/pedonc/CPCT/projects/MCF7_NO-202204/BAM_rmdup/RNAseq_2311x2_hg38/samtools_flagstat/', True)
    s = FlagStats('RNAseq_2311x2', '/mnt/storage/dept/pedonc/CPCT/projects/MCF7_NO-202204/', 'hg38')
    s.parse_flagstat()
    print(s)

if __name__ == '__main__':
    main()
