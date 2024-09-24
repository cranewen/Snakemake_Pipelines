import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
from collections import defaultdict
import yaml
import argparse
import time
from metaparser.dirbuilder import DirBuilder
from metaparser.sampleparser import SampleParser
from metaparser.yamlhandle import YamlHandle
from metaparser.metaconf import BedtoolSummaryColorConfig
from db.dbconnection import DBconnector
from db.cpctdb import DatasetSeqRunSample, DatasetSeqRun, Sample, SampleGroup, SampleAttribute
from datetime import datetime

# parsing flagstat from sorted, rmdup with both hg38 and dm6 data
class SamtoolsFlagstat:
    def __init__(self, dataset_name, ref_name):
        self.db_connection = DBconnector()
        self.dataset = dataset_name
        self.ref_name = ref_name
        self.dirs = DirBuilder(dataset_name)
        self.dirs.build_chipseq_dirs(ref_name)
        self.sorted_bam_flagstat_dir = self.dirs.sorted_bam_flagstat_dir
        self.rmdup_bam_flagstat_dir = self.dirs.rmdup_bam_flagstat_dir
        self.dirs.build_chipseq_spikein_dirs('dm6')
        self.sorted_bam_flagstat_dir_dm6 = self.dirs.sorted_bam_flagstat_dir
        self.rmdup_bam_flagstat_dir_dm6 = self.dirs.rmdup_bam_flagstat_dir
        sample_parser = SampleParser()
        sample_parser.get_sample_list(dataset_name)
        self.samples = sample_parser.sample_list
        sample_parser.get_sample_input_dict(dataset_name)
        self.sample_input_dict = sample_parser.sample_input_dict
        # {input:sample} as the group
        self.groups = defaultdict()
        self.flagstat_dict = defaultdict()
        # where the tss_all_samples.csv is at
        self.bedtools_coverage_dataset_rmdup_dir = self.dirs.bedtools_coverage_dataset_rmdup_dir
        self.tss_all_samples_path = self.bedtools_coverage_dataset_rmdup_dir + '/tss_all_samples.csv'
        self.bedtool_summary_config = self.dirs.bedtool_summary_config
        self.bedtool_summary_strategy = self.bedtool_summary_config['strategy']
        # color dict from yaml directly
        self.bedtool_summary_strategy_color = self.bedtool_summary_config['color']
        # the reformed color dict for final writing the spreadsheet
        # it's a dict of {each sample: color}
        self.color_dict = defaultdict()
        # RGB color dict
        self.color_config_dict = defaultdict()
        # groups within the inputs
        self.input_groups = self.bedtool_summary_config['input_group']
        self.raw_counts_df = None
        # {sample_cpct_name with S#: modified sample name just for the excel summary}
        self.sample_name_map = defaultdict()
        self.sample_list = [] # sample simplified names with the final order putting into spreadsheet
        # row names as the keys
        self.sample_reads_dict = defaultdict()
        # column names as the keys
        self.sample_reads_df_dict = defaultdict()
        # for saving to spreadsheet
        self.sample_reads_df = None


        self.rpk_count_ratios_df = pd.DataFrame()
        self.rpk_count_ratios_per_gene_df = None

        self.flyreads = None
        self.flyreads_only = None

        self.input_list = []
        self.rpk_count_ratios_colnames_dict = defaultdict()
        self.rpk_count_ratios_colnames_list = []

        self.per_gene_dict = defaultdict()
        
        self.sample_reads_np = None
        # numpy members, they map to any member with the same name but with df instead
        self.raw_counts_np = None   # paired with self.raw_counts_df
        self.raw_counts_matrix = None
        self.rpk_count_ratio_matrix = None
        self.rpk_count_ratios_per_gene_matrix = None

        # for coloring
        self.rpk_name_map = defaultdict()
        self.rpk_per_gene_name_map = defaultdict()
        self.rpk_count_ratio_color_dict = defaultdict()
        # per_gene color dictionary is the same as rpk_count_ratio_color_dict
        # self.rpk_count_ratio_per_gene_color_dict = defaultdict()
        self.sample_reads_color_dict = defaultdict()
        self.raw_counts_color_dict = defaultdict()
        
        self.sample_sheet_tab_dir = self.dirs.project_path
        self.summary_dir = self.dirs.bedtools_coverage_dataset_dir

        
    # parse flagstat and put it in a dict {sample_name:['paired in sequencing', 'properly paired']
    def parse_flagstat(self):
        self.flagstat_dict['sorted'] = defaultdict()
        self.flagstat_dict['rmdup'] = defaultdict()
        # for dm6 there is only 1 element in the array, because we only need the spike-in properly paired data
        # so it will be {sample_name:['properly paired']}
        self.flagstat_dict['sorted_dm6'] = defaultdict()
        self.flagstat_dict['rmdup_dm6'] = defaultdict()

        self.raw_counts_df = pd.read_csv(self.tss_all_samples_path, sep=',')
        sample_names = self.raw_counts_df.columns.values[5:]

        for k,v in self.bedtool_summary_strategy.items():
            self.sample_list.append(k)
            if not isinstance(v, str):
                for vi in v:
                    self.sample_list.append(vi)
                    
        
        for sample in sample_names:
            for k,v in self.bedtool_summary_strategy.items():
                if sample.split('_')[-1] == k.split('_')[-1]:
                    self.sample_name_map[sample] = k
                for vi in v:
                    if sample.split('_')[-1] == vi.split('_')[-1]:
                        self.sample_name_map[sample] = vi
 
        # {sample:input} as the dict for calculating "rmdup input normalization factor (6a)" etc.
        self.groups = {self.sample_name_map[k]:self.sample_name_map[v] for k,v in self.sample_input_dict.items()}

        # parsing only properly paired because spike-in share all paired reads with the original alignment
        def prop_paired(self, filepath, sample):
            sample_stats_dict = defaultdict()
            try:
                with open(filepath) as f:
                    for line in f:
                        if 'properly paired' in line:
                            sample_stats_dict[self.sample_name_map[sample]] = [int(line.split(' ')[0])]
            except FileNotFoundError:
                return sample_stats_dict
            return sample_stats_dict

        # parsing properly paired and paired in sequencing
        def both_paired(self, filepath, sample):
            with open(filepath) as f:
                sample_stats_dict = defaultdict()
                sample_stats_dict[self.sample_name_map[sample]] = []
                for line in f:
                    if 'properly paired' in line or 'paired in sequencing' in line:
                        sample_stats_dict[self.sample_name_map[sample]].append(int(line.split(' ')[0]))
            return sample_stats_dict

        for s in self.samples:
            # sorted_bam hg38 or mm10 and dm6
            sorted_flagstat = self.sorted_bam_flagstat_dir + s + '_' + self.ref_name + '_sorted_readgps_samtools_flagstat.txt'
            sorted_flagstat_dm6 = self.sorted_bam_flagstat_dir_dm6 + s + '_dm6_sorted_readgps_samtools_flagstat.txt'
            self.flagstat_dict['sorted'] = {**both_paired(self, sorted_flagstat, s), **self.flagstat_dict['sorted']}
            if self.flyreads or self.flyreads_only:
                self.flagstat_dict['sorted_dm6'] = {**prop_paired(self, sorted_flagstat_dm6, s), **self.flagstat_dict['sorted_dm6']}

            # rmdup_bam hg38 or mm10 and dm6
            rmdup_flagstat = self.rmdup_bam_flagstat_dir + s + '_' + self.ref_name + '_sorted_readgps_rmdup_samtools_flagstat.txt'
            rmdup_flagstat_dm6 = self.rmdup_bam_flagstat_dir_dm6 + s + '_dm6_sorted_readgps_rmdup_samtools_flagstat.txt'
            self.flagstat_dict['rmdup'] = {**both_paired(self, rmdup_flagstat, s), **self.flagstat_dict['rmdup']}
            if self.flyreads or self.flyreads_only:
                self.flagstat_dict['rmdup_dm6'] = {**prop_paired(self, rmdup_flagstat_dm6, s), **self.flagstat_dict['rmdup_dm6']}



    # prepare all the properties for sample_reads (e.g. sample_reads_dict, sample_reads_df...)            
    def prep_sample_reads(self): 
        rowname_list = ['rmdup paired reads (1)', 'rmdup % total', 'rmdup paired read ratio (2)', 
                        'rmdup ' + self.ref_name + ' aligned reads (1)', 'rmdup ' + self.ref_name + ' aligned read ratio (4)', 
                        'rmdup dm6 aligned reads (1)', 'rmdup all aligned reads (3)', 'rmdup % aligned', 
                        'rmdup spikein reads per million (5)', 'rmdup spikein rpm ratio (5b)', 
                        'rmdup input normalization factor (6a)', 'rmdup input-normalized ' + self.ref_name + ' reads (7a)', 
                        'rmdup input-normalized read ratio (8a)', 'rmdup read-adj input-norm read ratio (9a)', 
                        'rmdup read-adj spikein read ratio (9b)', 'rmdup final norm factor (input-norm) (10a)', 
                        'rmdup expected ' +  self.ref_name + ' aligned reads (input-norm) (11a)', 
                        'rmdup expected ' + self.ref_name + ' aligned reads (spikein rpm) (11b)', 'dupseq paired reads (1)', 
                        'dupseq paired read ratio (2)', 'dupseq ' + self.ref_name + ' aligned reads (1)', 'dupseq ' + self.ref_name + ' aligned read ratio (4)', 
                        'dupseq dm6 aligned reads (1)', 'dupseq all aligned reads (3)', 'dupseq % aligned', 
                        'dupseq spikein reads per million (5)', 'dupseq spikein rpm ratio (5b)', 'dupseq input normalization factor (6a)', 
                        'dupseq input-normalized ' + self.ref_name + ' reads (7a)', 'dupseq input-normalized read ratio (8a)', 
                        'dupseq read-adj input-norm read ratio (9a)', 'dupseq read-adj spikein read ratio (9b)', 
                        'dupseq final norm factor (input-norm) (10a)', 'dupseq expected ' + self.ref_name + ' aligned reads (input-norm) (11a)', 
                        'dupseq expected ' + self.ref_name + ' aligned reads (spikein rpm) (11b)'] 

        rowname_n = len(rowname_list)
        rowname_dict = {i+1:rowname_list[i] for i in range(rowname_n)}

        # notes: position 1, 4, 6, 19, 21, 23 are raw data
        self.sample_reads_dict = {r:{v:None for k,v in self.sample_name_map.items()} for r in rowname_list}

        self.sample_reads_df_dict = {v:[0]*rowname_n for v in self.sample_list}

                

        # To form a sample_reads df, calculate values by groups(strategy)

        if self.flyreads or self.flyreads_only:
            for k,v in self.bedtool_summary_strategy.items():
                if v != 'input':
                    self.sample_reads_dict[rowname_dict[1]][k] = self.flagstat_dict['rmdup'][k][0]
                    self.sample_reads_dict[rowname_dict[4]][k] = self.flagstat_dict['rmdup'][k][1]
                    self.sample_reads_dict[rowname_dict[19]][k] = self.flagstat_dict['sorted'][k][0]
                    self.sample_reads_dict[rowname_dict[21]][k] = self.flagstat_dict['sorted'][k][1]
                    self.sample_reads_dict[rowname_dict[6]][k] = self.flagstat_dict['rmdup_dm6'][k][0]
                    self.sample_reads_dict[rowname_dict[23]][k] = self.flagstat_dict['sorted_dm6'][k][0]
                    self.sample_reads_dict[rowname_dict[2]][k] = self.sample_reads_dict[rowname_dict[1]][k] / self.sample_reads_dict[rowname_dict[19]][k] * 100
                    self.sample_reads_dict[rowname_dict[7]][k] = self.sample_reads_dict[rowname_dict[4]][k] + self.sample_reads_dict[rowname_dict[6]][k]
                    self.sample_reads_dict[rowname_dict[8]][k] = self.sample_reads_dict[rowname_dict[7]][k] / self.sample_reads_dict[rowname_dict[1]][k] *100
                    self.sample_reads_dict[rowname_dict[9]][k] = self.sample_reads_dict[rowname_dict[6]][k] / self.sample_reads_dict[rowname_dict[7]][k] *1000000
                    self.sample_reads_dict[rowname_dict[11]][k] = self.sample_reads_dict[rowname_dict[9]][k] / self.sample_reads_dict[rowname_dict[9]][self.groups[k]]
                    self.sample_reads_dict[rowname_dict[12]][k] = self.sample_reads_dict[rowname_dict[4]][k] / self.sample_reads_dict[rowname_dict[11]][k]
                    self.sample_reads_dict[rowname_dict[17]][k] = self.sample_reads_dict[rowname_dict[4]][k]
                    self.sample_reads_dict[rowname_dict[18]][k] = self.sample_reads_dict[rowname_dict[4]][k]
                    self.sample_reads_dict[rowname_dict[24]][k] = self.sample_reads_dict[rowname_dict[21]][k] + self.sample_reads_dict[rowname_dict[23]][k]
                    self.sample_reads_dict[rowname_dict[25]][k] = self.sample_reads_dict[rowname_dict[24]][k] / self.sample_reads_dict[rowname_dict[19]][k] * 100
                    self.sample_reads_dict[rowname_dict[26]][k] = self.sample_reads_dict[rowname_dict[23]][k] / self.sample_reads_dict[rowname_dict[24]][k] * 1000000
                    self.sample_reads_dict[rowname_dict[28]][k] = self.sample_reads_dict[rowname_dict[26]][k] / self.sample_reads_dict[rowname_dict[26]][self.groups[k]]
                    self.sample_reads_dict[rowname_dict[29]][k] = self.sample_reads_dict[rowname_dict[21]][k] / self.sample_reads_dict[rowname_dict[28]][k]
                    self.sample_reads_dict[rowname_dict[34]][k] = self.sample_reads_dict[rowname_dict[21]][k]
                    self.sample_reads_dict[rowname_dict[35]][k] = self.sample_reads_dict[rowname_dict[21]][k]
                    
                    for vi in v:
                        self.sample_reads_dict[rowname_dict[1]][vi] = self.flagstat_dict['rmdup'][vi][0]
                        self.sample_reads_dict[rowname_dict[4]][vi] = self.flagstat_dict['rmdup'][vi][1]
                        self.sample_reads_dict[rowname_dict[19]][vi] = self.flagstat_dict['sorted'][vi][0]
                        self.sample_reads_dict[rowname_dict[21]][vi] = self.flagstat_dict['sorted'][vi][1]
                        self.sample_reads_dict[rowname_dict[6]][vi] = self.flagstat_dict['rmdup_dm6'][vi][0]
                        self.sample_reads_dict[rowname_dict[23]][vi] = self.flagstat_dict['sorted_dm6'][vi][0]
                        self.sample_reads_dict[rowname_dict[2]][vi] = self.sample_reads_dict[rowname_dict[1]][vi] / self.sample_reads_dict[rowname_dict[19]][vi] * 100
                        self.sample_reads_dict[rowname_dict[3]][vi] = self.sample_reads_dict[rowname_dict[1]][vi] / self.sample_reads_dict[rowname_dict[1]][k]
                        self.sample_reads_dict[rowname_dict[5]][vi] = self.sample_reads_dict[rowname_dict[4]][vi] / self.sample_reads_dict[rowname_dict[4]][k]
                        
                        self.sample_reads_dict[rowname_dict[7]][vi] = self.sample_reads_dict[rowname_dict[4]][vi] + self.sample_reads_dict[rowname_dict[6]][vi]

                        self.sample_reads_dict[rowname_dict[8]][vi] = self.sample_reads_dict[rowname_dict[7]][vi] / self.sample_reads_dict[rowname_dict[1]][vi] *100
                        self.sample_reads_dict[rowname_dict[9]][vi] = self.sample_reads_dict[rowname_dict[6]][vi] / self.sample_reads_dict[rowname_dict[7]][vi] *1000000
                        self.sample_reads_dict[rowname_dict[10]][vi] = self.sample_reads_dict[rowname_dict[9]][vi] / self.sample_reads_dict[rowname_dict[9]][k]
                        self.sample_reads_dict[rowname_dict[11]][vi] = self.sample_reads_dict[rowname_dict[9]][vi] / self.sample_reads_dict[rowname_dict[9]][self.groups[vi]]
                        self.sample_reads_dict[rowname_dict[12]][vi] = self.sample_reads_dict[rowname_dict[4]][vi] / self.sample_reads_dict[rowname_dict[11]][vi]
                        self.sample_reads_dict[rowname_dict[13]][vi] = self.sample_reads_dict[rowname_dict[12]][vi] / self.sample_reads_dict[rowname_dict[12]][k]
                        self.sample_reads_dict[rowname_dict[14]][vi] = self.sample_reads_dict[rowname_dict[13]][vi] / self.sample_reads_dict[rowname_dict[3]][vi]
                        self.sample_reads_dict[rowname_dict[15]][vi] = (1 / self.sample_reads_dict[rowname_dict[10]][vi]) / self.sample_reads_dict[rowname_dict[3]][vi]
                        self.sample_reads_dict[rowname_dict[16]][vi] = self.sample_reads_dict[rowname_dict[14]][vi] / self.sample_reads_dict[rowname_dict[5]][vi]
                        self.sample_reads_dict[rowname_dict[17]][vi] = self.sample_reads_dict[rowname_dict[4]][vi] * self.sample_reads_dict[rowname_dict[16]][vi]
                        self.sample_reads_dict[rowname_dict[18]][vi] = self.sample_reads_dict[rowname_dict[4]][vi] * self.sample_reads_dict[rowname_dict[15]][vi]
                        self.sample_reads_dict[rowname_dict[20]][vi] = self.sample_reads_dict[rowname_dict[19]][vi] / self.sample_reads_dict[rowname_dict[19]][k]
                        self.sample_reads_dict[rowname_dict[22]][vi] = self.sample_reads_dict[rowname_dict[21]][vi] / self.sample_reads_dict[rowname_dict[21]][k]
                        self.sample_reads_dict[rowname_dict[24]][vi] = self.sample_reads_dict[rowname_dict[21]][vi] + self.sample_reads_dict[rowname_dict[23]][vi]
                        self.sample_reads_dict[rowname_dict[25]][vi] = self.sample_reads_dict[rowname_dict[24]][vi] / self.sample_reads_dict[rowname_dict[19]][vi] * 100
                        self.sample_reads_dict[rowname_dict[26]][vi] = self.sample_reads_dict[rowname_dict[23]][vi] / self.sample_reads_dict[rowname_dict[24]][vi] * 1000000
                        self.sample_reads_dict[rowname_dict[27]][vi] = self.sample_reads_dict[rowname_dict[26]][vi] / self.sample_reads_dict[rowname_dict[26]][k]
                        self.sample_reads_dict[rowname_dict[28]][vi] = self.sample_reads_dict[rowname_dict[26]][vi] / self.sample_reads_dict[rowname_dict[26]][self.groups[vi]]
                        self.sample_reads_dict[rowname_dict[29]][vi] = self.sample_reads_dict[rowname_dict[21]][vi] / self.sample_reads_dict[rowname_dict[28]][vi]
                        self.sample_reads_dict[rowname_dict[30]][vi] = self.sample_reads_dict[rowname_dict[29]][vi] / self.sample_reads_dict[rowname_dict[29]][k]
                        self.sample_reads_dict[rowname_dict[31]][vi] = self.sample_reads_dict[rowname_dict[30]][vi] / self.sample_reads_dict[rowname_dict[20]][vi]
                        self.sample_reads_dict[rowname_dict[32]][vi] = (1 / self.sample_reads_dict[rowname_dict[27]][vi]) / self.sample_reads_dict[rowname_dict[20]][vi]
                        self.sample_reads_dict[rowname_dict[33]][vi] = self.sample_reads_dict[rowname_dict[31]][vi] / self.sample_reads_dict[rowname_dict[22]][vi]
                        self.sample_reads_dict[rowname_dict[34]][vi] = self.sample_reads_dict[rowname_dict[21]][vi] * self.sample_reads_dict[rowname_dict[33]][vi]
                        self.sample_reads_dict[rowname_dict[35]][vi] = self.sample_reads_dict[rowname_dict[21]][vi] * self.sample_reads_dict[rowname_dict[32]][vi]

                else:
                    self.sample_reads_dict[rowname_dict[1]][k] = self.flagstat_dict['rmdup'][k][0]
                    self.sample_reads_dict[rowname_dict[4]][k] = self.flagstat_dict['rmdup'][k][1]
                    self.sample_reads_dict[rowname_dict[19]][k] = self.flagstat_dict['sorted'][k][0]
                    self.sample_reads_dict[rowname_dict[21]][k] = self.flagstat_dict['sorted'][k][1]
                    self.sample_reads_dict[rowname_dict[6]][k] = self.flagstat_dict['rmdup_dm6'][k][0]
                    self.sample_reads_dict[rowname_dict[23]][k] = self.flagstat_dict['sorted_dm6'][k][0]
                    self.sample_reads_dict[rowname_dict[2]][k] = self.sample_reads_dict[rowname_dict[1]][k] / self.sample_reads_dict[rowname_dict[19]][k] * 100
                    self.sample_reads_dict[rowname_dict[7]][k] = self.sample_reads_dict[rowname_dict[6]][k] + self.sample_reads_dict[rowname_dict[4]][k]
                    self.sample_reads_dict[rowname_dict[8]][k] = self.sample_reads_dict[rowname_dict[7]][k] / self.sample_reads_dict[rowname_dict[1]][k] * 100
                    self.sample_reads_dict[rowname_dict[9]][k] = self.sample_reads_dict[rowname_dict[6]][k] / self.sample_reads_dict[rowname_dict[7]][k] * 1000000
                    if self.input_groups.get(k):
                        self.sample_reads_dict[rowname_dict[3]][k] = self.sample_reads_dict[rowname_dict[1]][k] / self.sample_reads_dict[rowname_dict[1]][self.input_groups[k]]
                        self.sample_reads_dict[rowname_dict[5]][k] = self.sample_reads_dict[rowname_dict[4]][k] / self.sample_reads_dict[rowname_dict[4]][self.input_groups[k]]
                        self.sample_reads_dict[rowname_dict[10]][k] = self.sample_reads_dict[rowname_dict[9]][k] / self.sample_reads_dict[rowname_dict[9]][self.input_groups[k]]
                        self.sample_reads_dict[rowname_dict[15]][k] = 1 / (self.sample_reads_dict[rowname_dict[10]][k] * self.sample_reads_dict[rowname_dict[3]][k])
                        self.sample_reads_dict[rowname_dict[18]][k] = self.sample_reads_dict[rowname_dict[4]][k] * self.sample_reads_dict[rowname_dict[15]][k]
                        self.sample_reads_dict[rowname_dict[20]][k] = self.sample_reads_dict[rowname_dict[19]][k] / self.sample_reads_dict[rowname_dict[19]][self.input_groups[k]]
                        self.sample_reads_dict[rowname_dict[22]][k] = self.sample_reads_dict[rowname_dict[21]][k] / self.sample_reads_dict[rowname_dict[21]][self.input_groups[k]]
                    else:
                        self.sample_reads_dict[rowname_dict[3]][k] = None
                        self.sample_reads_dict[rowname_dict[5]][k] = None
                        self.sample_reads_dict[rowname_dict[10]][k] = None
                        self.sample_reads_dict[rowname_dict[18]][k] = self.sample_reads_dict[rowname_dict[4]][k]
                    self.sample_reads_dict[rowname_dict[24]][k] = self.sample_reads_dict[rowname_dict[21]][k] + self.sample_reads_dict[rowname_dict[23]][k]
                    self.sample_reads_dict[rowname_dict[25]][k] = self.sample_reads_dict[rowname_dict[24]][k] / self.sample_reads_dict[rowname_dict[19]][k] * 100
                    self.sample_reads_dict[rowname_dict[26]][k] = self.sample_reads_dict[rowname_dict[23]][k] / self.sample_reads_dict[rowname_dict[24]][k] * 1000000
                    
                    if self.input_groups.get(k):
                        self.sample_reads_dict[rowname_dict[27]][k] = self.sample_reads_dict[rowname_dict[26]][k] / self.sample_reads_dict[rowname_dict[26]][self.input_groups[k]]
                        self.sample_reads_dict[rowname_dict[32]][k] = 1 / (self.sample_reads_dict[rowname_dict[27]][k] * self.sample_reads_dict[rowname_dict[20]][k])
                        self.sample_reads_dict[rowname_dict[35]][k] = self.sample_reads_dict[rowname_dict[21]][k] * self.sample_reads_dict[rowname_dict[32]][k]
                    else:
                        self.sample_reads_dict[rowname_dict[27]][k] = None
                  
                        self.sample_reads_dict[rowname_dict[35]][k] = self.sample_reads_dict[rowname_dict[21]][k]
        else:
            for k,v in self.bedtool_summary_strategy.items():
                if v != 'input':
                    self.sample_reads_dict[rowname_dict[1]][k] = self.flagstat_dict['rmdup'][k][0]
                    self.sample_reads_dict[rowname_dict[4]][k] = self.flagstat_dict['rmdup'][k][1]
                    self.sample_reads_dict[rowname_dict[19]][k] = self.flagstat_dict['sorted'][k][0]
                    self.sample_reads_dict[rowname_dict[21]][k] = self.flagstat_dict['sorted'][k][1]
                    self.sample_reads_dict[rowname_dict[2]][k] = self.sample_reads_dict[rowname_dict[1]][k] / self.sample_reads_dict[rowname_dict[19]][k] * 100
                    self.sample_reads_dict[rowname_dict[7]][k] = self.sample_reads_dict[rowname_dict[4]][k]
                    self.sample_reads_dict[rowname_dict[8]][k] = self.sample_reads_dict[rowname_dict[7]][k] / self.sample_reads_dict[rowname_dict[1]][k] *100
                    self.sample_reads_dict[rowname_dict[24]][k] = self.sample_reads_dict[rowname_dict[21]][k]
                    self.sample_reads_dict[rowname_dict[25]][k] = self.sample_reads_dict[rowname_dict[24]][k] / self.sample_reads_dict[rowname_dict[19]][k] * 100
                    
                    for vi in v:
                        self.sample_reads_dict[rowname_dict[1]][vi] = self.flagstat_dict['rmdup'][vi][0]
                        self.sample_reads_dict[rowname_dict[4]][vi] = self.flagstat_dict['rmdup'][vi][1]
                        self.sample_reads_dict[rowname_dict[19]][vi] = self.flagstat_dict['sorted'][vi][0]
                        self.sample_reads_dict[rowname_dict[21]][vi] = self.flagstat_dict['sorted'][vi][1]
                        self.sample_reads_dict[rowname_dict[2]][vi] = self.sample_reads_dict[rowname_dict[1]][vi] / self.sample_reads_dict[rowname_dict[19]][vi] * 100
                        self.sample_reads_dict[rowname_dict[3]][vi] = self.sample_reads_dict[rowname_dict[1]][vi] / self.sample_reads_dict[rowname_dict[1]][k]
                        self.sample_reads_dict[rowname_dict[5]][vi] = self.sample_reads_dict[rowname_dict[4]][vi] / self.sample_reads_dict[rowname_dict[4]][k]
                        self.sample_reads_dict[rowname_dict[7]][vi] = self.sample_reads_dict[rowname_dict[4]][vi]
                        self.sample_reads_dict[rowname_dict[8]][vi] = self.sample_reads_dict[rowname_dict[7]][vi] / self.sample_reads_dict[rowname_dict[1]][vi] *100
                        self.sample_reads_dict[rowname_dict[20]][vi] = self.sample_reads_dict[rowname_dict[19]][vi] / self.sample_reads_dict[rowname_dict[19]][k]
                        self.sample_reads_dict[rowname_dict[22]][vi] = self.sample_reads_dict[rowname_dict[21]][vi] / self.sample_reads_dict[rowname_dict[21]][k]
                        self.sample_reads_dict[rowname_dict[24]][vi] = self.sample_reads_dict[rowname_dict[21]][vi]
                        self.sample_reads_dict[rowname_dict[25]][vi] = self.sample_reads_dict[rowname_dict[24]][vi] / self.sample_reads_dict[rowname_dict[19]][vi] * 100

                else:
                    self.sample_reads_dict[rowname_dict[1]][k] = self.flagstat_dict['rmdup'][k][0]
                    self.sample_reads_dict[rowname_dict[4]][k] = self.flagstat_dict['rmdup'][k][1]
                    self.sample_reads_dict[rowname_dict[19]][k] = self.flagstat_dict['sorted'][k][0]
                    self.sample_reads_dict[rowname_dict[21]][k] = self.flagstat_dict['sorted'][k][1]
                    self.sample_reads_dict[rowname_dict[2]][k] = self.sample_reads_dict[rowname_dict[1]][k] / self.sample_reads_dict[rowname_dict[19]][k] * 100
                    self.sample_reads_dict[rowname_dict[7]][k] = self.sample_reads_dict[rowname_dict[4]][k]
                    self.sample_reads_dict[rowname_dict[8]][k] = self.sample_reads_dict[rowname_dict[7]][k] / self.sample_reads_dict[rowname_dict[1]][k] * 100
                    self.sample_reads_dict[rowname_dict[24]][k] = self.sample_reads_dict[rowname_dict[21]][k]
                    self.sample_reads_dict[rowname_dict[25]][k] = self.sample_reads_dict[rowname_dict[24]][k] / self.sample_reads_dict[rowname_dict[19]][k] * 100
                    

        rounding_factors_decimal = [2, 3, 5, 8, 10, 11, 13, 14, 15, 16, 20, 22, 25, 27, 28, 30, 31, 32, 33]
        range_n = range(1, rowname_n+1)
        rounding_factors_others = list(set(range_n) - (set(range_n) & set(rounding_factors_decimal)))
        rounding_factors_int = [9, 26]
        for i in range(rowname_n):
            for k in self.sample_reads_df_dict.keys():
                self.sample_reads_df_dict[k][i] = self.sample_reads_dict[rowname_dict[i+1]][k]
        #print(self.sample_reads_df_dict)
        # using ' ' as the key just for printing purpose
        self.sample_reads_df_dict = {**{' ': rowname_list}, **self.sample_reads_df_dict}
        self.sample_reads_df = pd.DataFrame(self.sample_reads_df_dict)
        sample_reads_matrix = self.sample_reads_df.to_numpy().T
        n_srm = len(sample_reads_matrix)
        for i in range(1, n_srm):
            for j in rounding_factors_decimal:
                if not np.isnan(sample_reads_matrix[i][j-1]):
                    sample_reads_matrix[i][j-1] = np.round(sample_reads_matrix[i][j-1], 2)
            for k in rounding_factors_others:
                if not np.isnan(sample_reads_matrix[i][k-1]):
                    sample_reads_matrix[i][k-1] = np.round(sample_reads_matrix[i][k-1])
                
        
        n_sample = len(self.sample_list)
        sample_reads_df_dict = {self.sample_list[i]:sample_reads_matrix[i+1] for i in range(n_sample)}
        sample_reads_df_dict = {**{' ': rowname_list}, **sample_reads_df_dict}
        
        self.sample_reads_df = pd.DataFrame(sample_reads_df_dict)
        #self.sample_reads_df.to_csv('test.csv', index=False)
        
    '''        
    def get_raw_counts_df(self):
        self.raw_counts_df = pd.read_csv(tss_all_samples_path, sep=',', header=0)
        self.raw_counts_sample_names = self.raw_counts_df.columns.values[5:]
    '''        
        
    def reorder_raw_counts_df(self):
        self.raw_counts_df.rename(self.sample_name_map, axis=1, inplace=True)
        self.raw_counts_df['Gene'] = [x[:-(len(x.split('-')[-1])+1)] for x in self.raw_counts_df['Interval']] 
        column_new_order = ['Interval', 'Gene', 'Chr', 'Start', 'End'] + list(self.sample_list)
        self.raw_counts_df = self.raw_counts_df[column_new_order]
        #self.raw_counts_df.to_csv('tab2.csv', sep=',', index=False)
        
    # if there is flyreads included
    def set_flyreads(self, flyreads=None, flyreads_only=None):
        self.flyreads = flyreads
        self.flyreads_only = flyreads_only

    # making the rpk_count_ratios tab dataframe
    def make_rpk_count_ratios(self):
        self.rpk_count_ratios_df = self.raw_counts_df[['Interval', 'Gene', 'Chr', 'Start', 'End']]
        self.rpk_count_ratios_df['Length'] = [4000] * len(self.raw_counts_df['Gene'])
        # for input in self.input_list:
        #     self.rpk_count_ratios_df[input + '_rpk'] = self.raw_counts_df[input].to_numpy() / (self.rpk_count_ratios_df['Length'].to_numpy() / 1000)
        if self.flyreads_only:
            for k,v in self.bedtool_summary_strategy.items():
                self.rpk_name_map[k + '_rpk'] = k
                if v == 'input':
                    self.rpk_count_ratios_df[k + '_rpk'] = self.raw_counts_df[k].to_numpy() / (self.rpk_count_ratios_df['Length'].to_numpy() / 1000)
                else:
                    self.rpk_count_ratios_df[k + '_rpk'] = self.raw_counts_df[k].to_numpy() / (self.rpk_count_ratios_df['Length'].to_numpy() / 1000)
                    for vi in v:
                        self.rpk_count_ratios_df[vi + '_rpk'] = self.raw_counts_df[vi].to_numpy() / (self.rpk_count_ratios_df['Length'].to_numpy() / 1000)
                        self.rpk_count_ratios_df[vi + '_rpk_ratio'] = (self.rpk_count_ratios_df[vi + '_rpk'].to_numpy() + 1) / (self.rpk_count_ratios_df[k + '_rpk'].to_numpy() + 1)
                        dmso_rpk_bi = self.rpk_count_ratios_df[k + '_rpk'] >= 10
                        sample_rpk_bi = self.rpk_count_ratios_df[vi + '_rpk'] >= 10
                        self.rpk_count_ratios_df[vi + '_rpk_10rpkfloor'] = np.logical_or(dmso_rpk_bi, sample_rpk_bi)
                        self.rpk_count_ratios_df[vi + '_rpk_10rpkfloor'] = np.where(self.rpk_count_ratios_df[vi + '_rpk_10rpkfloor'] == True, 'Y', 'N')
                        self.rpk_count_ratios_df[vi + '_rpk_flynorm'] = self.rpk_count_ratios_df[vi + '_rpk'].to_numpy() * self.sample_reads_dict['rmdup read-adj spikein read ratio (9b)'][vi]
                        self.rpk_count_ratios_df[vi + '_rpk_flynorm_ratio'] = (self.rpk_count_ratios_df[vi + '_rpk_flynorm'].to_numpy() + 1) / (self.rpk_count_ratios_df[k + '_rpk'].to_numpy() + 1)
                        sample_rpk_flynorm_bi = self.rpk_count_ratios_df[vi + '_rpk_flynorm'] >= 10
                        self.rpk_count_ratios_df[vi + '_rpk_flynorm_10rpkfloor'] = np.logical_or(dmso_rpk_bi, sample_rpk_flynorm_bi)
                        self.rpk_count_ratios_df[vi + '_rpk_flynorm_10rpkfloor'] = np.where(self.rpk_count_ratios_df[vi + '_rpk_flynorm_10rpkfloor'] == True, 'Y', 'N')
                        self.rpk_name_map[vi + '_rpk'] = vi
                        self.rpk_name_map[vi + '_rpk_ratio'] = vi
                        self.rpk_name_map[vi + '_rpk_10rpkfloor'] = vi
                        self.rpk_name_map[vi + '_rpk_flynorm'] = vi
                        self.rpk_name_map[vi + '_rpk_flynorm_ratio'] = vi
                        self.rpk_name_map[vi + '_rpk_flynorm_10rpkfloor'] = vi

        elif self.flyreads:
            for k,v in self.bedtool_summary_strategy.items():
                self.rpk_name_map[k + '_rpk'] = k
                if v == 'input':
                    self.rpk_count_ratios_df[k + '_rpk'] = self.raw_counts_df[k].to_numpy() / (self.rpk_count_ratios_df['Length'].to_numpy() / 1000)
                else:
                    self.rpk_count_ratios_df[k + '_rpk'] = self.raw_counts_df[k].to_numpy() / (self.rpk_count_ratios_df['Length'].to_numpy() / 1000)
                    for vi in v:
                        self.rpk_count_ratios_df[vi + '_rpk'] = self.raw_counts_df[vi].to_numpy() / (self.rpk_count_ratios_df['Length'].to_numpy() / 1000)
                        self.rpk_count_ratios_df[vi + '_rpk_ratio'] = (self.rpk_count_ratios_df[vi + '_rpk'] + 1) / (self.rpk_count_ratios_df[k + '_rpk'] + 1)
                        dmso_rpk_bi = self.rpk_count_ratios_df[k + '_rpk'] >= 10
                        sample_rpk_bi = self.rpk_count_ratios_df[vi + '_rpk'] >= 10
                        self.rpk_count_ratios_df[vi + '_rpk_10rpkfloor'] = np.logical_or(dmso_rpk_bi, sample_rpk_bi)
                        self.rpk_count_ratios_df[vi + '_rpk_10rpkfloor'] = np.where(self.rpk_count_ratios_df[vi + '_rpk_10rpkfloor'] == True, 'Y', 'N')
                        self.rpk_count_ratios_df[vi + '_rpk_readnorm'] = self.rpk_count_ratios_df[vi + '_rpk'] / self.sample_reads_dict['rmdup ' + self.ref_name + ' aligned read ratio (4)'][vi]
                        self.rpk_count_ratios_df[vi + '_rpk_readnorm_ratio'] = (self.rpk_count_ratios_df[vi + '_rpk_readnorm'] + 1) / (self.rpk_count_ratios_df[k + '_rpk'] + 1)
                        sample_rpk_readnorm_bi = self.rpk_count_ratios_df[vi + '_rpk_readnorm'] >= 10
                        self.rpk_count_ratios_df[vi + '_rpk_readnorm_10rpkfloor'] = np.logical_or(dmso_rpk_bi, sample_rpk_readnorm_bi)
                        self.rpk_count_ratios_df[vi + '_rpk_readnorm_10rpkfloor'] = np.where(self.rpk_count_ratios_df[vi + '_rpk_readnorm_10rpkfloor'] == True, 'Y', 'N')
                        self.rpk_count_ratios_df[vi + '_rpk_flynorm'] = self.rpk_count_ratios_df[vi + '_rpk'] * self.sample_reads_dict['rmdup read-adj spikein read ratio (9b)'][vi]
                        self.rpk_count_ratios_df[vi + '_rpk_flynorm_ratio'] = (self.rpk_count_ratios_df[vi + '_rpk_flynorm'] + 1) / (self.rpk_count_ratios_df[k + '_rpk'] + 1)
                        sample_rpk_flynorm_bi = self.rpk_count_ratios_df[vi + '_rpk_flynorm'] >= 10
                        self.rpk_count_ratios_df[vi + '_rpk_flynorm_10rpkfloor'] = np.logical_or(dmso_rpk_bi, sample_rpk_flynorm_bi)
                        self.rpk_count_ratios_df[vi + '_rpk_flynorm_10rpkfloor'] = np.where(self.rpk_count_ratios_df[vi + '_rpk_flynorm_10rpkfloor'] == True, 'Y', 'N')
                        self.rpk_name_map[vi + '_rpk'] = vi
                        self.rpk_name_map[vi + '_rpk_ratio'] = vi
                        self.rpk_name_map[vi + '_rpk_10rpkfloor'] = vi
                        self.rpk_name_map[vi + '_rpk_readnorm'] = vi
                        self.rpk_name_map[vi + '_rpk_readnorm_ratio'] = vi
                        self.rpk_name_map[vi + '_rpk_readnorm_10rpkfloor'] = vi
                        self.rpk_name_map[vi + '_rpk_flynorm'] = vi
                        self.rpk_name_map[vi + '_rpk_flynorm_ratio'] = vi
                        self.rpk_name_map[vi + '_rpk_flynorm_10rpkfloor'] = vi

        else:
            for k,v in self.bedtool_summary_strategy.items():
                self.rpk_name_map[k + '_rpk'] = k
                if v == 'input':
                    self.rpk_count_ratios_df[k + '_rpk'] = self.raw_counts_df[k].to_numpy() / (self.rpk_count_ratios_df['Length'].to_numpy() / 1000)
                else:
                    self.rpk_count_ratios_df[k + '_rpk'] = self.raw_counts_df[k].to_numpy() / (self.rpk_count_ratios_df['Length'].to_numpy() / 1000)
                    for vi in v:
                        self.rpk_count_ratios_df[vi + '_rpk'] = self.raw_counts_df[vi].to_numpy() / (self.rpk_count_ratios_df['Length'].to_numpy() / 1000)
                        self.rpk_count_ratios_df[vi + '_rpk_ratio'] = (self.rpk_count_ratios_df[vi + '_rpk'].to_numpy() + 1) / (self.rpk_count_ratios_df[k + '_rpk'].to_numpy() + 1)
                        dmso_rpk_bi = self.rpk_count_ratios_df[k + '_rpk'] >= 10
                        sample_rpk_bi = self.rpk_count_ratios_df[vi + '_rpk'] >= 10
                        self.rpk_count_ratios_df[vi + '_rpk_10rpkfloor'] = np.logical_or(dmso_rpk_bi, sample_rpk_bi)
                        self.rpk_count_ratios_df[vi + '_rpk_10rpkfloor'] = np.where(self.rpk_count_ratios_df[vi + '_rpk_10rpkfloor'] == True, 'Y', 'N')
                        self.rpk_count_ratios_df[vi + '_rpk_readnorm'] = self.rpk_count_ratios_df[vi + '_rpk'] / self.sample_reads_dict['rmdup ' + self.ref_name + ' aligned read ratio (4)'][vi]
                        self.rpk_count_ratios_df[vi + '_rpk_readnorm_ratio'] = (self.rpk_count_ratios_df[vi + '_rpk_readnorm'].to_numpy() + 1) / (self.rpk_count_ratios_df[k + '_rpk'].to_numpy() + 1)
                        sample_rpk_readnorm_bi = self.rpk_count_ratios_df[vi + '_rpk_readnorm'] >= 10
                        self.rpk_count_ratios_df[vi + '_rpk_readnorm_10rpkfloor'] = np.logical_or(dmso_rpk_bi, sample_rpk_readnorm_bi)
                        self.rpk_count_ratios_df[vi + '_rpk_readnorm_10rpkfloor'] = np.where(self.rpk_count_ratios_df[vi + '_rpk_readnorm_10rpkfloor'] == True, 'Y', 'N')
                        self.rpk_name_map[vi + '_rpk'] = vi
                        self.rpk_name_map[vi + '_rpk_ratio'] = vi
                        self.rpk_name_map[vi + '_rpk_10rpkfloor'] = vi
                        self.rpk_name_map[vi + '_rpk_readnorm'] = vi
                        self.rpk_name_map[vi + '_rpk_readnorm_ratio'] = vi
                        self.rpk_name_map[vi + '_rpk_readnorm_10rpkfloor'] = vi

        #self.rpk_count_ratios_df.to_csv('tab3.csv', sep=',')
 
    def get_per_gene_dict(self):
        n = len(self.rpk_count_ratios_df['Gene'])
        gene_list = list(set(self.rpk_count_ratios_df['Gene']))
        self.per_gene_dict = {x:[] for x in gene_list}

        for i in range(n):
            self.per_gene_dict[self.rpk_count_ratios_df['Gene'][i]].append(i)


    # making the rpk_count_ratios_per_gene tab
    def make_rpk_count_ratios_per_gene(self):
        dmso_samples = [k+'_rpk' for k,v in self.bedtool_summary_strategy.items()]
        dmso_df = self.rpk_count_ratios_df[['Interval', 'Gene', 'Chr', 'Start', 'End'] + dmso_samples]
        dmso_df_matrix = dmso_df.to_numpy()[:,5:]
        dmso_samples_max_mean_dict = defaultdict()
        dmso_samples_max_mean_list = []
        for k,v in self.per_gene_dict.items():
            dmso_samples_max_mean_list.append(v[np.argmax(np.sum(dmso_df_matrix[self.per_gene_dict[k]], axis=1))])

        self.rpk_count_ratios_per_gene_df = self.rpk_count_ratios_df.iloc[dmso_samples_max_mean_list]
        #self.rpk_count_ratios_per_gene_df.to_csv('tab4.csv', sep=',')

    def excel_cell_num(self, num):
        c = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        results = []
        if num <= 26**2:
            for i in range(num):
                l1 = i // 26
                l2 = i % 26
                if l1 == 0:
                    results.append(c[l2] + "1")
                else:
                    results.append(c[l1 - 1] + c[l2] + "1")
        return results

    # picking color gradient steps
    # return a list of indice of the sub colors
    def color_picking(self, num):
        if num == 1:
            return 7 // 2
        return [x * (7 // num) for x in range(1, num + 1)]
        

    def color_reform(self):
        bedtool_summary_color_config_handle = YamlHandle(BedtoolSummaryColorConfig.BEDTOOL_SUMMARY_COLOR_CONFIG_YAML_PATH.value).read_yaml()
        self.color_config_dict = next(bedtool_summary_color_config_handle)

        input_dict = {v:[] for k,v in self.input_groups.items()}
        for k,v in self.input_groups.items():
            input_dict[v].append(k)
        for k,v in self.bedtool_summary_strategy.items():
            if not v == 'input':
                n = len(v)
                color_list = self.color_picking(n)
                self.color_dict[k] = self.color_config_dict[self.bedtool_summary_strategy_color[k]]['base']
                for i in range(n):
                    self.color_dict[v[i]] = self.color_config_dict[self.bedtool_summary_strategy_color[k]]['sub'][i]
        for k,v in input_dict.items():
            n = len(v)
            color_list = self.color_picking(n)
            self.color_dict[k] = self.color_config_dict[self.bedtool_summary_strategy_color[k]]['base']
            for i in range(n):
                self.color_dict[v[i]] = self.color_config_dict[self.bedtool_summary_strategy_color[k]]['sub'][i]


    # mapping color to specific cell. e.g. A1:color
    def get_cell_color_dict(self):
        rpk_count_ratio_df_cols = list(self.rpk_count_ratios_df.columns)
        n = len(rpk_count_ratio_df_cols)
        rpk_count_ratio_cell = self.excel_cell_num(n)
        for i in range(n):
            try:
                #[cell, color]
                self.rpk_count_ratio_color_dict[rpk_count_ratio_df_cols[i]] = [rpk_count_ratio_cell[i], self.color_dict[self.rpk_name_map[rpk_count_ratio_df_cols[i]]]]
            except KeyError:
                pass

        # for the 1st tab(sample_reads) coloring
        sample_reads_columns = list(self.sample_reads_df.columns)
        m = len(sample_reads_columns)
        sample_reads_cell = self.excel_cell_num(m)
        for i in range(m):
            try:
                self.sample_reads_color_dict[sample_reads_columns[i]] = [sample_reads_cell[i], self.color_dict[sample_reads_columns[i]]]
            except KeyError:
                pass

        # for the 2nd tab(raw_counts) coloring
        raw_counts_columns = list(self.raw_counts_df.columns)
        k = len(raw_counts_columns)
        raw_counts_cell = self.excel_cell_num(k)
        for i in range(k):
            try:
                self.raw_counts_color_dict[raw_counts_columns[i]] = [raw_counts_cell[i], self.color_dict[raw_counts_columns[i]]]
            except KeyError:
                pass

    def write_excel(self):
        date_today = datetime.today().strftime('%Y%m%d')
        sample_sheet_file_name = self.dataset + '_normalization_factors_' + date_today + '.xlsx'
        if self.flyreads:
            sample_sheet_file_name = self.dataset + '_normalization_factors_fly-' + date_today + '.xlsx'
        if self.flyreads_only:
            sample_sheet_file_name = self.dataset + '_flyonly-' + date_today + '.xlsx'

        summary_file_name = self.dataset + '_bedtools_coverage_TSS_onepergene_signal_ratio_summarys-' + date_today + '.xlsx'
        
        with pd.ExcelWriter(self.sample_sheet_tab_dir + sample_sheet_file_name, engine='xlsxwriter') as writer1:
            self.sample_reads_df.to_excel(writer1, sheet_name='sample_reads', index=False)
            workbook = writer1.book
            worksheet = writer1.sheets['sample_reads']

            a1_format = workbook.add_format()
            worksheet.write('A1', None, a1_format)
            digit_format = workbook.add_format({'num_format': '#,##0'})
            decimal_format = workbook.add_format({'num_format': '##0.00'})
            for k,v in self.sample_reads_color_dict.items():
                color_format = workbook.add_format({'bg_color': v[1]})
                worksheet.write(v[0], k, color_format)
            last_cell = self.excel_cell_num(len(self.sample_reads_df.columns))[-1][:-1]
            cell_range = 'B2:' + last_cell + str(len(self.sample_reads_df) + 1)
            worksheet.conditional_format(cell_range, {'type': 'cell', 'criteria': '>', 'value': 100, 'format': digit_format})
            worksheet.conditional_format(cell_range, {'type': 'cell', 'criteria': '<=', 'value': 100, 'format': decimal_format})

        print('===== writing normalization factor summary done =====')
             
        with pd.ExcelWriter(self.summary_dir + summary_file_name, engine='xlsxwriter') as writer:
            self.sample_reads_df.to_excel(writer, sheet_name='sample_reads', index=False)

            # change color section
            workbook = writer.book
            # adding colors to the 1st row on 1st tab
            worksheet = writer.sheets['sample_reads']

            a1_format = workbook.add_format()
            worksheet.write('A1', None, a1_format)
            digit_format = workbook.add_format({'num_format': '#,##0'})
            decimal_format = workbook.add_format({'num_format': '##0.00'})
            for k,v in self.sample_reads_color_dict.items():
                color_format = workbook.add_format({'bg_color': v[1]})
                worksheet.write(v[0], k, color_format)
            last_cell = self.excel_cell_num(len(self.sample_reads_df.columns))[-1][:-1]
            cell_range = 'B2:' + last_cell + str(len(self.sample_reads_df) + 1)
            worksheet.conditional_format(cell_range, {'type': 'cell', 'criteria': '>', 'value': 100, 'format': digit_format})
            worksheet.conditional_format(cell_range, {'type': 'cell', 'criteria': '<=', 'value': 100, 'format': decimal_format})
                
            print('===== writing sample_sheet tab done =====')

            # adding colors to the 1st row on 2nd tab
            self.raw_counts_df.to_excel(writer, sheet_name='raw_counts', index=False)
            worksheet = writer.sheets['raw_counts']
            for k,v in self.raw_counts_color_dict.items():
                color_format = workbook.add_format({'bg_color': v[1]})
                worksheet.write(v[0], k, color_format)

            print('===== writing raw_counts tab done =====')

            # adding colors to the 1st row on 3rd tab
            self.rpk_count_ratios_df.to_excel(writer, sheet_name='rpk_count_ratios', index=False)
            worksheet = writer.sheets['rpk_count_ratios']
            format1 = workbook.add_format({'num_format': '#,##0.00'})
            worksheet.set_column(7,len(self.rpk_count_ratios_df.columns),10,format1)
            for k,v in self.rpk_count_ratio_color_dict.items():
                color_format = workbook.add_format({'bg_color': v[1]})
                worksheet.write(v[0], k, color_format)

            print('===== writing rpk_count_ratios tab done =====')

            # adding colors to the 1st row on 4th tab
            self.rpk_count_ratios_per_gene_df.to_excel(writer, sheet_name='rpk_count_ratios_per_gene', index=False)
            worksheet = writer.sheets['rpk_count_ratios_per_gene']
            worksheet.set_column(7,len(self.rpk_count_ratios_per_gene_df.columns),10,format1)
            for k,v in self.rpk_count_ratio_color_dict.items():
                color_format = workbook.add_format({'bg_color': v[1]})
                worksheet.write(v[0], k, color_format)

            print('===== writing rpk_count_ratios_per_gene tab done =====')
            print('===== done =====')

class ConfigWriter:
    def __init__(self, dataset):
        self.dataset = dataset
        self.db_connection = DBconnector()

    def write_config(self):
        # match with bedtool_summary_config.yaml
        strategy_dict = defaultdict()
        input_group_dict = defaultdict()
        color_dict = defaultdict()
        with self.db_connection.create_session() as session:
            sample_group = session.query(SampleGroup)
            sample_attr = session.query(DatasetSeqRunSample, DatasetSeqRun, Sample, SampleAttribute, SampleGroup).outerjoin(SampleGroup, SampleGroup.sampleGroupID == DatasetSeqRunSample.sampleGroupID).filter(DatasetSeqRun.datasetSeqRunID == DatasetSeqRunSample.datasetSeqRunID).filter(DatasetSeqRun.datasetSeqRunDirectory == self.dataset).filter(SampleGroup.sampleGroupID == DatasetSeqRunSample.sampleGroupID).filter(Sample.sampleID == DatasetSeqRunSample.sampleID).filter(SampleAttribute.sampleID == Sample.sampleID).filter(SampleAttribute.sampleAttributeTypeID == 16).order_by(DatasetSeqRunSample.sampleGroupID, DatasetSeqRunSample.sampleGroupPosition).all()
            colors = session.query(SampleGroup)
            # a dict from SampleGroup {sampleGroupID: baseColor}
            color = defaultdict()
            for c in colors:
                color[c.sampleGroupID] = c.sampleGroupBaseColor
            color['input'] = color[1]
                
            temp = defaultdict()
            temp_inner = defaultdict()
            t_i = ''
            input_group_ids = []
            temp_color = None
            for s in sample_attr:
                len_s = len(s[3].sampleAttributeValue.split('_')[-1])
                if temp.get(s[0].sampleGroupID) is None:
                    temp_inner = defaultdict()
                # checking S#
                if s[3].sampleAttributeValue[-len_s:] != 'S' + str(s[0].fastqSampleNumber):
                    print(s[3].sampleAttributeValue + " SampleAttributeValue S# is wrong, it should be " + str(s[0].fastqSampleNumber))
                if s[0].chIPInputSampleFlag == 'Y':
                    input_group_ids.append(s[0].sampleGroupID)
                    if s[0].sampleGroupPosition == 1:
                        #t_i = s[3].sampleAttributeValue
                        t_i = s[3].sampleAttributeValue[:-len_s] + 'S' + str(s[0].fastqSampleNumber)
                    else:
                        input_group_dict[s[3].sampleAttributeValue[:-len_s] + 'S' + str(s[0].fastqSampleNumber)] = t_i


                if s[0].sampleGroupPosition == 1:
                    if s[0].sampleGroupID == 1 and s[0].chIPInputSampleFlag != 'Y':
                        temp_color = {k-1:v for (k,v) in color.items() if isinstance(k, int)}
                        temp_color['input'] = color[1]
                        color = temp_color
                    if s[0].chIPInputSampleFlag == 'Y':
                        color_dict[s[3].sampleAttributeValue[:-len_s] + 'S' + str(s[0].fastqSampleNumber)] = color['input']
                    else:
                        color_dict[s[3].sampleAttributeValue[:-len_s] + 'S' + str(s[0].fastqSampleNumber)] = color[s[0].sampleGroupID]

                temp_inner[s[0].sampleGroupPosition] = s[3].sampleAttributeValue[:-len_s] + 'S' + str(s[0].fastqSampleNumber)
                temp[s[0].sampleGroupID] = temp_inner
            input_group_ids = set(input_group_ids)
            for k,v in temp.items():
                if not k in input_group_ids:
                    for ki,vi in v.items():
                        if ki == 1:
                            strategy_dict[vi] = []
                        else:
                            strategy_dict[v[1]].append(vi)
                else:
                    for ki,vi in v.items():
                        strategy_dict[vi] = 'input'
        with open('metaconfig/bedtool_summary_config.yaml', 'a') as f:
            f.write(self.dataset + ':\n')
            f.write('  strategy:\n')
            for k,v in strategy_dict.items():
                f.write('    ' + k + ':\n')
                if isinstance(v, str):
                    f.write('      ' + '\"input\"\n')
                else:
                    for vi in v:
                        f.write('      - \"' + vi + '\"\n')
            
            f.write('  input_group:\n')
            for k,v in input_group_dict.items():
                f.write('    ' + k + ':\n')   
                f.write('      \"' + v + '\"\n')
            
            f.write('  color:\n')
            for k,v in color_dict.items():
                f.write('    ' + k + ':\n')   
                f.write('      \"' + v + '\"\n')
            f.write('\n\n')
    


def main(dataset, ref_genome, norm, tss):
    cw = ConfigWriter(dataset)
    cw.write_config()
    print('===== writing config file done =====')
    
    #sf = SamtoolsFlagstat('ChIPseq_221108_DW', 'hg38')
    sf = SamtoolsFlagstat(dataset, ref_genome)
    if norm == 'reads':
        sf.set_flyreads()
    if norm == 'fly':
        sf.set_flyreads(flyreads_only = 1)
    if norm == 'reads_fly':
        sf.set_flyreads(flyreads = 1)
    if tss:
        sf.tss_all_samples_path = sf.bedtools_coverage_dataset_rmdup_dir + '/' + tss
    sf.parse_flagstat()
    sf.prep_sample_reads()
    sf.reorder_raw_counts_df()
    sf.make_rpk_count_ratios()
    sf.get_per_gene_dict()
    sf.make_rpk_count_ratios_per_gene()
    sf.color_reform()
    sf.get_cell_color_dict()
    sf.write_excel()
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Bedtools summary')
    parser.add_argument('-d', '--dataset', type=str, required=True)
    parser.add_argument('-norm', '--normalization', type=str)
    parser.add_argument('-ref', '--ref_genome', type=str)
    parser.add_argument('-tss', '--tss_file', type=str)
    args = parser.parse_args()
    main(args.dataset, args.ref_genome, args.normalization, args.tss_file)

