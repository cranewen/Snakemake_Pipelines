from collections import defaultdict
import pandas as pd
import numpy as np
import sys, os
#sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
sys.path.insert(0, '/home/yw900/lab_pipelines')
import yaml
from metaparser.dirbuilder import DirBuilder
from metaparser.sampleparser import SampleParser
from metaparser.yamlhandle import YamlHandle
from metaparser.metaconf import BedtoolSummaryColorConfig
from db.dbconnection import DBconnector
from db.cpctdb import DatasetSeqRunSample, DatasetSeqRun, SampleGroup, SampleAttribute, SampleAttributeType, Sample as dbSample, Dataset as dbDataset, Project, Scientist, Protocol as dbProtocol
from datetime import datetime
from data import *
from utility import YamlReader
from typing import Dict
import random
from icecream import ic

# any function starts with '_' is the database version
class Database:
    def __init__(self, dataset):
        self.dataset: Dataset = dataset # passing a Dataset object when initialize Database
        self.db_connection: DBconnector = DBconnector()
        self.sample_attr: db.cpctdb = None
        # a dict {sample short name: Sample object}
        self.sample_meta_shortname: Dict = defaultdict()
        self.sample_meta: Dict = defaultdict() # key is cpctName with S#
        self.dataset_run: DatasetRun = None
        if 'RNAseq' in self.dataset.dataset_name:
            self.dataset.protocol = 'RNAseq'
        if 'ChIPseq' in self.dataset.dataset_name:
            self.dataset.protocol = 'ChIPseq'
        print(f'This is the class for processing all the data from MySQL')

    def _protocol_check(self):
        with self.db_connection.create_session() as session:
            protocol = session.query(dbProtocol, DatasetSeqRun).filter(DatasetSeqRun.protocolID == dbProtocol.protocolID).filter(DatasetSeqRun.datasetSeqRunDirectory == self.dataset.dataset_name).one()
            self.dataset.protocol = protocol[0].protocolName
            print(f'protocol ======= {protocol[0].protocolName}')
        

    def _get_sample_info_rna(self):
        data_p = DataParser(self.dataset)

        with self.db_connection.create_session() as session:
            sample_group = session.query(SampleGroup).all()
            all_colors = [s.sampleGroupBaseColor for s in sample_group]
            reordered_color_dict = defaultdict()
            n = len(all_colors)
            for i in range(n-1):
                reordered_color_dict[all_colors[i]] = all_colors[i+1]
            self.sample_attr = session.query(DatasetSeqRunSample, DatasetSeqRun, dbSample, SampleAttribute, SampleGroup).filter(dbSample.sampleID == DatasetSeqRunSample.sampleID).filter(DatasetSeqRun.datasetSeqRunID == DatasetSeqRunSample.datasetSeqRunID).filter(DatasetSeqRun.datasetSeqRunDirectory == self.dataset.dataset_name).filter(SampleGroup.sampleGroupID == DatasetSeqRunSample.sampleGroupID).filter(dbSample.sampleID == SampleAttribute.sampleID).filter(SampleAttribute.sampleAttributeTypeID == 16).order_by(DatasetSeqRunSample.sampleGroupID, DatasetSeqRunSample.sampleGroupPosition).all()
            dataset_project = session.query(dbDataset, DatasetSeqRun, Project, Scientist).filter(dbDataset.datasetID == DatasetSeqRun.datasetID).filter(DatasetSeqRun.scientistID == Scientist.scientistID).filter(Scientist.scientistID == Project.scientistID).filter(Project.projectID == dbDataset.projectID).filter(DatasetSeqRun.datasetSeqRunDirectory == self.dataset.dataset_name).one()
            ic(self.sample_attr)        
            # !!! needs self defined exception
            for s in self.sample_attr:
                s_n = int(s[3].sampleAttributeValue.split('_')[-1].split('S')[-1])
                if s_n == s[0].fastqSampleNumber:
                    fullname = s[2].sampleCPCTName + '_S' + str(s[0].fastqSampleNumber)
                    try:
                        sample = Sample(self.dataset.dataset_name, s[2].sampleCPCTName, s[0].fastqSampleNumber, sample_input_dict[fullname], s[3].sampleAttributeValue, s[0].sampleGroupID, s[0].sampleGroupPosition)
                    except KeyError:
                        sample = Sample(self.dataset.dataset_name, s[2].sampleCPCTName, s[0].fastqSampleNumber, None, s[3].sampleAttributeValue, s[0].sampleGroupID, s[0].sampleGroupPosition)
                        
                    #if self.dataset.protocol == 'RNAseq':
                    #    sample.color = reordered_color_dict[s[4].sampleGroupBaseColor]
                    # define a color by is_input
                    if s[0].chIPInputSampleFlag == 'Y':
                        sample.color = all_colors[0]
                    else:
                        #print(sample.is_input)
                        sample.color = s[4].sampleGroupBaseColor
                        
                    self.sample_meta_shortname[s[3].sampleAttributeValue] = sample
                    self.sample_meta[fullname] = sample
                else:
                    print(f'Please check sample {s[2].sampleCPCTName} \'s S #! It doesn\'t match fastqSampleNumber: {s[0].fastqSampleNumber}')
                    break
            def _set_datasetRun_info(self):
                f_s = {k[0] : k[1] for k in zip(self.sample_meta, self.sample_meta_shortname)}
                s_f = {k[1] : k[0] for k in zip(self.sample_meta, self.sample_meta_shortname)}
                self.dataset_run = DatasetRun(dataset = self.dataset.dataset_name, f_s = f_s, s_f = s_f, project = dataset_project[2].projectName, scientist = dataset_project[3].scientistFirstName)
            _set_datasetRun_info(self)

    # getting sample attributes, groups, colors, cpct_name
    # also set a DatasetRun object
    def _get_sample_info(self):
        data_p = DataParser(self.dataset)
        sample_input_dict = None
        if self.dataset.protocol == 'RNAseq':
            sample_input_dict = defaultdict() 
        else:
            sample_input_dict = data_p.parse_sample_input(self.dataset.dataset_name)

        with self.db_connection.create_session() as session:
            sample_group = session.query(SampleGroup).all()
            all_colors = [s.sampleGroupBaseColor for s in sample_group]
            reordered_color_dict = defaultdict()
            n = len(all_colors)
            for i in range(n-1):
                reordered_color_dict[all_colors[i]] = all_colors[i+1]
            self.sample_attr = session.query(DatasetSeqRunSample, DatasetSeqRun, dbSample, SampleAttribute, SampleGroup).filter(dbSample.sampleID == DatasetSeqRunSample.sampleID).filter(DatasetSeqRun.datasetSeqRunID == DatasetSeqRunSample.datasetSeqRunID).filter(DatasetSeqRun.datasetSeqRunDirectory == self.dataset.dataset_name).filter(SampleGroup.sampleGroupID == DatasetSeqRunSample.sampleGroupID).filter(dbSample.sampleID == SampleAttribute.sampleID).filter(SampleAttribute.sampleAttributeTypeID == 16).order_by(DatasetSeqRunSample.sampleGroupID, DatasetSeqRunSample.sampleGroupPosition).all()
            dataset_project = session.query(dbDataset, DatasetSeqRun, Project, Scientist).filter(dbDataset.datasetID == DatasetSeqRun.datasetID).filter(DatasetSeqRun.scientistID == Scientist.scientistID).filter(Scientist.scientistID == Project.scientistID).filter(Project.projectID == dbDataset.projectID).filter(DatasetSeqRun.datasetSeqRunDirectory == self.dataset.dataset_name).one()
            ic(self.sample_attr)        
            # !!! needs self defined exception
            for s in self.sample_attr:
                s_n = int(s[3].sampleAttributeValue.split('_')[-1].split('S')[-1])
                if s_n == s[0].fastqSampleNumber:
                    fullname = s[2].sampleCPCTName + '_S' + str(s[0].fastqSampleNumber)
                    try:
                        sample = Sample(self.dataset.dataset_name, s[2].sampleCPCTName, s[0].fastqSampleNumber, sample_input_dict[fullname], s[3].sampleAttributeValue, s[0].sampleGroupID, s[0].sampleGroupPosition)
                    except KeyError:
                        sample = Sample(self.dataset.dataset_name, s[2].sampleCPCTName, s[0].fastqSampleNumber, None, s[3].sampleAttributeValue, s[0].sampleGroupID, s[0].sampleGroupPosition)
                        
                    #if self.dataset.protocol == 'RNAseq':
                    #    sample.color = reordered_color_dict[s[4].sampleGroupBaseColor]
                    # define a color by is_input
                    if s[0].chIPInputSampleFlag == 'Y':
                        sample.color = all_colors[0]
                    else:
                        #print(sample.is_input)
                        sample.color = s[4].sampleGroupBaseColor
                        
                    self.sample_meta_shortname[s[3].sampleAttributeValue] = sample
                    self.sample_meta[fullname] = sample
                else:
                    print(f'Please check sample {s[2].sampleCPCTName} \'s S #! It doesn\'t match fastqSampleNumber: {s[0].fastqSampleNumber}')
                    break
            def _set_datasetRun_info(self):
                f_s = {k[0] : k[1] for k in zip(self.sample_meta, self.sample_meta_shortname)}
                s_f = {k[1] : k[0] for k in zip(self.sample_meta, self.sample_meta_shortname)}
                self.dataset_run = DatasetRun(dataset = self.dataset.dataset_name, f_s = f_s, s_f = s_f, project = dataset_project[2].projectName, scientist = dataset_project[3].scientistFirstName)
            _set_datasetRun_info(self)

                
        
# An instance of any data needs to be parsed from files or file system(file path)
class DataParser:
    def __init__(self, dataset: Dataset):
        self.d: Directory = None # flagstat directories info stored in a Directory object
        self.dataset: Dataset = dataset
        #self.tss_csv = tss_csv
        self.strategy: Strategy = None
        self.tss_data: pd.DataFrame = None
        print(f'This is the class for parsing all the data and metadata!')
        self.sample_meta: Dict = defaultdict()
        
    def parse_dirs(self, tss_all: str = None):
        print(f'self.dataset.protocol : === {self.dataset.protocol}')
        dirs = DirBuilder(self.dataset.dataset_name)
        # for spikein
        dirs_sp = DirBuilder(self.dataset.dataset_name)
        self.d = Directory()
        d = None
        if self.dataset.spikein_ref_genome:
            match self.dataset.protocol:
                case 'ChIPseq':
                    print(f'parse_dirs =================== chipseq')
                    dirs.build_chipseq_dirs(self.dataset.ref_genome)
                    dirs_sp.build_chipseq_spikein_dirs(self.dataset.spikein_ref_genome)
                    d = Directory()
                    d.tss_all_samples_path = tss_all
                    d.sorted_bam_flagstat_dir = dirs.sorted_bam_flagstat_dir
                    d.rmdup_bam_flagstat_dir = dirs.rmdup_bam_flagstat_dir
                    d.sorted_bam_flagstat_dir_spikein = dirs_sp.sorted_bam_flagstat_dir
                    d.rmdup_bam_flagstat_dir_spikein = dirs_sp.rmdup_bam_flagstat_dir
                    d.bedtools_coverage_dataset_rmdup_dir = dirs.bedtools_coverage_dataset_rmdup_dir
                    # tss_all_samples_path default has been set as the regular tss_all_samples.csv, unless specify your customized csv path
                    if not d.tss_all_samples_path:
                        d.tss_all_samples_path = dirs.bedtools_coverage_dataset_rmdup_dir + '/tss_all_samples.csv'
                    else:
                        d.tss_all_samples_path = tss_all
                case 'RNAseq':
                    disr.build_rnaseq_dirs(self.dataset.ref_genome)
                    d = Directory()
                    d.tss_all_samples_path = tss_all
                    d.sorted_bam_flagstat_dir = dirs.sorted_bam_flagstat_dir
                    d.rmdup_bam_flagstat_dir = dirs.rmdup_bam_flagstat_dir
                    if not d.tss_all_samples_path:
                        d.tss_all_samples_path = dirs.bedtools_coverage_dataset_rmdup_dir + '/tss_all_samples.csv'
                    else:
                        d.tss_all_samples_path = tss_all
            #d.strategy_path = 'metaconfig/bedtool_summary_config.yaml'
            self.d = d
        else:
            match self.dataset.protocol:
                case "ChIPseq":
                    print(f'parse_dirs =================== chipseq')
                    dirs.build_chipseq_dirs(self.dataset.ref_genome)
                    d = Directory()
                    d.tss_all_samples_path = tss_all
                    d.sorted_bam_flagstat_dir = dirs.sorted_bam_flagstat_dir
                    d.rmdup_bam_flagstat_dir = dirs.rmdup_bam_flagstat_dir
                    d.sorted_bam_flagstat_dir_spikein = None
                    d.rmdup_bam_flagstat_dir_spikein = None
                    d.bedtools_coverage_dataset_rmdup_dir = dirs.bedtools_coverage_dataset_rmdup_dir
                    if not d.tss_all_samples_path:
                        d.tss_all_samples_path = dirs.bedtools_coverage_dataset_rmdup_dir + '/tss_all_samples.csv'
                    else:
                        d.tss_all_samples_path = tss_all
                case "RNAseq":
                    dirs.build_rnaseq_dirs(self.dataset.ref_genome)
                    d = Directory()
                    d.tss_all_samples_path = tss_all
                    d.sorted_bam_flagstat_dir = dirs.sorted_bam_flagstat_dir
                    d.rmdup_bam_flagstat_dir = dirs.rmdup_bam_flagstat_dir
                    if not d.tss_all_samples_path:
                        d.tss_all_samples_path = dirs.bedtools_coverage_dataset_rmdup_dir + '/tss_all_samples.csv'
                    else:
                        d.tss_all_samples_path = tss_all
            #d.strategy_path = 'metaconfig/bedtool_summary_config.yaml'
            self.d = d


        
        #return d 

    # adding flagstat directories to Sample object
    def add_dirs(self, sample_meta: Dict):
        for k,v in sample_meta.items():
            sample_meta[k].sorted_flagstat = self.d.sorted_bam_flagstat_dir + sample_meta[k].sample_cpct_name + '_S' + str(sample_meta[k].s_n) + '_' + self.dataset.ref_genome + '_sorted_readgps_samtools_flagstat.txt'
            sample_meta[k].rmdup_flagstat = self.d.rmdup_bam_flagstat_dir + sample_meta[k].sample_cpct_name + '_S' + str(sample_meta[k].s_n) + '_' + self.dataset.ref_genome + '_sorted_readgps_rmdup_samtools_flagstat.txt'
            if self.dataset.spikein_ref_genome:
                sample_meta[k].sorted_flagstat_spikein = self.d.sorted_bam_flagstat_dir_spikein + sample_meta[k].sample_cpct_name + '_S' + str(sample_meta[k].s_n) + '_' + self.dataset.spikein_ref_genome + '_sorted_readgps_samtools_flagstat.txt'
                sample_meta[k].rmdup_flagstat_spikein = self.d.rmdup_bam_flagstat_dir_spikein + sample_meta[k].sample_cpct_name + '_S' + str(sample_meta[k].s_n) + '_' + self.dataset.spikein_ref_genome + '_sorted_readgps_rmdup_samtools_flagstat.txt'

   
    # adding group leader to samples
    def add_group_info(self, sample_meta: Dict) -> Dict:
        # group leader bin, a list stores all the leader key
        gl = []
        # non group leader bin
        none_gl = []
        for k,v in sample_meta.items():
            if v.group_position == 1:
                gl.append(k)
            else:
                none_gl.append(k)
        
        for g in gl:
            for n_g in none_gl:
                if sample_meta[g].group_id == sample_meta[n_g].group_id:
                    sample_meta[n_g].group_leader = g
                sample_meta[n_g]._update_group()
            sample_meta[g]._update_group()
        return sample_meta
        

    # parse flagstat and update data directly on sample object
    def _parse_flagstat(self, sample: Sample, flyreads: bool):
        if flyreads:
            try:
                with open(sample.sorted_flagstat_spikein) as f:
                    for line in f:
                        if 'properly paired' in line:
                            sample.prop_paired_sorted_sp = int(line.split(' ')[0])
                with open(sample.rmdup_flagstat_spikein) as f:
                    for line in f:
                        if 'properly paired' in line:
                            sample.prop_paired_rmdup_sp = int(line.split(' ')[0])
            except FileNotFoundError:
                return None

        try:
            with open(sample.sorted_flagstat) as f:
                for line in f:
                    if 'paired in sequencing' in line:
                        sample.seq_paired_sorted =  int(line.split(' ')[0])
                    if 'properly paired' in line:
                        sample.prop_paired_sorted = int(line.split(' ')[0]) 
            with open(sample.rmdup_flagstat) as f:
                for line in f:
                    if 'paired in sequencing' in line:
                        sample.seq_paired_rmdup =  int(line.split(' ')[0])
                    if 'properly paired' in line:
                        sample.prop_paired_rmdup = int(line.split(' ')[0]) 
        except FileNotFoundError:
            return None

        return sample
 


    def parse_strategy(self) -> Dict:
        y = YamlReader(self.d.strategy_path)
        for i in y.read_yaml():
            self.strategy = i[self.dataset.dataset_name]
        return self.strategy

    # parsing data from sample_meta.yaml
    def parse_sample_meta(self) -> Dict:
        sm = YamlReader('metaconfig/sample_meta.yaml')
        sm_yml = defaultdict()
        try:
            for i in sm.read_yaml():
                sm_yml = i[self.dataset.dataset_name]
        except KeyError:
            print(f'{self.dataset.dataset_name} isn\'t in sample_meta.yaml!')
            return 0
        for k,v in sm_yml.items():
            s = Sample(sample_cpct_name = v['sample_cpct_name'], s_n = v['s_n'])
            s.input_sample = v['input_sample']
            s.is_input = v['is_input']
            s.short_name = v['short_name']
            s.group_id = v['group_id']
            s.group_position = v['group_position']
            s.group_leader = v['group_leader']
            s.is_group_leader = v['is_group_leader']
            s.color = v['color']
            s.seq_paired_sorted = v['seq_paired_sorted']
            s.prop_paired_sorted = v['prop_paired_sorted']
            s.seq_paired_rmdup = v['seq_paired_rmdup']
            s.prop_paired_rmdup = v['prop_paired_rmdup']
            s.prop_paired_sorted_sp = v['prop_paired_sorted_sp']
            s.prop_paired_rmdup_sp = v['prop_paired_rmdup_sp']
            self.sample_meta[k] = s

        #print(self.sample_meta)
        return self.sample_meta
 
    
        
    def parse_samples(self, dataset: str):
        sample = SampleParser()
        sample.get_sample_list(dataset)
        return sample.sample_list

    # A {input:sample}
    def parse_sample_input(self, dataset: str):
        sample = SampleParser()
        sample.get_sample_input_dict(dataset)
        return sample.sample_input_dict

    def parse_tss(self, col_names) -> pd.DataFrame:
        tss = pd.read_csv(self.d.tss_all_samples_path)
        info_cols = ['Interval', 'Gene', 'Chr', 'Start', 'End']
        tss['Gene'] = [x[:-(len(x.split('-')[-1])+1)] for x in tss['Interval']]
        new_order = info_cols + col_names
        tss = tss[new_order]
        return tss
        

class DataParserCustom(DataParser):
    def __init__(self, dataset: Dataset):
        super().__init__(dataset)
        '''
        self.d: Directory = None # flagstat directories info stored in a Directory object
        self.dataset: Dataset = dataset
        #self.tss_csv = tss_csv
        self.strategy: Strategy = None
        self.tss_data: pd.DataFrame = None
        print(f'This is the class for parsing all the data and metadata!')
        self.sample_meta: Dict = defaultdict()
        '''
 
    def parse_dirs(self, tss_all: str = None):
        print(f'self.dataset.protocol : === {self.dataset.protocol}')
        self.d = Directory()
        self.d.tss_all_samples_path = tss_all
        self.d.strategy_path = 'metaconfig/bedtool_summary_config.yaml'
        d = None


class DataConfig:
    # read strategy yaml and UPDATE sample object's group related info & type(input) & color
    def __init__(self, s_d: Dict, strategy: Strategy): # a dictionary of sample objects {shortname: sample obj}
        self.s_d = s_d
        
    def update_sample(self):
        print(f'placeholder')
        
class DataWriter:
    def __init__(self, sample_meta, dataset):
        self.db_connection: DBconnector = DBconnector()
        self.attribute_type = defaultdict() # SampleAttributeType table to Dict
        self.sample_attributes = defaultdict() # final results
        self.unique_samples = None
        self.sample_info = defaultdict()

        self.sample_meta: Dict = sample_meta
        self.dataset: Dataset = dataset

    def write_attribute(self):
        with self.db_connection.create_session() as session:
            self.unique_samples = session.query(SampleAttribute.sampleID).distinct().all()
            self.attributes = session.query(SampleAttribute).order_by(SampleAttribute.sampleID).all()


        #with self.db_connection.create_session() as session:
            sample_info = session.query(dbSample, DatasetSeqRunSample).filter(dbSample.sampleID == DatasetSeqRunSample.sampleID).all()
            print(len(sample_info))
            test = len(set([s[0].sampleID for s in sample_info]))
            print(test)
            for s in sample_info:
                self.sample_info[s[0].sampleID] = s[0].sampleCPCTName + '_S' + str(s[1].fastqSampleNumber)

            
        #with self.db_connection.create_session() as session:
            attribute_type = session.query(SampleAttributeType).all()
            self.attribute_type = {a.sampleAttributeTypeID:a.sampleAttributeTypeName for a in attribute_type}

            self.sample_attributes = {self.sample_info[s[0]]:defaultdict() for s in self.unique_samples}
                
            for a in self.attributes:
                try:
                    self.sample_attributes[self.sample_info[a.sampleID]][self.attribute_type[a.sampleAttributeTypeID]] = a.sampleAttributeValue
                except KeyError:
                    pass
                    print(f'check sample id : {a.sampleAttributeTypeID}')

            with open('metaconfig/attributes.yaml', 'a') as f:
                for k,v in self.sample_attributes.items():
                    f.write(k + ':\n')
                    for ki,vi in v.items():
                        f.write('  ' + ki + ':\n')
                        f.write('    \"' + vi + '\"\n')
            print(self.sample_attributes)
            
    def write_strategy(self, config_path):
        group_samples = defaultdict()
        strategy_yml = YamlReader(config_path)
        # create group_samples
        for k,s in self.sample_meta.items():
            if not s.is_group_leader:
                if group_samples.get(s.group_leader):
                    group_samples[s.group_leader].append(k)
                else:
                    group_samples[s.group_leader] = [] # if not initialized group leader, do it
                    group_samples[s.group_leader].append(k)
        
        # writing strategy file
        for sy in strategy_yml.read_yaml():
            try:
                if sy[self.dataset.dataset_name]:
                    print(f'{self.dataset.dataset_name} exists in the file already!!!')
                    return 0
            except KeyError:
                print("good to go")

        with open(config_path, 'a') as f:
            f.write(self.dataset.dataset_name + ':\n')
            for k,s in self.sample_meta.items():
                if group_samples.get(k) or s.is_group_leader:
                    f.write('  ' + s.cpct_full_name + ':\n')
                    f.write('    type:\n')
                    if s.is_input:
                        f.write('        \"input\"\n')  
                    else:
                        f.write('        \"sample\"\n')  
                    f.write('    group_samples:\n')
                    try:
                        for gs in group_samples[k]:
                            f.write('        ' + '- \"' + gs + '\"\n')
                    except KeyError:
                        f.write('        \"0\"\n')
                    f.write('    color:\n')
                    f.write('        \"' + s.color + '\"\n')
                
            f.write('\n\n')
                    
    
    def write_sample_meta(self, extra:bool = False, extra_samples:Dict = None):
        color_conf = YamlReader('metaconfig/bedtool_summary_color_config.yaml')
        group_samples = defaultdict()
        color_dict = defaultdict() # assign color to each sample(all)
        # create group_samples
        for k,s in self.sample_meta.items():
            if not s.is_group_leader:
                if group_samples.get(s.group_leader):
                    group_samples[s.group_leader].append(k)
                else:
                    group_samples[s.group_leader] = [] # if not initialized group leader, do it
                    group_samples[s.group_leader].append(k)
        
        # create color_dict            
        for k,s in self.sample_meta.items(): 
            if not color_dict.get(k):
                if group_samples.get(k):
                    n = len(group_samples[k])
                    for c in color_conf.read_yaml():
                        color_dict[k] = c[s.color]['base']
                    sub_ix = color_picking(n)
                    for i in range(n):
                        for c in color_conf.read_yaml():
                            color_dict[group_samples[k][i]] = c[s.color]['sub'][sub_ix[i]]
                else:
                    for c in color_conf.read_yaml():
                        color_dict[k] = c[s.color]['base']

        sample_meta_check = YamlReader('metaconfig/sample_meta.yaml')
        for smc in sample_meta_check.read_yaml():
            try:
                if smc.get(self.dataset.dataset_name):
                    print(f'{self.dataset.dataset_name} exists!')
                    return 0
            except KeyError:
                print(f'Writing {self.dataset.dataset_name} sample meta!')

        with open('metaconfig/sample_meta.yaml', 'a') as f:
            f.write(self.dataset.dataset_name + ':\n')
            for k,v in self.sample_meta.items():
                f.write('  ' + k + ':\n')
                f.write('    sample_cpct_name:\n' + '      \"' + v.sample_cpct_name + '\"\n')
                f.write('    s_n:\n' + '      ' + str(v.s_n) + '\n')
                try:
                    f.write('    input_sample:\n' + '      \"' + v.input_sample + '\"\n')
                except TypeError:
                    f.write('    input_sample:\n' + '      null\n')
                f.write('    is_input:\n' + '      ' + str(v.is_input).lower() + '\n')
                f.write('    short_name:\n' + '      \"' + v.short_name + '\"\n')
                f.write('    group_id:\n' + '      ' + str(v.group_id) + '\n')
                f.write('    group_position:\n' + '      ' + str(v.group_position) + '\n')
                try:
                    f.write('    group_leader:\n' + '      \"' + v.group_leader + '\"\n')
                except TypeError:
                    f.write('    group_leader:\n' + '      null\n')
                f.write('    color:\n' + '      \"' + color_dict[k] + '\"\n')
                f.write('    is_group_leader:\n' + '      ' + str(v.is_group_leader).lower() + '\n')
                f.write('    seq_paired_sorted:\n' + '      ' + str(v.seq_paired_sorted) + '\n')
                f.write('    prop_paired_sorted:\n' + '      ' + str(v.prop_paired_sorted) + '\n')
                f.write('    seq_paired_rmdup:\n' + '      ' + str(v.seq_paired_rmdup) + '\n')
                f.write('    prop_paired_rmdup:\n' + '      ' + str(v.prop_paired_rmdup) + '\n')
                if v.prop_paired_sorted_sp:
                    f.write('    prop_paired_sorted_sp:\n' + '      ' + str(v.prop_paired_sorted_sp) + '\n')
                else:
                    f.write('    prop_paired_sorted_sp:\n' + '      null\n')
                if  v.prop_paired_rmdup_sp:
                    f.write('    prop_paired_rmdup_sp:\n' + '      ' + str(v.prop_paired_rmdup_sp) + '\n')
                else:
                    f.write('    prop_paired_rmdup_sp:\n' + '      null\n')
            f.write('\n\n')
                
            


        

# a dataframe class stores the results' dataframe
class DF:
    # passing a sample object dictionary "sample_meta"
    def __init__(self, sample_meta, dataset):
        self.sample_reads_df = None
        self.raw_counts_df = None
        self.rpk_count_ratios_df = None
        self.rpk_count_ratios_per_gene_df = None
        self.sample_meta: Dict = sample_meta
        self.dataset: Dataset = dataset
        # strategy config, dict of strategies
        self.strategy: Dict = None

    def add_color(self, cell_list, sample_names):
        color_conf = YamlReader('metaconfig/bedtool_summary_color_config.yaml')
        color_dict = defaultdict()
        color_results = defaultdict() 
            
        for k,s in self.sample_meta.items():
            if color_dict.get(s.color):
                color_dict[s.color] += 1
            else:
                color_dict[s.color] = 1
        
        color_list = []
        if len(sample_names) > len(self.sample_meta):
            for k,v in color_dict.items():
                for c in color_conf.read_yaml():
                    color_list.append(c[k]['base'])
                    for i in color_picking(v):
                        color_list.append(c[k]['sub'][i])
                        if self.dataset.flyreads:
                            for j in range(8):
                                color_list.append(c[k]['sub'][i])
                        if not self.dataset.flyreads:
                            for j in range(5):
                                color_list.append(c[k]['sub'][i])
            color_results = {c[0]: [c[1], c[2]] for c in zip(cell_list, color_list, sample_names)}
        else:

            for k,v in color_dict.items():
                for c in color_conf.read_yaml():
                    color_list.append(c[k]['base'])
                    for i in color_picking(v):
                        color_list.append(c[k]['sub'][i])
            color_results = {c[0]: [c[1], c[2]] for c in zip(cell_list, color_list, sample_names)}
        return color_results
        
        
            
        
    # A unitility tool sets a dict {'r1': 1, 'r2': 2 ...} for mapping the rows with r1, r2 names 
    # to 'rmdup % total', etc
    def set_rowname_dict(self):
        return {'r'+str(r): r for r in range(1, 36)}

    def read_txt_file(self, filepath):
        with open(filepath) as f:
            for line in f:
                yield line

    def read_csv_to_df(self, filepath):
        return pd.read_csv(filepath)

    def get_sample_reads(self):
        print(f'form a sample reads dictionary(1st tab)')


    ################## row reference ######################
    '''
    r1 = 'rmdup paired reads (1)'
    r2 = 'rmdup % total'
    r3 = 'rmdup paired read ratio (2)'
    r4 = 'rmdup ' + ref_name + ' aligned reads (1)'
    r5 = 'rmdup ' + ref_name + ' aligned read ratio (4)'
    r6 = 'rmdup dm6 aligned reads (1)'
    r7 = 'rmdup all aligned reads (3)'
    r8 = 'rmdup % aligned'
    r9 = 'rmdup spikein reads per million (5)'
    r10 = 'rmdup spikein rpm ratio (5b)'
    r11 = 'rmdup input normalization factor (6a)'
    r12 = 'rmdup input-normalized ' + ref_name + ' reads (7a)'
    r13 = 'rmdup input-normalized read ratio (8a)'
    r14 = 'rmdup read-adj input-norm read ratio (9a)'
    r15 = 'rmdup read-adj spikein read ratio (9b)'
    r16 = 'rmdup final norm factor (input-norm) (10a)'
    r17 = 'rmdup expected ' + ref_name + ' aligned reads (input-norm) (11a)'
    r18 = 'rmdup expected ' + ref_name + ' aligned reads (spikein rpm) (11b)'
    r19 = 'dupseq paired reads (1)'
    r20 = 'dupseq paired read ratio (2)'
    r21 = 'dupseq ' + ref_name + ' aligned reads (1)'
    r22 = 'dupseq ' + ref_name + ' aligned read ratio (4)'
    r23 = 'dupseq dm6 aligned reads (1)'
    r24 = 'dupseq all aligned reads (3)'
    r25 = 'dupseq % aligned'
    r26 = 'dupseq spikein reads per million (5)'
    r27 = 'dupseq spikein rpm ratio (5b)'
    r28 = 'dupseq input normalization factor (6a)'
    r29 = 'dupseq input-normalized ' + ref_name + ' reads (7a)'
    r30 = 'dupseq input-normalized read ratio (8a)'
    r31 = 'dupseq read-adj input-norm read ratio (9a)'
    r32 = 'dupseq read-adj spikein read ratio (9b)'
    r33 = 'dupseq final norm factor (input-norm) (10a)'
    r34 = 'dupseq expected ' + ref_name + ' aligned reads (input-norm) (11a)'
    r35 = 'dupseq expected ' + ref_name + ' aligned reads (spikein rpm) (11b)'
    '''

    #######################################################

    # passing a dictionary of sample objects, aka sample_meta
    # r# + 1 is the cell# in spreadsheets
    def calculate_sample_stats(self, samples: Dict):
        s_dict = defaultdict() # Dict for storing all the sample objs
        # assign the non-calculation data 
        for k,s in samples.items():
            ss = SampleStat(s = s)

            ss.r1 = s.seq_paired_rmdup
            ss.r4 = s.prop_paired_rmdup
            ss.r19 = s.seq_paired_sorted
            ss.r21 = s.prop_paired_sorted
            ss.r6 = s.prop_paired_rmdup_sp
            ss.r23 = s.prop_paired_sorted_sp
            s_dict[k] = ss

        # looping through NEW s_dict, so it can calculate across different samples
        for k,ss in s_dict.items():
            #print(f' Leader : {ss.s.group_leader}')
        
            ss.r2 = ss.r1 / ss.r19 * 100
            if not ss.s.is_group_leader:
                ss.r3 = ss.r1 / s_dict[ss.s.group_leader].r1
                ss.r5 = ss.r4 / s_dict[ss.s.group_leader].r4
                ss.r20 = ss.r19 / s_dict[ss.s.group_leader].r19
                ss.r22 = ss.r21 / s_dict[ss.s.group_leader].r21
            else:
                ss.r3 = None
                ss.r5 = None
                ss.r20 = None
                ss.r22 = None

            if ss.r6 and ss.r23:
                ss.r7 = ss.r4 + ss.r6
                ss.r8 = ss.r7 / ss.r1 * 100
                ss.r9 = ss.r6 / (ss.r7 / 1000000)
                if not ss.s.is_group_leader:
                    ss.r10 = ss.r9 / s_dict[ss.s.group_leader].r9

                if not ss.s.is_input:
                    ss.r11 = ss.r9 / s_dict[ss.s.input_sample].r9
                    ss.r12 = ss.r4 / ss.r11

                if not ss.s.is_input and not ss.s.is_group_leader:
                    ss.r13 = ss.r12 / s_dict[ss.s.group_leader].r12
                    ss.r14 = ss.r13
                    ss.r16 = ss.r14 / ss.r5

                if not ss.s.is_group_leader:
                    ss.r15 = 1 / ss.r10 
    
                if not ss.s.is_input:
                    try:
                        ss.r17 = ss.r4 * ss.r16
                        ss.r18 = ss.r4 * ss.r15
                    # for none input leader
                    except TypeError:
                        ss.r17 = ss.r4
                        ss.r18 = ss.r4

                '''
                if not ss.s.is_group_leader:
                    ss.r20 = ss.r19 / s_dict[ss.s.group_leader].r19
                    ss.r22 = ss.r21 / s_dict[ss.s.group_leader].r21
                '''
                    
                ss.r24 = ss.r21 + ss.r23
                ss.r25 = ss.r24 / ss.r19 * 100
                ss.r26 = ss.r23 / (ss.r24 / 1000000)

                if not ss.s.is_group_leader:
                    ss.r27 = ss.r26 / s_dict[ss.s.group_leader].r26

                if not ss.s.is_input:
                    ss.r28 = ss.r26 / s_dict[ss.s.input_sample].r26
                    ss.r29 = ss.r21 / ss.r28

                if not ss.s.is_input and not ss.s.is_group_leader:
                    ss.r30 = ss.r29 / s_dict[ss.s.group_leader].r29
                    ss.r31 = ss.r30
                    ss.r33 = ss.r31 / ss.r22

                if not ss.s.is_group_leader:
                    ss.r32 = 1 / ss.r27

                if not ss.s.is_input and not ss.s.is_group_leader:
                    ss.r34 = ss.r21 * ss.r33
                    ss.r35 = ss.r21 * ss.r32
                if not ss.s.is_input and ss.s.is_group_leader:
                    ss.r34 = ss.r21
                    ss.r35 = ss.r21
                
                #s_dict[k] = ss
            else:
                ss.r7 = ss.r4
                ss.r24 = ss.r21
                ss.r8 = ss.r7 / ss.r1 * 100
                ss.r25 = ss.r24 / ss.r19 * 100

        return s_dict
            
    # passing a SampleStat object Dictionary (from calculate_sample_stats())
    def calculate_rpk_count_ratio(self, stats: Dict, tss_data: pd.DataFrame, rpk_length: int = 4000, method = 'a'):
        # {sample: SampleRPK object}
        rpk_samples = defaultdict()
        n_row = len(tss_data['Interval'])
        info_colnames = ['Interval', 'Gene', 'Chr', 'Start', 'End']
        # result dictionary as pandas df 
        rpk_df = defaultdict()
        # rpk_length default has been set to 4000
        rpk_df['Interval'] = tss_data['Interval']
        rpk_df['Gene'] = [x[:-(len(x.split('-')[-1])+1)] for x in rpk_df['Interval']]
        rpk_df['Chr'] = tss_data['Chr']
        rpk_df['Start'] = tss_data['Start']
        rpk_df['End'] = tss_data['End']
        rpk_df['Lenth'] = np.array([rpk_length] * n_row)

        for k,ss in stats.items():
            s_rpk = SampleRPK(s=ss.s)
            s_rpk.rpk = np.array(tss_data[k]) / (rpk_df['Lenth'] / 1000)
            
            rpk_samples[k] = s_rpk

        # here r is the "s_rpk" previously
        for k,r in rpk_samples.items():
            if not r.s.is_input and not r.s.is_group_leader:
                r.rpk_ratio = np.array(r.rpk + 1) / np.array(rpk_samples[r.s.group_leader].rpk + 1)
                r.rpk_10k_floor = np.where(np.logical_or(r.rpk >= 10, rpk_samples[r.s.group_leader].rpk >= 10) == True, 'Y', 'N')
                r.rpk_readnorm = np.array(r.rpk) / stats[k].r5
                #???#
                #r.rpk_readnorm_ratio = np.array(r.rpk_readnorm + 1) / np.array(r.rpk + 1)
                #???#
                r.rpk_readnorm_ratio = np.array(r.rpk_readnorm + 1) / np.array(rpk_samples[r.s.group_leader].rpk+ 1)
                r.rpk_readnorm_10k_floor = np.where(np.logical_or(r.rpk_readnorm >= 10, rpk_samples[r.s.group_leader].rpk >= 10) == True, 'Y', 'N')
                '''
                if self.dataset.flyreads:
                    if method == 'a':
                        r.rpk_flynorm = np.array(r.rpk) * stats[k].r16
                    if method == 'b':
                        r.rpk_flynorm = np.array(r.rpk) * stats[k].r15
                    #r.rpk_flynorm_ratio = np.array(r.rpk_flynorm + 1) / np.array(r.rpk + 1)
                    #???#
                    r.rpk_flynorm_ratio = np.array(r.rpk_flynorm + 1) / np.array(rpk_samples[r.s.group_leader].rpk + 1)
                    #???#
                    r.rpk_flynorm_10k_floor = np.where(np.logical_or(r.rpk_flynorm >= 10, rpk_samples[r.s.group_leader].rpk >= 10) == True, 'Y', 'N')

                '''
                if stats[k].r16 is not None and method == 'a':
                    r.rpk_flynorm = np.array(r.rpk) * stats[k].r16
                    r.rpk_flynorm_ratio = np.array(r.rpk_flynorm + 1) / np.array(rpk_samples[r.s.group_leader].rpk + 1)
                    r.rpk_flynorm_10k_floor = np.where(np.logical_or(r.rpk_flynorm >= 10, rpk_samples[r.s.group_leader].rpk >= 10) == True, 'Y', 'N')
                if stats[k].r15 is not None and method == 'b':
                    r.rpk_flynorm = np.array(r.rpk) * stats[k].r15
                    r.rpk_flynorm_ratio = np.array(r.rpk_flynorm + 1) / np.array(rpk_samples[r.s.group_leader].rpk + 1)
                    r.rpk_flynorm_10k_floor = np.where(np.logical_or(r.rpk_flynorm >= 10, rpk_samples[r.s.group_leader].rpk >= 10) == True, 'Y', 'N')

        # To make the rpk_df just like what the spreadsheet shows
        #ic(rpk_samples)
        for k,v in rpk_samples.items():
            rpk_df[v.rpk_name] = v.rpk
            if not v.s.is_input and not v.s.is_group_leader:
                rpk_df[v.rpk_ratio_name] = v.rpk_ratio
                rpk_df[v.rpk_10k_floor_name] = v.rpk_10k_floor
                rpk_df[v.rpk_readnorm_name] = v.rpk_readnorm
                rpk_df[v.rpk_readnorm_ratio_name] = v.rpk_readnorm_ratio
                rpk_df[v.rpk_readnorm_10k_floor_name] = v.rpk_readnorm_10k_floor                
                #ic(v.rpk_flynorm)
                if v.rpk_flynorm is not None:
                #if self.dataset.flyreads:
                    rpk_df[v.rpk_flynorm_name] = v.rpk_flynorm
                    rpk_df[v.rpk_flynorm_ratio_name] = v.rpk_flynorm_ratio
                    rpk_df[v.rpk_flynorm_10k_floor_name] = v.rpk_flynorm_10k_floor

        
        self.rpk_count_ratios_df = rpk_df

        # just add all the rpk data get the max of each gene
        def rpk_ratio_per_gene(self):
            selected_ind = []
            gene_sum = np.array([0] * n_row)
            for k,s in rpk_samples.items():
                gene_sum = np.add(gene_sum, s.rpk)
            self.rpk_count_ratios_per_gene_df = gene_sum
            sorted_index = gene_sum.argsort()[::-1][:n_row]
            gene_dict = defaultdict()
            for i in sorted_index:
                if not gene_dict.get(rpk_df['Gene'][i]):
                    gene_dict[rpk_df['Gene'][i]] = i
                    selected_ind.append(i)
            
            rpk_pdf = pd.DataFrame(rpk_df)
            self.rpk_count_ratios_per_gene_df = rpk_pdf.loc[selected_ind]
            print(self.rpk_count_ratios_per_gene_df)
                

        rpk_ratio_per_gene(self)
        
        return rpk_samples
                       

def color_picking(num):
    if num == 1:
        return [7 // 2]
    return [x * (7 // num) for x in range(1, num + 1)]

def color_cell(n):
    c = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    results = []
    if n <= 26**2:
        for i in range(n):
            l1 = i // 26
            l2 = i % 26
            if l1 == 0:
                results.append(c[l2] + "1")
            else:
                results.append(c[l1 - 1] + c[l2] + "1")
    return results
        


def main():
    #dset = Dataset('ChIPseq_221118_JR', 'hg38')
    #dset = Dataset('ChIPseq_221114', 'mm10', spikein_ref_genome = 'dm6')
    dset = Dataset('ChIPseq_230809', 'hg38', spikein_ref_genome = 'dm6')
    #dset = Dataset('ChIPseq_230113', 'mm10')
    datab = Database(dset)
    datab._get_sample_info()
    datab._protocol_check()
    p = DataParser(dset)
    #p.parse_dirs('TSS_ALL_CUSTOMIZED')
    p.parse_dirs()
    #print(p.d)
    #p.parse_strategy()
    p.parse_sample_meta()
    sample_meta_with_group = p.add_group_info(datab.sample_meta)
    p.add_dirs(datab.sample_meta)
    for k,v in datab.sample_meta.items():
        datab.sample_meta[k] = p._parse_flagstat(v, dset.flyreads)
    tss = p.parse_tss(list(datab.sample_meta.keys()))
    #print(datab.sample_meta)
    #print(tss)
    #print(sample_meta_with_group)
    #print(p.strategy)
    df = DF(datab.sample_meta, dset)

    raw_stats = df.calculate_sample_stats(datab.sample_meta)
    for k,v in raw_stats.items():
        v.update_rs()
        raw_stats[k] = v
    #print(raw_stats)

    rpk = df.calculate_rpk_count_ratio(raw_stats, tss)
    #print(rpk)
    #print(df.rpk_count_ratios_per_gene_df)
    #print(df.rpk_count_ratios_df)

    sample_reads_df = defaultdict()
    t = Table(ref_name='hg38')
    sample_reads_df['index'] = [v for k,v in t.sample_reads_index.items()]
    for k,v in raw_stats.items():
        n_k = datab.sample_meta[k].short_name
        sample_reads_df[n_k] = v.r_np
    sample_reads_df_out = pd.DataFrame(sample_reads_df)

    print(sample_reads_df_out)

    tss_out = pd.DataFrame(tss)
    print(tss_out)
    '''
    print(list(df.rpk_count_ratios_df.keys()))

    # for c1=tab1,c2=tabs and c3=tab3,4
    c1 = df.add_color(color_cell(len(datab.sample_meta) + 1)[1:], list(sample_reads_df.keys())[1:])
    c2 = df.add_color(color_cell(len(datab.sample_meta) + 5)[5:], list(sample_reads_df.keys())[1:])
    c3 = df.add_color(color_cell(len(list(df.rpk_count_ratios_df.keys())[6:]) + 6)[6:], list(df.rpk_count_ratios_df.keys())[6:])
    print(c3)

    '''

    dw = DataWriter(datab.sample_meta, dset)
    #dw.write_strategy('metaconfig/bedtool_summary_config.yaml')
    #dw.write_attribute()
    dw.write_sample_meta()
    
    '''
    with pd.ExcelWriter('xlsx_test/ChIPseq_230113_mm10_all.xlsx', engine='xlsxwriter') as writer:
        sample_reads_df_out.to_excel(writer, sheet_name='sample_reads', index=False)
        workbook = writer.book
        worksheet = writer.sheets['sample_reads']
        a1_format = workbook.add_format()
        worksheet.write('A1', None, a1_format)
        digit_format = workbook.add_format({'num_format': '#,##0'})
        decimal_format = workbook.add_format({'num_format': '##0.00'})
        for k,v in c1.items():
            color_format = workbook.add_format({'bg_color': v[0]})
            worksheet.write(k, v[1], color_format)
        
        tss_out.to_excel(writer, sheet_name='raw_counts', index=False)
        workbook = writer.book
        worksheet = writer.sheets['raw_counts']
        digit_format = workbook.add_format({'num_format': '#,##0'})
        decimal_format = workbook.add_format({'num_format': '##0.00'})
        for k,v in c2.items():
            color_format = workbook.add_format({'bg_color': v[0]})
            worksheet.write(k, v[1], color_format)

        pd.DataFrame(df.rpk_count_ratios_df).to_excel(writer, sheet_name='rpk_count_ratios', index=False)
        workbook = writer.book
        worksheet = writer.sheets['rpk_count_ratios']
        digit_format = workbook.add_format({'num_format': '#,##0'})
        decimal_format = workbook.add_format({'num_format': '##0.00'})
        for k,v in c3.items():
            color_format = workbook.add_format({'bg_color': v[0]})
            worksheet.write(k, v[1], color_format)

        pd.DataFrame(df.rpk_count_ratios_per_gene_df).to_excel(writer, sheet_name='rpk_count_ratios_per_gene', index=False)
        workbook = writer.book
        worksheet = writer.sheets['rpk_count_ratios_per_gene']
        digit_format = workbook.add_format({'num_format': '#,##0'})
        decimal_format = workbook.add_format({'num_format': '##0.00'})
        for k,v in c3.items():
            color_format = workbook.add_format({'bg_color': v[0]})
            worksheet.write(k, v[1], color_format)

    '''


    #print(color_cell(len(datab.sample_meta) + 2)[1:-1])

    '''


    test_df = defaultdict()
    for k,v in raw_stats.items():
        n_k = datab.sample_meta[k].short_name
        test_df[n_k] = v.r_np
    print(test_df)
    out_df = pd.DataFrame(test_df)
    with pd.ExcelWriter('xlsx_test/chipseq_221118_jr.xlsx', engine='xlsxwriter') as writer:
        out_df.to_excel(writer, sheet_name='sample_reads', index=False)
        workbook = writer.book
        worksheet = writer.sheets['sample_reads']
    '''
    
    #print(datab.sample_meta)
    
   

if __name__ == '__main__':
    main()

