from dataclasses import dataclass, is_dataclass, make_dataclass, fields, field, asdict, astuple
from typing import Dict, List
from pprint import pprint
from collections import defaultdict
import numpy as np
import pandas as pd


@dataclass
class Directory:
    #dataset: str
    #ref_name: str
    sorted_bam_flagstat_dir: str = None
    rmdup_bam_flagstat_dir: str = None
    sorted_bam_flagstat_dir_spikein: str = None
    rmdup_bam_flagstat_dir_spikein: str = None
    bedtools_coverage_dataset_dir: str = None
    bedtools_coverage_dataset_rmdup_dir: str = None
    tss_all_samples_path: str = None # can be substitued with customized path
    strategy_path: str = None

    def __repr__(self):
        print(f'===== Directory object members =====')
        pprint(self.__dict__)
        return f'===== Directory object members ====='


@dataclass
class Dataset:
    dataset_name: str = None
    ref_genome: str = None
    reads_norm: bool = True
    flyreads: bool = field(default=False, init=False) 
    spikein_ref_genome: str = None
    protocol: str = None

    def __post_init__(self):
        if self.spikein_ref_genome:
            self.flyreads = True
        #_protocol_check_(self, dataset_name):

    def __repr__(self):
        print(f'===== Dataset object members =====')
        pprint(self.__dict__)
        return f'===== Dataset object members ====='

@dataclass
class DatasetRun:
    dataset: str = None # dataset_name
    # a dict of {full_cpct_name: shortname}
    f_s: Dict = None
    # a dict of {shortname: full_cpct_name}
    s_f: Dict = None
    project: str = None
    scientist: str = None
    date:str = None # can be substituted with Datetime type

    def __repr__(self):
        print(f'===== DatasetRun object members =====')
        pprint(self.__dict__)
        return f'===== DatasetRun object members ====='

@dataclass
class Files:
    flagstat_files_path: str = None
    tss_all_samples_path: str = None

    def __repr__(self):
        print(f'===== Files object members =====')
        pprint(self.__dict__)
        return f'===== Files object members ====='


@dataclass
class Protocol:
    protocol_name: str = None
    protocol_id: int = None

    def __repr__(self):
        print(f'===== Protocol object members =====')
        pprint(self.__dict__)
        return f'===== Protocol object members ====='

@dataclass
class Sample:
    dataset: str = None # dataset_name from Dataset
    sample_cpct_name: str = None
    s_n: int = None # S number/fastq number
    input_sample: str = None # if None means it's not an input
    is_input: bool = field(default=False, init=False)
    short_name: str = None # sampleAttributeValue from SampleAttribute table
    group_id: int = None # from SampleGroup table
    group_position: int = None
    group_leader: str = None
    is_group_leader: bool = field(default=False, init=False)
    color: str = None
    color_count: int = None # for rpk excel coloring
    cpct_full_name: str = field(default=None, init=False)

    # This section is the FULL FILE PATH for all the flagstats
    sorted_flagstat: str = None  
    rmdup_flagstat: str = None
    sorted_flagstat_spikein: str = None  
    rmdup_flagstat_spikein: str = None

    # flagstat sorted
    seq_paired_sorted: int = None
    prop_paired_sorted: int = None
    # flagstat rmdup
    seq_paired_rmdup: int = None
    prop_paired_rmdup: int = None
    # spike_in flagstat, spike_in shares "seq paired" with non-spikein data, for example: hg38 or mm10 alignment
    prop_paired_sorted_sp: int = None
    prop_paired_rmdup_sp: int = None
    
    
    
    def __post_init__(self):
        self.cpct_full_name = self.sample_cpct_name + '_S' + str(self.s_n)
        if not self.input_sample:
            self.is_input = True
        else:
            self.is_input = False

    def _update_group(self):
        if self.group_leader is None:
            self.is_group_leader = True
        else:
            self.is_group_leader = False
            
    def _update_rna_input(self):
        self.is_input = False        


    def __repr__(self):
        print(f'===== Sample object members =====')
        pprint(self.__dict__)
        return f'===== Sample object members ====='
 

# sample_reads stats
@dataclass
class SampleStat:
    # s is an object of Sample class
    s: Sample = None
    r1: float = None
    r2: float = None
    r3: float = None
    r4: float = None
    r5: float = None
    r6: float = None
    r7: float = None
    r8: float = None
    r9: float = None
    r10: float = None
    r11: float = None
    r12: float = None
    r13: float = None
    r14: float = None
    r15: float = None
    r16: float = None
    r17: float = None
    r18: float = None
    r19: float = None
    r20: float = None
    r21: float = None
    r22: float = None
    r23: float = None
    r24: float = None
    r25: float = None
    r26: float = None
    r27: float = None
    r28: float = None
    r29: float = None
    r30: float = None
    r31: float = None
    r32: float = None
    r33: float = None
    r34: float = None
    r35: float = None
    r_np: np.array = field(default=None, init=False)
    def __post_init__(self):
        self.r_np = np.array([self.r1, self.r2, self.r3, self.r4, self.r5, self.r6, self.r7, self.r8 ,self.r9 ,self.r10, self.r11, self.r12, self.r13, self.r14, self.r15, self.r16, self.r17, self.r18, self.r19, self.r20, self.r21, self.r22, self.r23, self.r24, self.r25, self.r26, self.r27, self.r28, self.r29, self.r30, self.r31, self.r32, self.r33, self.r34, self.r35])


    def __repr__(self):
        print(f'===== Sample Stats dictionary represents sample_reads(1st tab) data =====')
        pprint(self.__dict__)
        return f'===== Sample Stats dictionary represents sample_reads(1st tab) data ====='

    def update_rs(self):
        self.r_np = np.array([self.r1, self.r2, self.r3, self.r4, self.r5, self.r6, self.r7, self.r8 ,self.r9 ,self.r10, self.r11, self.r12, self.r13, self.r14, self.r15, self.r16, self.r17, self.r18, self.r19, self.r20, self.r21, self.r22, self.r23, self.r24, self.r25, self.r26, self.r27, self.r28, self.r29, self.r30, self.r31, self.r32, self.r33, self.r34, self.r35])

# each sample includes rpk (fly or no fly) data with all genes
@dataclass
class SampleRPK:
    s: Sample = None
    rpk_length: int = 4000 # default is 4000, formula = raw_counts / (rpk_length / 1000)
    rpk: np.array = None
    rpk_ratio: np.array = None # rpk / input rpk
    rpk_10k_floor: np.array = None #Boolean type
    rpk_readnorm: np.array = None
    rpk_readnorm_ratio: np.array = None
    rpk_readnorm_10k_floor: np.array = None
    rpk_flynorm: np.array = None
    rpk_flynorm_ratio: np.array = None
    rpk_flynorm_10k_floor: np.array = None

    def __post_init__(self):
        '''
        if self.s.is_input or self.s.is_group_leader:
            self.rpk_name = self.s.short_name + '_rpk'
        elif not self.s.is_input and not self.s.is_group_leader:
        '''
        self.rpk_name = self.s.short_name + '_rpk'
        self.rpk_ratio_name = self.s.short_name + '_rpk_ratio'
        self.rpk_10k_floor_name = self.s.short_name + '_rpk_10rpkfloor'
        self.rpk_readnorm_name = self.s.short_name + '_rpk_readnorm'
        self.rpk_readnorm_ratio_name = self.s.short_name + '_rpk_readnorm_ratio'
        self.rpk_readnorm_10k_floor_name = self.s.short_name + '_rpk_readnorm_10rpkfloor'
        self.rpk_flynorm_name = self.s.short_name + '_rpk_flynorm'
        self.rpk_flynorm_ratio_name = self.s.short_name + '_rpk_flynorm_ratio'
        self.rpk_flynorm_10k_floor_name = self.s.short_name + '_rpk_flynorm_10rpkfoor'
            

'''
@dataclass
class SampleAttribute:
    sample: Sample = None
    
'''
       
@dataclass
class SampleList:
    samples: list = None
    
    def __repr__(self):
        print(f'===== SampleList object members =====')
        pprint(self.__dict__)
        return f'===== SampleList object members ====='

# SampleInput is specific for {input:sample} objects
@dataclass
class SampleInput:
    sample_input_dict: Dict = None

    def __repr__(self):
        print(f'===== SampleInput object members =====')
        pprint(self.__dict__)
        return f'===== SampleInput object members ====='


@dataclass
class Stats:
    stat_dict: Dict = None # flagstat {sample: data}
    zscore_dict: Dict = None # {sample: zscore}

    def __repr__(self):
        print(f'===== Stats object members =====')
        pprint(self.__dict__)
        return f'===== Stats object members ====='

@dataclass
class Strategy:
    # all using cpct_full_name
    s: Sample = None
    name: str = None
    s_type: str = None
    is_group_sample: bool = None
    group_samples: List = None

    def __repr__(self):
        print(f'===== Strategy object members =====')
        pprint(self.__dict__)
        return f'===== Strategy object members ====='

@dataclass
class Table:
    ref_name: str
    sp_ref_name: str
    rowname_map: Dict = field(init = False, repr = True) 
    def __post_init__(self):
        self.sample_reads_index: Dict = {1: 'rmdup paired reads (1)',
                        2: 'rmdup % total',
                        3: 'rmdup paired read ratio (2)',
                        4: 'rmdup ' + self.ref_name + ' aligned reads (1)',
                        5: 'rmdup ' + self.ref_name + ' aligned read ratio (4)',
                        6: 'rmdup ' + self.sp_ref_name + ' aligned reads (1)',
                        7: 'rmdup all aligned reads (3)',
                        8: 'rmdup % aligned',
                        9: 'rmdup spikein reads per million (5)',
                        10: 'rmdup spikein rpm ratio (5b)',
                        11: 'rmdup input normalization factor (6a)',
                        12: 'rmdup input-normalized ' + self.ref_name + ' reads (7a)',
                        13: 'rmdup input-normalized read ratio (8a)',
                        14: 'rmdup input-normalized read ratio (9a)',
                        15: 'rmdup inverted spikein read ratio (9b)',
                        16: 'rmdup final norm factor (input-norm) (10a)',
                        17: 'rmdup expected ' + self.ref_name + ' aligned reads (input-norm) (11a)',
                        18: 'rmdup expected ' + self.ref_name + ' aligned reads (spikein rpm) (11b)',
                        19: 'dupseq paired reads (1)',
                        20: 'dupseq paired read ratio (2)',
                        21: 'dupseq ' + self.ref_name + ' aligned reads (1)',
                        22: 'dupseq ' + self.ref_name + ' aligned read ratio (4)',
                        23: 'dupseq ' + self.sp_ref_name + ' aligned reads (1)',
                        24: 'dupseq all aligned reads (3)',
                        25: 'dupseq % aligned',
                        26: 'dupseq spikein reads per million (5)',
                        27: 'dupseq spikein rpm ratio (5b)',
                        28: 'dupseq input normalization factor (6a)',
                        29: 'dupseq input-normalized ' + self.ref_name + ' reads (7a)',
                        30: 'dupseq input-normalized read ratio (8a)',
                        31: 'dupseq input-normalized read ratio (9a)',
                        32: 'dupseq inverted spikein read ratio (9b)',
                        33: 'dupseq final norm factor (input-norm) (10a)',
                        34: 'dupseq expected ' + self.ref_name + ' aligned reads (input-norm) (11a)',
                        35: 'dupseq expected ' + self.ref_name + ' aligned reads (spikein rpm) (11b)'}


    def __repr__(self):
        print(f'===== Table object members =====')
        pprint(self.__dict__)
        return f'===== Table object members ====='



def main():
    dirs = Directory('Chipseq_test', 'hg38')
    print(asdict(dirs))
    d = Dataset('ChIPseq_221118_JR', 'hg38', True)
    print(asdict(d))
    print(d.ref_genome)
    '''
    d.set_dataset_name('ChIPseq_221118_JR')
    d.set_ref_genome('hg38')
    d.set_spikein_attr(reads_norm = True, flyreads = True)
    '''
    print(d)
    f = Files('file path test')
    print(asdict(f))

    st = Stats()
    print(st)

    '''
    tb = Table('hg38')
    print(tb.set_rowname_dict())
    print(tb)
    '''

    np_test = SampleRPK(s = Sample(dataset='test_dataset_np'), rpk=np.array([1,2,3]))
    print(np_test.s.dataset)
    print(np_test)

if __name__ == '__main__':
    main()
