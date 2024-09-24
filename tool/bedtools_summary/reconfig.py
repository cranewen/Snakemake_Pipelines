from collections import defaultdict
import sys, os
sys.path.insert(0, '/home/yw900/lab_pipelines')
import yaml
from metaparser.metaconf import BedtoolSummaryColorConfig
from metaparser.sampleparser import SampleParser
from utility import YamlReader
from data import *
from typing import Dict
from process import color_picking, color_cell

class Reconf:
    def __init__(self, new_dataset: str = None, datasets: list = None):
        self.new_dataset: str = new_dataset
        self.datasets: list = datasets
        self.datasets_dict: Dict = defaultdict() # parsing multiple original datasets of sample_meta {dataset: dict info}
        self.sample_meta: Dict = defaultdict() # sample meta of datasets
        self.sample_input_dict: Dict = defaultdict()

    def parse_sample_input(self):
        sample = SampleParser()
        sample.get_sample_input_dict(self.new_dataset)
        self.sample_input_dict = sample.sample_input_dict
        return sample.sample_input_dict

    def get_orig_sample_meta(self):
        sm = YamlReader('metaconfig/sample_meta.yaml')
        sm_yml = defaultdict()
        for d in self.datasets:
            try:
                for i in sm.read_yaml():
                    sm_yml = i[d]
            except KeyError:
                print(f'{d} isn\'t in sample_meta.yaml!')
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
        return self.sample_meta

    def get_new_strategy(self):
        st = YamlReader('metaconfig/bedtool_summary_config.yaml')
        new_st_yml = defaultdict()
        try:
            for i in st.read_yaml():
                new_st_yml = i[self.new_dataset]
        except KeyError:
            print(f'{self.new_dataset} isn\'t in bedtool_summary_config.yaml')
        return new_st_yml

    def rewrite_sample_meta(self, new_st):
        new_sample_meta = defaultdict()
        group_leader = []
        group_id = 1
        group_position = 1
        for k,v in new_st.items():
            if not self.sample_meta.get(k):
                s_tail = k.split('_')[-1]
                l_t = len(s_tail) # S# tail length
                l = len(k) # length of cpct full name
                for k1,v1 in self.sample_meta.items():
                    if k[:l - l_t - 1] == v1.sample_cpct_name:
                        new_sample_meta[k] = self.sample_meta[k1]
                        new_sample_meta[k].group_id = group_id
                        new_sample_meta[k].group_position = group_position
                        new_sample_meta[k].s_n = int(s_tail.split('S')[1])
                        new_sample_meta[k].cpct_full_name = new_sample_meta[k].sample_cpct_name + '_S' + str(new_sample_meta[k].s_n)
                        new_sample_meta[k].is_group_leader = True
                        new_sample_meta[k].short_name = '_'.join(new_sample_meta[k].short_name.split('_')[:-1]) + '_' + s_tail
                        new_sample_meta[k].color = v['color']
                        group_position += 1
                        try:
                            new_sample_meta[k].input_sample = self.sample_input_dict[k]
                        except KeyError:
                            pass
            else:
                new_sample_meta[k] = self.sample_meta[k]
                new_sample_meta[k].color = v['color']
                
            for v_1 in v['group_samples']:
                if not self.sample_meta.get(v_1):
                    s_tail = v_1.split('_')[-1]
                    l_t = len(s_tail) # S# tail length
                    l = len(v_1) # length of cpct full name
                    for k1,v1 in self.sample_meta.items():
                        if v_1[:l - l_t - 1] == v1.sample_cpct_name:
                            new_sample_meta[v_1] = self.sample_meta[k1]
                            new_sample_meta[v_1].group_id = group_id
                            new_sample_meta[v_1].group_position = group_position
                            new_sample_meta[v_1].s_n = int(s_tail.split('S')[1])
                            new_sample_meta[v_1].cpct_full_name = new_sample_meta[v_1].sample_cpct_name + '_S' + str(new_sample_meta[v_1].s_n)
                            new_sample_meta[v_1].group_leader = k
                            new_sample_meta[v_1].short_name = '_'.join(new_sample_meta[v_1].short_name.split('_')[:-1]) + '_' + s_tail
                            new_sample_meta[v_1].color = v['color']
                            group_position += 1
                            try:
                                new_sample_meta[v_1].input_sample = self.sample_input_dict[v_1]
                            except KeyError:
                                pass
                else:
                    new_sample_meta[v_1] = self.sample_meta[v_1]
                    new_sample_meta[v_1].color = v['color']

                
            group_position = 1
            group_id += 1

        # this section is for reorder the new_sample_meta, because later the calculate needs the inputs being calculated first,
        # otherwise, calculation could have NoneType in the process
        leader_input_bin = []
        leader_bin = []
        input_bin = []
        sample_bin = []
        for k,v in new_sample_meta.items():
            if v.is_group_leader and v.is_input:
                leader_input_bin.append(k)
            elif v.is_input:
                input_bin.append(k)
            elif v.is_group_leader:
                leader_bin.append(k)
            elif not v.is_group_leader and not v.is_input:
                sample_bin.append(k)
        new_order = leader_input_bin + leader_bin + input_bin + sample_bin
        print(new_order)
        new_sample_meta_copy = new_sample_meta.copy()
        new_sample_meta = defaultdict()
        new_sample_meta = {k:new_sample_meta_copy[k] for k in new_order}
        #print(new_sample_meta)

        color_conf = YamlReader('metaconfig/bedtool_summary_color_config.yaml')
        group_samples = defaultdict()
        color_dict = defaultdict() # assign color to each sample(all)
        for k,s in new_sample_meta.items():
            if not s.is_group_leader:
                if group_samples.get(s.group_leader):
                    group_samples[s.group_leader].append(k)
                else:
                    group_samples[s.group_leader] = [] # if not initialized group leader, do it
                    group_samples[s.group_leader].append(k)

        for k,s in new_sample_meta.items():
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

        #print(color_dict)

        sample_meta_check = YamlReader('metaconfig/sample_meta.yaml')
        for smc in sample_meta_check.read_yaml():
            try:
                if smc.get(self.new_dataset):
                    print(f'{self.new_dataset} exists!')
                    return 0
            except KeyError:
                print(f'Writing {self.new_dataset} sample meta!')
        
        with open('metaconfig/sample_meta.yaml', 'a') as f:
            f.write(self.new_dataset + ':\n')
            for k,v in new_sample_meta.items():
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

        return new_sample_meta
                            
                        

            
def main():
    reconf = Reconf('ChIPseq_230209_hs_combine_test', ['ChIPseq_230209_hs', 'ChIPseq_230127_M7'])
    reconf.get_orig_sample_meta()
    #print(len(reconf.sample_meta))
    #print(reconf.rewrite_sample_meta(reconf.get_new_strategy()))
    reconf.rewrite_sample_meta(reconf.get_new_strategy())

if __name__ == '__main__':
    main()
        
