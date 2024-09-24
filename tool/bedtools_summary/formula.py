from collections import defaultdict
import pandas as pd
import numpy as np
from process import *
from data import *
from metaparser.sampleparser import SampleParser
import argparse
import sys
parser = argparse.ArgumentParser(description='Bedtool Summary Pipeline with Formulas')
#parser.add_argument('help', nargs='+')
parser.add_argument('-d', '--dataset', type=str)
parser.add_argument('-ref', '--ref_genome', type=str)
parser.add_argument('-ref_sp', '--ref_genome_spikein', type=str)
parser.add_argument('-tss', '--tss_file', type=str) # raw_counts second tab
parser.add_argument('-m', '--method', type=str)
parser.add_argument('-p', '--project', type=str)
parser.add_argument('-o', '--output', type=str)
parser.add_argument('--only_sample_reads', action=argparse.BooleanOptionalAction)
parser.add_argument('-wsm', '--write_sample_meta', type=str)
parser.add_argument('-prot', '--protocol', nargs='?', const=1, type=int)
parser.add_argument('-rna', '--rna_uni', type=str) # mark it as rna universal summary
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()
from icecream import ic

# given a number(distance) from A, return an Excel cell letter(s)
# e.g. if you want to get B, input 1, because it's 1 distance from A
# because the latter code has constrained n with 65 as the start point, so here the 'n' 
# is the distance, !!!form_cell(n) in general, the n>=65!!!
def form_cell(n):
    A = 65 # A ascii, just a note for where to start
    Z = 90 # Z ascii
    R = 26 # loop 26 letters
    result = None
    # l1, l2, l3 are the letters from right to left => l3l2l1 or l2l1 or l1
    l1 = '' 
    l2 = ''
    l3 = ''
    # make sure it starts at 'A'
    if n < 65:
        print(f'Please enter a number is >= 65')
        return False
    if n <= 90:
        return(chr(n))
    # 766 = 26*26+90
    if 90 < n <= 766:
        l1 = chr((n - 91) % 26 + 65) 
        l2 = chr((n - 91) // 26 + 65)
        result = l2+l1
    if 767 <= n <= 17667:
        l1 = chr((n - 91) % 26 + 65) 
        l2 = chr(((n - 91) // 26) % 26 + 65)
        l3 = chr((n - 676) // 676 + 65)
        result = l3+l2+l1
    return(result)
        
def write_down_sample_meta(dset, s_stats):
    dw = DataWriter(None, None)
    dw.write_downsample_meta(dset, s_stats)
    

# write sample_meta.yaml for RNAseq
def write_meta_rna(dataset, ref = 'hg38', ref_sp = None):
    dset = Dataset(dataset, ref, spikein_ref_genome = ref_sp)
    db = Database(dset)
    #db._protocol_check()
    db._get_sample_info_rna()
    #ic(db.sample_meta)
    print(f'protocol:=======> {db.dataset.protocol}')
    data_p = DataParser(dset)
    data_p.parse_dirs()
    ic(db.sample_meta)
    data_p.add_group_info(db.sample_meta)
    data_p.add_dirs(db.sample_meta)
    if ref_sp:
        for k,v in db.sample_meta.items():
            db.sample_meta[k] = data_p._parse_flagstat(v, True)
    else:
        for k,v in db.sample_meta.items():
            db.sample_meta[k] = data_p._parse_flagstat(v, False)
    ic(db.sample_meta)
    dw = DataWriter(db.sample_meta, dset)
    dw.write_rna_sample_meta()

# write sample_meta.yaml for ChIPseq
def write_meta(dataset, ref = 'hg38', ref_sp = None):
    dset = Dataset(dataset, ref, spikein_ref_genome = ref_sp)
    db = Database(dset)
    #ic(db)
    db.set_protocol()
    #ic(dset)
    db._get_sample_info()
    data_p = DataParser(dset)
    data_p.parse_dirs()
    #ic(db.sample_meta)
    data_p.add_group_info(db.sample_meta)
    data_p.add_dirs(db.sample_meta)
    if ref_sp:
        for k,v in db.sample_meta.items():
            print(k)
            db.sample_meta[k] = data_p._parse_flagstat(v, True)
    else:
        for k,v in db.sample_meta.items():
            db.sample_meta[k] = data_p._parse_flagstat(v, False)

    print(f'sorted_bam_flagstat:== {data_p.d.sorted_bam_flagstat_dir}')
    print(f'sorted_rmdup_bam_flagstat:== {data_p.d.rmdup_bam_flagstat_dir}')
    ic(db.sample_meta)
    dw = DataWriter(db.sample_meta, dset)
    dw.write_sample_meta()

# standalone version that doesn't rely on database
def make_formula_rna(dataset, flagstat_path, ref, ref_sp, output):
    s = SampleParser()
    s.get_sample_list(dataset)
    sample_names = s.sample_list
    stat_dict = defaultdict()

    t = None
    if ref_sp:
        t = Table(ref_name=ref, sp_ref_name=ref_sp)
    else:
        t = Table(ref_name=ref, sp_ref_name='dm6')
    index_col = [v for k,v in t.sample_reads_index.items()]
        
    for s in sample_names:
        stat_dict[s] = defaultdict()
        if ref_sp:
            try:
                with open(f'{flagstat_path}BAM_sorted/{dataset}_{ref_sp}/samtools_flagstat/{s}_{ref_sp}_sorted_readgps_samtools_flagstat.txt') as f:
                    for line in f:
                        if 'properly paired' in line:
                            stat_dict[s]['sorted_spikein'] = int(line.split(' ')[0])
                with open(f'{flagstat_path}BAM_rmdup/{dataset}_{ref_sp}/samtools_flagstat/{s}_{ref_sp}_sorted_readgps_rmdup_samtools_flagstat.txt') as f:
                    for line in f:
                        if 'properly paired' in line:
                            stat_dict[s]['rmdup_spikein'] = int(line.split(' ')[0])
            except FileNotFoundError:
                return None
        try:
            with open(f'{flagstat_path}BAM_sorted/{dataset}_{ref}/samtools_flagstat/{s}_{ref}_sorted_readgps_samtools_flagstat.txt') as f:
                for line in f:
                    if 'paired in sequencing' in line:
                        stat_dict[s]['seq_paired_sorted'] = int(line.split(' ')[0])
                    if 'properly paired' in line:
                        stat_dict[s]['prop_paired_sorted'] = int(line.split(' ')[0])
            with open(f'{flagstat_path}BAM_rmdup/{dataset}_{ref}/samtools_flagstat/{s}_{ref}_sorted_readgps_rmdup_samtools_flagstat.txt') as f:
                for line in f:
                    if 'paired in sequencing' in line:
                        stat_dict[s]['seq_paired_rmdup'] = int(line.split(' ')[0])
                    if 'properly paired' in line:
                        stat_dict[s]['prop_paired_rmdup'] = int(line.split(' ')[0])
        except FileNotFoundError:
            return None

    col_letter = defaultdict()
    s_c = 1 #sample_count
    for k in sample_names:
        col_letter[k] = form_cell(65 + s_c)
        s_c += 1
    print(col_letter)
    m = len(index_col)
    readcounts_df = defaultdict()
    readcounts_df['Index'] = index_col
    for s in sample_names:
        readcounts_df[s] = [None]*m
        readcounts_df[s][0] = stat_dict[s]['seq_paired_rmdup']
        readcounts_df[s][3] = stat_dict[s]['prop_paired_rmdup']
        readcounts_df[s][6] = stat_dict[s]['prop_paired_rmdup']
        readcounts_df[s][18] = stat_dict[s]['seq_paired_sorted']
        readcounts_df[s][20] = stat_dict[s]['prop_paired_sorted']
        readcounts_df[s][23] = stat_dict[s]['prop_paired_sorted']
        readcounts_df[s][1] = f'={col_letter[s]}2/{col_letter[s]}20*100'
        readcounts_df[s][7] = f'={col_letter[s]}8/{col_letter[s]}2*100'
        readcounts_df[s][24] = f'={col_letter[s]}25/{col_letter[s]}20*100'

    readcounts_df = pd.DataFrame(readcounts_df)
        
    # writing file part (an independent version for parsing flagstat from file system, rnaseq_meta.yaml and sample.yaml are required)
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        readcounts_df.to_excel(writer, sheet_name='sample_reads', index=False)
        workbook = writer.book
        worksheet = writer.sheets['sample_reads']
        digit_format = workbook.add_format({'num_format': '#,##0'})
        decimal_format = workbook.add_format({'num_format': '##0.00'})
        cell_range = 'B2:' + form_cell(66 + m) + str(37)
        worksheet.conditional_format(cell_range, {'type': 'cell', 'criteria': '>', 'value': 100, 'format': digit_format})
        worksheet.conditional_format(cell_range, {'type': 'cell', 'criteria': '<=', 'value': 100, 'format': decimal_format})
        border_line_range = 'A20:' + form_cell(66 + m) + '20'
        border_format = workbook.add_format()
        border_format.set_top(2)
        worksheet.conditional_format(border_line_range, {'type': 'cell', 'criteria': '>', 'value': 0, 'format': border_format})
    print(readcounts_df)
            
    #ic(stat_dict)

# for RNAseq summary (sample reads 1st tab only) 
def rna_sample_reads_formula(dataset, ref, ref_sp, output):
    col_start = 66 # ascii of B, from A:65 to Z:90
    dset = Dataset(dataset, ref, spikein_ref_genome = ref_sp)
    data_p = DataParser(dset)
    smd = data_p.parse_sample_meta() # sample_meta dict from sample_meta.yaml
    t = None
    if ref_sp:
        t = Table(ref_name=ref, sp_ref_name=ref_sp)
    else:
        t = Table(ref_name=ref, sp_ref_name='dm6')
    index_col = [v for k,v in t.sample_reads_index.items()]

    sample_names = [*smd]
    short_names = [smd[s].short_name for s in sample_names]
    col_letter = defaultdict()
    s_c = 1 #sample_count
    color_dict = defaultdict() # a color dict for the first row {sample full name: [cell, color]}
    for k in sample_names:
        col_letter[k] = form_cell(65 + s_c)
        s_c += 1
        color_dict[k] = [f'{col_letter[k]}1', smd[k].color, smd[k].short_name]

    ic(color_dict)
    m = len(index_col)
    readcounts_df = defaultdict()
    readcounts_df['Index'] = index_col
    for s in sample_names:
        readcounts_df[s] = [None]*m
        readcounts_df[s][0] = smd[s].seq_paired_rmdup
        readcounts_df[s][3] = smd[s].prop_paired_rmdup
        readcounts_df[s][6] = smd[s].prop_paired_rmdup
        readcounts_df[s][18] = smd[s].seq_paired_sorted
        readcounts_df[s][20] = smd[s].prop_paired_sorted
        readcounts_df[s][23] = smd[s].prop_paired_sorted
        readcounts_df[s][1] = f'={col_letter[s]}2/{col_letter[s]}20*100'
        readcounts_df[s][7] = f'={col_letter[s]}8/{col_letter[s]}2*100'
        readcounts_df[s][24] = f'={col_letter[s]}25/{col_letter[s]}20*100'

    readcounts_df = pd.DataFrame(readcounts_df)
    '''
    ic(smd)
    ic(smd['MCF7_RNA_PF9363_50nM_48h_REP1_230525_GGTAGC_S7'])
    ic(index_col)
    ic(readcounts_df)
    '''

    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        readcounts_df.to_excel(writer, sheet_name='sample_reads', index=False)
        workbook = writer.book
        worksheet = writer.sheets['sample_reads']
        for s in sample_names:
            color_format = workbook.add_format({'bg_color': color_dict[s][1]})
            worksheet.write(color_dict[s][0], color_dict[s][2], color_format)
        digit_format = workbook.add_format({'num_format': '#,##0'})
        decimal_format = workbook.add_format({'num_format': '##0.00'})
        cell_range = 'B2:' + form_cell(66 + m) + str(37)
        worksheet.conditional_format(cell_range, {'type': 'cell', 'criteria': '>', 'value': 100, 'format': digit_format})
        worksheet.conditional_format(cell_range, {'type': 'cell', 'criteria': '<=', 'value': 100, 'format': decimal_format})
        border_line_range = 'A20:' + form_cell(66 + m) + '20'
        border_format = workbook.add_format()
        border_format.set_top(2)
        worksheet.conditional_format(border_line_range, {'type': 'cell', 'criteria': '>', 'value': 0, 'format': border_format})


    
# doing the majority of the things 
def make_formula(dataset, ref = 'hg38', ref_sp = None, tss_dir = None, method = 'a', output = None, only_sample_reads = False):
    col_start = 66 # ascii of B, from A:65 to Z:90
    # a newly reformed dictionary with keys: "B,C,D,E..." 
    # the values are a list of selected attributes from sample_meta_dict
    f_d = defaultdict() 

    method_cell_n = 0
    if method == 'a':
        method_cell_n = 17
    elif method == 'b':
        method_cell_n = 16
    else:
        method_cell_n = None
        print('Please check your method input! Make sure it should be either a or b.')
        

    formula_dict = defaultdict() # the entire sample_reads tab 
   

    ######## data from mysql for writing sample_meta.yaml ########
    dset = Dataset(dataset, ref, spikein_ref_genome = ref_sp)
    print(f'dset is :== {dset}')
    db = Database(dset)
    db.set_protocol()
    data_p = DataParser(dset)
    data_p.parse_dirs()


    ######### write sample_meta ########

    smd = data_p.parse_sample_meta()

    table = None
    if dset.spikein_ref_genome:
        table = Table(ref_name = dset.ref_genome , sp_ref_name = dset.spikein_ref_genome)
    else:
        table = Table(ref_name = dset.ref_genome , sp_ref_name = 'dm6')
    
    sample_reads_first_col = {f'A{k+1}':v for k,v in table.sample_reads_index.items()}

    # raw_counts from tss
    tss = None
    if tss_dir:
        tss = pd.read_csv(tss_dir)
    else:
        try:
            print(f'data_p.d.tss_all_samples_path ===> {data_p.d.tss_all_samples_path}')
            tss = pd.read_csv(data_p.d.tss_all_samples_path)
        except FileNotFoundError:
            tss = None
            print('There is no tss file!')
    info_cols = ['Interval', 'Gene', 'Chr', 'Start', 'End']
    if tss is not None:
        tss['Gene'] = [x[:-(len(x.split('-')[-1])+1)] for x in tss['Interval']]
        new_order = info_cols + list(smd.keys())
        tss = tss[new_order]

    col_letter = defaultdict()
    s_c = 1 #sample_count
    for k in list(smd.keys()):
        col_letter[k] = form_cell(65 + s_c)
        s_c += 1

    
    formula_dict['top_row'] = defaultdict()

    
    if dset.flyreads:
        for k,v in smd.items():
            c_l = col_letter[k] # column letter for current k
            # its group leader's column letter
            g_l = None
            i_l = None
            try:
                g_l = col_letter[v.group_leader]
            except KeyError:
                g_l = None
            # its input's column letter
            try:
                i_l = col_letter[v.input_sample]
            except KeyError:
                i_l = None
            #print(f'g_l==={g_l}   i_l==={i_l}')

            formula_dict['top_row'][c_l + '1'] = [f'{v.short_name}', v.color]
            if v.prop_paired_rmdup_sp is None and v.prop_paired_sorted_sp is None:
                formula_dict[c_l] = {c_l + '2': f'{v.seq_paired_rmdup}',
                                     c_l + '5': f'{v.prop_paired_rmdup}',
                                     c_l + '20': f'{v.seq_paired_sorted}',
                                     c_l + '22': f'{v.prop_paired_sorted}'}

                if v.is_group_leader and v.is_input:
                    formula_dict[c_l][c_l + '3'] = f'={c_l}2/{c_l}20*100'
                if v.is_input and not v.is_group_leader:
                    formula_dict[c_l][c_l + '3'] = f'={c_l}2/{c_l}20*100'
                    formula_dict[c_l][c_l + '4'] = f'={c_l}2/{g_l}2'
                    formula_dict[c_l][c_l + '6'] = f'={c_l}5/{g_l}5'
                if v.is_group_leader and not v.is_input:
                    formula_dict[c_l][c_l + '3'] = f'={c_l}2/{c_l}20*100'
                if not v.is_group_leader and not v.is_input:
                    formula_dict[c_l][c_l + '3'] = f'={c_l}2/{c_l}20*100'
                    formula_dict[c_l][c_l + '4'] = f'={c_l}2/{g_l}2'
                    formula_dict[c_l][c_l + '6'] = f'={c_l}5/{g_l}5'

            else:
                formula_dict[c_l] = {c_l + '2': f'{v.seq_paired_rmdup}',
                                     c_l + '5': f'{v.prop_paired_rmdup}',
                                     c_l + '7': f'{v.prop_paired_rmdup_sp}',
                                     c_l + '20': f'{v.seq_paired_sorted}',
                                     c_l + '22': f'{v.prop_paired_sorted}',
                                     c_l + '24': f'{v.prop_paired_sorted_sp}'}
                if v.is_group_leader and v.is_input:
                    formula_dict[c_l][c_l + '3'] = f'={c_l}2/{c_l}20*100'
                    formula_dict[c_l][c_l + '8'] = f'={c_l}5 + {c_l}7'
                    formula_dict[c_l][c_l + '9'] = f'={c_l}8/{c_l}2*100'
                    formula_dict[c_l][c_l + '10'] = f'={c_l}7/({c_l}8/1000000)' 
                    formula_dict[c_l][c_l + '19'] = f'={c_l}5'
                    formula_dict[c_l][c_l + '25'] = f'={c_l}22 + {c_l}24'
                    formula_dict[c_l][c_l + '26'] = f'={c_l}25/{c_l}20*100'
                    formula_dict[c_l][c_l + '27'] = f'={c_l}24/({c_l}25/1000000)'
                    formula_dict[c_l][c_l + '36'] = f'={c_l}22'

                if v.is_input and not v.is_group_leader:
                    formula_dict[c_l][c_l + '3'] = f'={c_l}2/{c_l}20*100'
                    formula_dict[c_l][c_l + '4'] = f'={c_l}2/{g_l}2'
                    formula_dict[c_l][c_l + '6'] = f'={c_l}5/{g_l}5'
                    formula_dict[c_l][c_l + '8'] = f'={c_l}5 + {c_l}7'
                    formula_dict[c_l][c_l + '9'] = f'={c_l}8/{c_l}2*100'
                    formula_dict[c_l][c_l + '10'] = f'={c_l}7/({c_l}8/1000000)' 
                    formula_dict[c_l][c_l + '11'] = f'={c_l}10/{g_l}10'
                    formula_dict[c_l][c_l + '16'] = f'=1/{c_l}11'
                    formula_dict[c_l][c_l + '19'] = f'={c_l}5'
                    formula_dict[c_l][c_l + '21'] = f'={c_l}20/{g_l}20'
                    formula_dict[c_l][c_l + '23'] = f'={c_l}22/{g_l}22'
                    formula_dict[c_l][c_l + '25'] = f'={c_l}22 + {c_l}24'
                    formula_dict[c_l][c_l + '26'] = f'={c_l}25/{c_l}20*100'
                    formula_dict[c_l][c_l + '27'] = f'={c_l}24/({c_l}25/1000000)'
                    formula_dict[c_l][c_l + '28'] = f'={c_l}27/{g_l}27'
                    formula_dict[c_l][c_l + '33'] = f'=1/{c_l}28'
                    formula_dict[c_l][c_l + '36'] = f'={c_l}22'
                    
                # control sample
                if v.is_group_leader and not v.is_input:
                    formula_dict[c_l][c_l + '3'] = f'={c_l}2/{c_l}20*100'
                    formula_dict[c_l][c_l + '8'] = f'={c_l}5 + {c_l}7'
                    formula_dict[c_l][c_l + '9'] = f'={c_l}8/{c_l}2*100'
                    formula_dict[c_l][c_l + '10'] = f'={c_l}7/({c_l}8/1000000)' 
                    formula_dict[c_l][c_l + '12'] = f'={c_l}10/{i_l}10' 
                    formula_dict[c_l][c_l + '13'] = f'={c_l}5/{c_l}12' 
                    formula_dict[c_l][c_l + '18'] = f'={c_l}5' 
                    formula_dict[c_l][c_l + '19'] = f'={c_l}5' 
                    formula_dict[c_l][c_l + '25'] = f'={c_l}22 + {c_l}24'
                    formula_dict[c_l][c_l + '26'] = f'={c_l}25/{c_l}20*100'
                    formula_dict[c_l][c_l + '27'] = f'={c_l}24/({c_l}25/1000000)' 
                    formula_dict[c_l][c_l + '29'] = f'={c_l}27/{i_l}27' 
                    formula_dict[c_l][c_l + '30'] = f'={c_l}22/{c_l}29' 
                    formula_dict[c_l][c_l + '35'] = f'={c_l}22' 
                    formula_dict[c_l][c_l + '36'] = f'={c_l}22' 

                # regular sample
                if not v.is_group_leader and not v.is_input:
                    formula_dict[c_l][c_l + '3'] = f'={c_l}2/{c_l}20*100'
                    formula_dict[c_l][c_l + '4'] = f'={c_l}2/{g_l}2'
                    formula_dict[c_l][c_l + '6'] = f'={c_l}5/{g_l}5'
                    formula_dict[c_l][c_l + '8'] = f'={c_l}5 + {c_l}7'
                    formula_dict[c_l][c_l + '9'] = f'={c_l}8/{c_l}2*100'
                    formula_dict[c_l][c_l + '10'] = f'={c_l}7/({c_l}8/1000000)' 
                    formula_dict[c_l][c_l + '11'] = f'={c_l}10/{g_l}10' 
                    formula_dict[c_l][c_l + '12'] = f'={c_l}10/{i_l}10' 
                    formula_dict[c_l][c_l + '13'] = f'={c_l}5/{c_l}12' 
                    formula_dict[c_l][c_l + '14'] = f'={c_l}13/{g_l}13' 
                    formula_dict[c_l][c_l + '15'] = f'={c_l}13/{g_l}13' 
                    formula_dict[c_l][c_l + '16'] = f'=1/{c_l}11'
                    formula_dict[c_l][c_l + '17'] = f'={c_l}15/{c_l}6'
                    formula_dict[c_l][c_l + '18'] = f'={c_l}5*{c_l}17'
                    formula_dict[c_l][c_l + '19'] = f'={c_l}5*{c_l}16'
                    formula_dict[c_l][c_l + '21'] = f'={c_l}20/{g_l}20'
                    formula_dict[c_l][c_l + '23'] = f'={c_l}22/{g_l}22'
                    formula_dict[c_l][c_l + '25'] = f'={c_l}22 + {c_l}24'
                    formula_dict[c_l][c_l + '26'] = f'={c_l}25/{c_l}20*100'
                    formula_dict[c_l][c_l + '27'] = f'={c_l}24/({c_l}25/1000000)' 
                    formula_dict[c_l][c_l + '28'] = f'={c_l}27/{g_l}27'
                    formula_dict[c_l][c_l + '29'] = f'={c_l}27/{i_l}27' 
                    formula_dict[c_l][c_l + '30'] = f'={c_l}22/{c_l}29' 
                    formula_dict[c_l][c_l + '31'] = f'={c_l}30/{g_l}30' 
                    formula_dict[c_l][c_l + '32'] = f'={c_l}30/{g_l}30' 
                    formula_dict[c_l][c_l + '33'] = f'=1/{c_l}28'
                    formula_dict[c_l][c_l + '34'] = f'={c_l}32/{c_l}23'
                    formula_dict[c_l][c_l + '35'] = f'={c_l}22*{c_l}34'
                    formula_dict[c_l][c_l + '36'] = f'={c_l}22*{c_l}33'
    else:
         for k,v in smd.items():
            c_l = col_letter[k] # column letter for current k
            # its group leader's column letter
            g_l = None
            i_l = None
            try:
                g_l = col_letter[v.group_leader]
            except KeyError:
                g_l = None
            # its input's column letter
            try:
                i_l = col_letter[v.input_sample]
            except KeyError:
                i_l = None
            #print(f'g_l==={g_l}   i_l==={i_l}')

            formula_dict['top_row'][c_l + '1'] = [f'{v.short_name}', v.color]
            formula_dict[c_l] = {c_l + '2': f'{v.seq_paired_rmdup}',
                                 c_l + '5': f'{v.prop_paired_rmdup}',
                                 c_l + '20': f'{v.seq_paired_sorted}',
                                 c_l + '22': f'{v.prop_paired_sorted}'}
            if v.is_group_leader and v.is_input:
                formula_dict[c_l][c_l + '3'] = f'={c_l}2/{c_l}20*100'
                #formula_dict[c_l][c_l + '19'] = f'={c_l}5'
                #formula_dict[c_l][c_l + '36'] = f'={c_l}22'

            if v.is_input and not v.is_group_leader:
                formula_dict[c_l][c_l + '3'] = f'={c_l}2/{c_l}20*100'
                formula_dict[c_l][c_l + '4'] = f'={c_l}2/{g_l}2'
                formula_dict[c_l][c_l + '6'] = f'={c_l}5/{g_l}5'
                #formula_dict[c_l][c_l + '19'] = f'={c_l}5'
                formula_dict[c_l][c_l + '21'] = f'={c_l}20/{g_l}20'
                formula_dict[c_l][c_l + '23'] = f'={c_l}22/{g_l}22'
                #formula_dict[c_l][c_l + '36'] = f'={c_l}22'
                
            # control sample
            if v.is_group_leader and not v.is_input:
                formula_dict[c_l][c_l + '3'] = f'={c_l}2/{c_l}20*100'
                #formula_dict[c_l][c_l + '18'] = f'={c_l}5' 
                #formula_dict[c_l][c_l + '19'] = f'={c_l}5' 
                #formula_dict[c_l][c_l + '35'] = f'={c_l}22' 
                #formula_dict[c_l][c_l + '36'] = f'={c_l}22' 

            # regular sample
            if not v.is_group_leader and not v.is_input:
                formula_dict[c_l][c_l + '3'] = f'={c_l}2/{c_l}20*100'
                formula_dict[c_l][c_l + '4'] = f'={c_l}2/{g_l}2'
                formula_dict[c_l][c_l + '6'] = f'={c_l}5/{g_l}5'
                formula_dict[c_l][c_l + '21'] = f'={c_l}20/{g_l}20'
                formula_dict[c_l][c_l + '23'] = f'={c_l}22/{g_l}22'

    raw_counts_top_row = defaultdict()
    r_top_start = 70 # raw_counts tab starts with F
    for k,v in formula_dict['top_row'].items():
        raw_counts_top_row[form_cell(r_top_start) + str(1)] = v
        r_top_start += 1


    n = len(tss['Gene']) # placeholder
    rpk_ratio_dict = defaultdict()
    rpk_s_c = 1
    for s in [*smd]:
        rpk_ratio_dict[s] = defaultdict()
        rpk_ratio_dict[s]['name'] = []
        rpk_ratio_dict[s]['s_cell'] = form_cell(65 + rpk_s_c) # here the 'cell' matches the sample_reads tab 
        rpk_ratio_dict[s]['r_cell'] = form_cell(69 + rpk_s_c) # here the 'cell' matches the sample_reads tab 
        rpk_ratio_dict[s]['cell'] = []
        rpk_s_c += 1
        

    rpk_tab_top_row_dict = {'group_leader': '_rpk', 
                            'input': '_rpk',
                            'input_leader': '_rpk',
                            'sample_': ['_rpk', '_rpk_ratio', '_rpk_10rpkfloor', '_rpk_readnorm', '_rpk_readnorm_ratio',
                                        '_rpk_readnorm_10rpkfloor'],
                            'sample_sp': ['_rpk_flynorm', '_rpk_flynorm_ratio', '_rpk_flynorm_10rpkfloor']}

    name_cell_dict = defaultdict()

    rpk_cell_start = 71  # G
    c_c = 1
    print(f'flyreads: y?n {dset.flyreads}')
    if dset.flyreads:
        for k,v in smd.items():
            if (v.prop_paired_rmdup_sp is None):
                print(f'k================={k}')
            if v.is_group_leader or v.is_input:
                rpk_ratio_dict[k]['name'].append(v.short_name + '_rpk')
                #rpk_ratio_dict[k]['cell'][form_cell(rpk_cell_start)] = 1
                rpk_ratio_dict[k]['cell'].append(form_cell(rpk_cell_start))
                name_cell_dict[v.short_name + '_rpk'] = form_cell(rpk_cell_start)
                rpk_cell_start += 1
            else:
                if (v.prop_paired_rmdup_sp is None):
                    print(f'(k================={k})')
                    for i in rpk_tab_top_row_dict['sample_']:
                        rpk_ratio_dict[k]['name'].append(v.short_name + i)
                        #rpk_ratio_dict[k]['cell'][form_cell(rpk_cell_start)] = c_c
                        rpk_ratio_dict[k]['cell'].append(form_cell(rpk_cell_start))
                        name_cell_dict[v.short_name + i] = form_cell(rpk_cell_start)
                        rpk_cell_start += 1
                        c_c += 1
                else:
                    for i in rpk_tab_top_row_dict['sample_']:
                        rpk_ratio_dict[k]['name'].append(v.short_name + i)
                        #rpk_ratio_dict[k]['cell'][form_cell(rpk_cell_start)] = c_c
                        rpk_ratio_dict[k]['cell'].append(form_cell(rpk_cell_start))
                        name_cell_dict[v.short_name + i] = form_cell(rpk_cell_start)
                        rpk_cell_start += 1
                        c_c += 1
                    for j in rpk_tab_top_row_dict['sample_sp']:
                        rpk_ratio_dict[k]['name'].append(v.short_name + j)
                        #rpk_ratio_dict[k]['cell'][form_cell(rpk_cell_start)] = c_c
                        rpk_ratio_dict[k]['cell'].append(form_cell(rpk_cell_start))
                        name_cell_dict[v.short_name + i] = form_cell(rpk_cell_start)
                        rpk_cell_start += 1
                        c_c += 1
            c_c = 1
    else:
        for k,v in smd.items():
            if v.is_group_leader or v.is_input:
                rpk_ratio_dict[k]['name'].append(v.short_name + '_rpk')
                #rpk_ratio_dict[k]['cell'][form_cell(rpk_cell_start)] = 1
                rpk_ratio_dict[k]['cell'].append(form_cell(rpk_cell_start))
                name_cell_dict[v.short_name + '_rpk'] = form_cell(rpk_cell_start)
                rpk_cell_start += 1
            else:
                for i in rpk_tab_top_row_dict['sample_']:
                    rpk_ratio_dict[k]['name'].append(v.short_name + i)
                    #rpk_ratio_dict[k]['cell'][form_cell(rpk_cell_start)] = c_c
                    rpk_ratio_dict[k]['cell'].append(form_cell(rpk_cell_start))
                    name_cell_dict[v.short_name + i] = form_cell(rpk_cell_start)
                    rpk_cell_start += 1
                    c_c += 1                
            c_c = 1

    m = rpk_cell_start - 71 + 1 # number of columns after adding rpk stuffs
        
        
    rpk_ratio_top_row_info = {'A1': 'Interval', 'B1': 'Gene', 'C1': 'Chr', 'D1': 'Start', 'E1': 'End', 'F1': 'Length'}
    rpk_ratio_top_row_sample = defaultdict()
    for k,v in rpk_ratio_dict.items():
        c_n = len(v['cell'])
        for i in range(c_n):
            rpk_ratio_top_row_sample[v['cell'][i] + str(1)] = [v['name'][i], smd[k].color]

    formula_rpk_dict = defaultdict()
    for info in info_cols:
        formula_rpk_dict[info] = tss[info]
    formula_rpk_dict['Length'] = [4000] * n 
    for k,v in rpk_ratio_dict.items():
        for v1 in v['cell']:
            formula_rpk_dict[v1] = [''] * n
    for k,v in rpk_ratio_dict.items():
        if smd[k].is_group_leader and smd[k].is_input:
            r = v['r_cell']
            for i in range(2, n+1):
                formula_rpk_dict[[*v['cell']][0]][i-2] = f'=raw_counts!{r}{i}/(F{i}/1000)'
        if smd[k].is_input and not smd[k].is_group_leader:
            r = v['r_cell']
            for i in range(2, n+1):
                formula_rpk_dict[[*v['cell']][0]][i-2] = f'=raw_counts!{r}{i}/(F{i}/1000)'
        if not smd[k].is_input and smd[k].is_group_leader:
            r = v['r_cell']
            for i in range(2, n+1):
                formula_rpk_dict[[*v['cell']][0]][i-2] = f'=raw_counts!{r}{i}/(F{i}/1000)'
        if not smd[k].is_input and not smd[k].is_group_leader:
            i_s = smd[k].input_sample
            l_s = smd[k].group_leader
            if len(v['cell']) == 9:
                cells = [*rpk_ratio_dict[k]['cell']]
                r = rpk_ratio_dict[k]['r_cell']
                l = [*rpk_ratio_dict[l_s]['cell']][0]
                rpk = cells[0]
                rpk_ratio = cells[1]
                s_read_ratio = rpk_ratio_dict[k]['s_cell']
                rpk_readnorm = cells[3]
                rpk_flynorm = cells[6]
                ### rpk
                for i in range(2, n+1):
                    formula_rpk_dict[v['cell'][0]][i-2] = f'=raw_counts!{r}{i}/(F{i}/1000)'
                ### rpk ratio
                for i in range(2, n+1):
                    formula_rpk_dict[v['cell'][1]][i-2] = f'=({rpk}{i}+1)/({l}{i}+1)'
                ### rpk 10rpkfloor \042 represents "
                for i in range(2, n+1):
                    formula_rpk_dict[v['cell'][2]][i-2] = f"=IF({l}{i}>=10,\042Y\042,IF({rpk}{i}>=10,\042Y\042,\042N\042))"
                ### rpk_readnorm
                for i in range(2, n+1):
                    formula_rpk_dict[v['cell'][3]][i-2] = f'={rpk}{i}/sample_reads!{s_read_ratio}6'
                ### rpk_readnorm_ratio
                for i in range(2, n+1):
                    formula_rpk_dict[v['cell'][4]][i-2] = f'=({rpk_readnorm}{i}+1)/({l}{i}+1)'
                ### rpk_readnorm_10k_floor
                for i in range(2, n+1):
                    formula_rpk_dict[v['cell'][5]][i-2] = f"=IF({l}{i}>=10,\042Y\042,IF({rpk_readnorm}{i}>=10,\042Y\042,\042N\042))"
                '''
                if dset.flyreads:
                '''
                ### rpk_flynorm
                for i in range(2, n+1):
                    formula_rpk_dict[v['cell'][6]][i-2] = f'={rpk}{i}*sample_reads!{s_read_ratio}{method_cell_n}'
                ### rpk_flynorm_ratio
                for i in range(2, n+1):
                    formula_rpk_dict[v['cell'][7]][i-2] = f'=({rpk_flynorm}{i}+1)/({l}{i}+1)'
                ### rpk_flynorm_10k_floor
                for i in range(2, n+1):
                    formula_rpk_dict[v['cell'][8]][i-2] =f"=IF({l}{i}>=10,\042Y\042,IF({rpk_flynorm}{i}>=10,\042Y\042,\042N\042))"

            else:
                print('Processing data without fly spikein!')
                cells = [*rpk_ratio_dict[k]['cell']]
                r = rpk_ratio_dict[k]['r_cell']
                l = [*rpk_ratio_dict[l_s]['cell']][0]
                rpk = cells[0]
                rpk_ratio = cells[1]
                s_read_ratio = rpk_ratio_dict[k]['s_cell']
                rpk_readnorm = cells[3]
                ### rpk
                for i in range(2, n+1):
                    formula_rpk_dict[v['cell'][0]][i-2] = f'=raw_counts!{r}{i}/(F{i}/1000)'
                ### rpk ratio
                for i in range(2, n+1):
                    formula_rpk_dict[v['cell'][1]][i-2] = f'=({rpk}{i}+1)/({l}{i}+1)'
                ### rpk 10rpkfloor \042 represents "
                for i in range(2, n+1):
                    formula_rpk_dict[v['cell'][2]][i-2] = f"=IF({l}{i}>=10,\042Y\042,IF({rpk}{i}>=10,\042Y\042,\042N\042))"
                ### rpk_readnorm
                for i in range(2, n+1):
                    formula_rpk_dict[v['cell'][3]][i-2] = f'={rpk}{i}/sample_reads!{s_read_ratio}6'
                ### rpk_readnorm_ratio
                for i in range(2, n+1):
                    formula_rpk_dict[v['cell'][4]][i-2] = f'=({rpk_readnorm}{i}+1)/({l}{i}+1)'
                ### rpk_readnorm_10k_floor
                for i in range(2, n+1):
                    formula_rpk_dict[v['cell'][5]][i-2] = f"=IF({l}{i}>=10,\042Y\042,IF({rpk_readnorm}{i}>=10,\042Y\042,\042N\042))"   
    
    formula_rpk_df = pd.DataFrame(formula_rpk_dict)
    
    ############ rpk_count_ratios_per_gene ############
    rpk_df = DF(smd, dset)
    stats = rpk_df.calculate_sample_stats(smd)
    write_down_sample_meta(dset, stats)
    #ic(stats)
    
    rpk_df.calculate_rpk_count_ratio(stats, tss, method=method)
    rpk_per_gene_df = rpk_df.rpk_count_ratios_per_gene_df
    rpk_per_gene_df = rpk_per_gene_df.sort_values(by=['Chr', 'Start'])
    #print(rpk_per_gene_df)
     
    
    ############ rpk_count_ratios_per_gene ############
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        if only_sample_reads:
            workbook = writer.book
            worksheet = workbook.add_worksheet('sample_reads')
            for k,v in formula_dict.items():
                if k == 'top_row':
                    for ki,vi in v.items():
                        color_format = workbook.add_format({'bg_color': vi[1]})
                        worksheet.write(ki, vi[0], color_format)
                else:
                    for ki,vi in v.items():
                        worksheet.write_dynamic_array_formula(ki, vi)
            for k,v in sample_reads_first_col.items():
                worksheet.write(k, v)
            if method == 'a':
                color_format = workbook.add_format({'bg_color': '#ffff00'})
                worksheet.write('A17', sample_reads_first_col['A17'], color_format)
            if method == 'b':
                color_format = workbook.add_format({'bg_color': '#ffff00'})
                worksheet.write('A16', sample_reads_first_col['A16'], color_format)
            digit_format = workbook.add_format({'num_format': '#,##0'})
            decimal_format = workbook.add_format({'num_format': '##0.00'})
            cell_range = 'B2:' + form_cell(66 + len(formula_dict)) + str(37)
            worksheet.conditional_format(cell_range, {'type': 'cell', 'criteria': '>', 'value': 100, 'format': digit_format})
            worksheet.conditional_format(cell_range, {'type': 'cell', 'criteria': '<=', 'value': 100, 'format': decimal_format})
            border_line_range = 'A20:' + form_cell(66 + len(formula_dict)) + '20'
            border_format = workbook.add_format()
            border_format.set_top(2)
            worksheet.conditional_format(border_line_range, {'type': 'cell', 'criteria': '>', 'value': 0, 'format': border_format})
        else:
            workbook = writer.book
            worksheet = workbook.add_worksheet('sample_reads')
            for k,v in formula_dict.items():
                if k == 'top_row':
                    for ki,vi in v.items():
                        color_format = workbook.add_format({'bg_color': vi[1]})
                        worksheet.write(ki, vi[0], color_format)
                else:
                    for ki,vi in v.items():
                        worksheet.write_dynamic_array_formula(ki, vi)
            for k,v in sample_reads_first_col.items():
                worksheet.write(k, v)
            if method == 'a':
                color_format = workbook.add_format({'bg_color': '#ffff00'})
                worksheet.write('A17', sample_reads_first_col['A17'], color_format)
            if method == 'b':
                color_format = workbook.add_format({'bg_color': '#ffff00'})
                worksheet.write('A16', sample_reads_first_col['A16'], color_format)
            digit_format = workbook.add_format({'num_format': '#,##0'})
            decimal_format = workbook.add_format({'num_format': '##0.00'})
            cell_range = 'B2:' + form_cell(66 + len(formula_dict)) + str(37)
            worksheet.conditional_format(cell_range, {'type': 'cell', 'criteria': '>', 'value': 100, 'format': digit_format})
            worksheet.conditional_format(cell_range, {'type': 'cell', 'criteria': '<=', 'value': 100, 'format': decimal_format})
            border_line_range = 'A20:' + form_cell(66 + len(formula_dict)) + '20'
            border_format = workbook.add_format()
            border_format.set_top(2)
            worksheet.conditional_format(border_line_range, {'type': 'cell', 'criteria': '>', 'value': 0, 'format': border_format})

            tss.to_excel(writer, sheet_name='raw_counts', index=False)
            workbook = writer.book
            worksheet = writer.sheets['raw_counts']
            for k,v in raw_counts_top_row.items():
                color_format = workbook.add_format({'bg_color': v[1]})
                worksheet.write(k, v[0], color_format)
            digit_format = workbook.add_format({'num_format': '#,##0'})
            decimal_format = workbook.add_format({'num_format': '##0.00'})

            formula_rpk_df.to_excel(writer, sheet_name='rpk_count_ratios', index=False)
            workbook = writer.book
            worksheet = writer.sheets['rpk_count_ratios']
            digit_format = workbook.add_format({'num_format': '#,##0'})
            decimal_format = workbook.add_format({'num_format': '##0.00'})
            worksheet.set_column(7, 7 + m, None, decimal_format)
            for k,v in rpk_ratio_top_row_info.items():
                worksheet.write(k, v)
            for k,v in rpk_ratio_top_row_sample.items():
                color_format = workbook.add_format({'bg_color': v[1]})
                worksheet.write(k, v[0], color_format)
        
            rpk_per_gene_df.to_excel(writer, sheet_name='rpk_count_ratios_per_gene', index=False)
            workbook = writer.book
            worksheet = writer.sheets['rpk_count_ratios_per_gene']
            digit_format = workbook.add_format({'num_format': '#,##0'})
            decimal_format = workbook.add_format({'num_format': '##0.00'})
            worksheet.set_column(7, 7 + m, None, decimal_format)
            for k,v in rpk_ratio_top_row_info.items():
                worksheet.write(k, v)
            for k,v in rpk_ratio_top_row_sample.items():
                color_format = workbook.add_format({'bg_color': v[1]})
                worksheet.write(k, v[0], color_format)


def main():
    if args.write_sample_meta:
        if 'RNAseq' in args.dataset:
            write_meta_rna(dataset = args.dataset, ref = args.ref_genome, ref_sp = args.ref_genome_spikein)
        else:
            write_meta(dataset = args.dataset, ref = args.ref_genome, ref_sp = args.ref_genome_spikein)
    else:
        if args.only_sample_reads:
            if 'RNAseq' in args.dataset:
                ## rna universal summary
                #make_formula_rna(dataset = args.dataset, flagstat_path = '/mnt/storage/dept/pedonc/CPCT/projects/MCF7_NO-202204/', ref = args.ref_genome, ref_sp = args.ref_genome_spikein, output = args.output)
                #make_formula_rna(dataset = args.dataset, flagstat_path = args.project, ref = args.ref_genome, ref_sp = args.ref_genome_spikein, output = args.output)
                rna_sample_reads_formula(dataset = args.dataset, ref = args.ref_genome, ref_sp = args.ref_genome_spikein, output = args.output)
            else:
                if args.method == 'a' or args.method == 'b':
                    print('Please check if you have dm6 input for method a or b')
                make_formula(dataset = args.dataset, ref = args.ref_genome, ref_sp = args.ref_genome_spikein, tss_dir = args.tss_file, method = args.method, output = args.output, only_sample_reads = args.only_sample_reads)
        else:
            if args.tss_file:
                if args.method == 'a' or args.method == 'b':
                    print('Please check if you have dm6 input for method a or b')
                make_formula(dataset = args.dataset, ref = args.ref_genome, ref_sp = args.ref_genome_spikein, tss_dir = args.tss_file, method = args.method, output = args.output)
            else:
                if args.method == 'a' or args.method == 'b':
                    print('Please check if you have dm6 input for method a or b')
                make_formula(dataset = args.dataset, ref = args.ref_genome, ref_sp = args.ref_genome_spikein, tss_dir = None, method = args.method, output = args.output)


if __name__ == "__main__":
    main()
    

