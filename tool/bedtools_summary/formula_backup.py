from collections import defaultdict
import pandas as pd
from process import *
from data import *
#from openpyxl import Workbook

def form_cell(n):
    a = 65
    z = 90
    R = 26
    result = None
    # l1, l2, l3 are the letters from right to left => l3l2l1 or l2l1 or l1
    l1 = '' 
    l2 = ''
    l3 = ''
    if n <= 90:
        return(chr(n))
    # 767 = 26*26+90+1
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
        

def make_formula():
    col_start = 66 # ascii of B, from A:65 to Z:90
    # a newly reformed dictionary with keys: "B,C,D,E..." 
    # the values are a list of selected attributes from sample_meta_dict
    f_d = defaultdict() 

    formula_dict = defaultdict() # the entire sample_reads tab 
   
    dset = Dataset('ChIPseq_230601', 'hg38', spikein_ref_genome = 'dm6')
    #dset = Dataset('ChIPseq_230524', 'hg38', spikein_ref_genome = 'dm6')
    #dset = Dataset('ChIPseq_230809', 'hg38', spikein_ref_genome = 'dm6')
    data_p = DataParser(dset)
    smd = data_p.parse_sample_meta()
    #print(smd)

    # raw_counts from tss
    tss = pd.read_csv('xlsx_test/chipseq_230601_tss.csv')
    info_cols = ['Interval', 'Gene', 'Chr', 'Start', 'End']
    tss['Gene'] = [x[:-(len(x.split('-')[-1])+1)] for x in tss['Interval']]
    new_order = info_cols + list(smd.keys())
    tss = tss[new_order]
    #print(tss.columns)

    col_letter = defaultdict()
    s_c = 1 #sample_count
    for k in list(smd.keys()):
        col_letter[k] = form_cell(65 + s_c)
        s_c += 1

    ########################################### delete later #################################################
    '''
    f_i = 0 # first letter incremental mark
    first_l = ''
    for k in list(smd.keys()):
        if col_start <= 90:
            col_letter[k] = first_l + chr(col_start)
            col_start += 1
        else:
            first_l = chr(65 + f_i)
            col_start = 65
            f_i += 1
    '''
    ########################################### delete later #################################################


    formula_dict['top_row'] = defaultdict()

    ####################################################################
    all_dict = defaultdict() # everything {sample : {'sample_reads': {'B': ...}, 'raw_counts': {'B':...}, ...}}

    tab_cell_dict = defaultdict() # {sample : {'sample_reads': 'cellxx', 'raw_counts': 'cellxx' ...}}
    tab_cell_dict['sample_reads'] = defaultdict()
    tab_cell_dict['raw_counts'] = defaultdict()
    tab_cell_dict['rpk_count_ratios'] = defaultdict()
    tab_cell_dict['rpk_count_ratios_per_gene'] = defaultdict()
    ####################################################################

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
        print(f'g_l==={g_l}   i_l==={i_l}')

        formula_dict['top_row'][c_l + '1'] = [f'{v.short_name}', v.color]
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
            formula_dict[c_l][c_l + '16'] = f'=(1/{c_l}11)/{c_l}4'
            formula_dict[c_l][c_l + '19'] = f'={c_l}5'
            formula_dict[c_l][c_l + '21'] = f'={c_l}20/{g_l}20'
            formula_dict[c_l][c_l + '23'] = f'={c_l}22/{g_l}22'
            formula_dict[c_l][c_l + '25'] = f'={c_l}22 + {c_l}24'
            formula_dict[c_l][c_l + '26'] = f'={c_l}25/{c_l}20*100'
            formula_dict[c_l][c_l + '27'] = f'={c_l}24/({c_l}25/1000000)'
            formula_dict[c_l][c_l + '28'] = f'={c_l}27/{g_l}27'
            formula_dict[c_l][c_l + '33'] = f'=(1/{c_l}28)/{c_l}21'
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
            formula_dict[c_l][c_l + '15'] = f'={c_l}14/{c_l}4' 
            formula_dict[c_l][c_l + '16'] = f'=(1/{c_l}11)/{c_l}4'
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
            formula_dict[c_l][c_l + '32'] = f'={c_l}31/{c_l}21' 
            formula_dict[c_l][c_l + '33'] = f'=(1/{c_l}28)/{c_l}21'
            formula_dict[c_l][c_l + '34'] = f'={c_l}32/{c_l}23'
            formula_dict[c_l][c_l + '35'] = f'={c_l}22*{c_l}34'
            formula_dict[c_l][c_l + '36'] = f'={c_l}22*{c_l}33'

    n = 'tss_column_len' # placeholder
    tss_data = '' # from tss data
    rpk_ratio_dict = defaultdict()

    print(formula_dict['top_row'])
    #for k,v in smb.items():
        
        
    raw_counts_top_row = defaultdict()
    r_top_start = 70 # raw_counts tab starts with F
    for k,v in formula_dict['top_row'].items():
        raw_counts_top_row[form_cell(r_top_start) + 1] = v
        r_top_start += 1
        
    print(f'raw_counts_top_row :=========== {raw_counts_top_row}')
        
        
        
        
        
    formula_rpk_dict = defaultdict()
    formula_rpk_dict['top_row'] = defaultdict()
    # forming a dict for rpk_count_ratio top row cells
    rpk_tab_top_row_dict = {'group_leader': '_rpk', 
                            'input': '_rpk',
                            'input_leader': '_rpk',
                            'sample_': ['_rpk', '_rpk_ratio', '_rpk_10rpkfloor', '_rpk_readnorm', '_rpk_readnorm_ratio',
                                        '_rpk_readnorm_10rpkfloor'],
                            'sample_sp': ['_rpk_flynorm', '_rpk_flynorm_ratio', '_rpk_flynorm_10rpkfloor']}



    print(formula_dict)

    #with pd.ExcelWriter('tool/bedtools_summary/formula_test/chipseq_230524_formula_test.xlsx', engine='xlsxwriter') as writer:
    #with pd.ExcelWriter('tool/bedtools_summary/formula_test/chIPseq_230809_formula_manual_test.xlsx', engine='xlsxwriter') as writer:

    '''
    with pd.ExcelWriter('tool/bedtools_summary/formula_test/chipseq_230601_formula_test_new.xlsx', engine='xlsxwriter') as writer:
        #test_df.to_excel(writer, sheet_name='sample_reads', index=False)
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

        tss.to_excel(writer, sheet_name='raw_counts', index=False)
        workbook = writer.book
        worksheet = writer.sheets['raw_counts']
        digit_format = workbook.add_format({'num_format': '#,##0'})
        decimal_format = workbook.add_format({'num_format': '##0.00'})
        #worksheet = workbook.add_worksheet('raw_counts')
    '''


make_formula()
