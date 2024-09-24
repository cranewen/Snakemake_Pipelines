from collections import defaultdict
import pandas as pd
from process import *
from data import *
#from openpyxl import Workbook

def make_formula():
    #input_group = {'B': ['C', 'D']}
    #sample_groups = {'E': ['F', 'G'], 'H':['I', 'J']}
    #input_dict = {'E': 'B', 'H': 'B', 'F': 'C', 'I': 'C', 'G': 'D', 'J': 'D'}

    ig = 'BCD'
    sg = 'EFG|HIJ'
    ind = 'BEH|CFI|DGJ'

    input_group = {ig[0]: [x for x in ig[1:]]}
    sample_groups = {g[0]: [x for x in g[1:]] for g in sg.split('|')}
    input_dict =  {x: y[0] for y in ind.split('|') for x in y[1:]}
        
    print(input_group)
    print(sample_groups)
    print(input_dict)

    col_start = 66 # ascii of B
    # a newly reformed dictionary with keys: "B,C,D,E..." 
    # the values are a list of selected attributes from sample_meta_dict
    f_d = defaultdict() 
    

    #input0s = list(input_group.keys())  # get all the input leaders
    input0s = [*input_group]  # get all the input leaders


    # Deep nests warning!!! Stupid but making perfect sense with the 3 dictionaries initials structure above
    formula_dict = defaultdict() # the entire sample_reads tab 
    for i0 in input0s:
        formula_dict[i0] = defaultdict()
        
        formula_dict[i0] = {i0 + '3': f'={i0}2/{i0}20*100',
                            i0 + '8': f'={i0}5 + {i0}7', 
                            i0 + '9': f'={i0}8/{i0}2*100',
                            i0 + '10': f'={i0}7/({i0}8/1000000)',
                            i0 + '19': f'={i0}5',
                            i0 + '25': f'={i0}22 + {i0}24',
                            i0 + '26': f'={i0}24/{i0}25*100',
                            i0 + '27': f'={i0}24/({i0}25/1000000)',
                            i0 + '36': f'={i0}22'}
        input1s = input_group[i0]
        for i1 in input1s:
            formula_dict[i1] = defaultdict()

            formula_dict[i1] = {i1 + '3': f'={i1}2/{i1}20*100',
                                i1 + '4': f'={i1}2/{i0}2', # input1, 2, 3...
                                i1 + '6': f'={i1}5/{i0}5', # input1, 2, 3...
                                i1 + '8': f'={i1}5 + {i1}7', 
                                i1 + '9': f'={i1}8/{i1}2*100',
                                i1 + '10': f'={i1}7/({i1}8/1000000)',
                                i1 + '11': f'={i1}10/{i0}10', # input1, 2, 3...
                                i1 + '16': f'=(1/{i1}11)/{i1}4', # input1, 2, 3...
                                i1 + '19': f'={i1}5*{i1}16',
                                i1 + '21': f'={i1}20/{i0}20', # input1, 2, 3...
                                i1 + '23': f'={i1}22/{i0}22', # input1, 2, 3...
                                i1 + '25': f'={i1}22 + {i1}24',
                                i1 + '26': f'={i1}24/{i1}25*100',
                                i1 + '27': f'={i1}24/({i1}25/1000000)',
                                i1 + '28': f'={i1}27/{i0}27', # input1, 2, 3...
                                i1 + '33': f'=(1/{i1}28)/{i1}21)', # input1, 2, 3...
                                i1 + '36': f'={i1}22'}
        control0s = [*sample_groups]
        for c0 in control0s:
            formula_dict[c0] = defaultdict()

            formula_dict[c0] = {c0 + '3': f'={c0}2/{c0}20*100',
                                c0 + '8': f'={c0}5 + {c0}7', 
                                c0 + '9': f'={c0}8/{c0}2*100',
                                c0 + '10': f'={c0}7/({c0}8/1000000)',
                                c0 + '12': f'={c0}10/{i0}10', # control0
                                c0 + '13': f'={c0}5/{c0}12', # control0
                                c0 + '18': f'={c0}5', # control0
                                c0 + '19': f'={c0}5',
                                c0 + '25': f'={c0}22 + {c0}24',
                                c0 + '26': f'={c0}24/{c0}25*100',
                                c0 + '27': f'={c0}24/({c0}25/1000000)',
                                c0 + '29': f'={c0}27/{i0}27', # control0
                                c0 + '30': f'={c0}22/{c0}29', # control0
                                c0 + '35': f'={c0}22', # control0
                                c0 + '36': f'={c0}22'}

            samples = sample_groups[c0]
            
            for s in samples:
                formula_dict[s] = defaultdict()

                formula_dict[s] = {s + '3': f'={s}2/{s}20*100',
                                   s + '4': f'={s}2/{c0}2',
                                   s + '6': f'={s}5/{c0}5',
                                   s + '8': f'={s}5 + {s}7', 
                                   s + '9': f'={s}8/{s}2*100',
                                   s + '10': f'={s}7/({s}8/1000000)',
                                   s + '11': f'={s}10/{c0}10',
                                   s + '12': f'={s}10/{input_dict[s]}10',
                                   s + '13': f'={s}5/{s}12',
                                   s + '14': f'={s}13/{c0}13',
                                   s + '15': f'={s}14/{s}4',
                                   s + '16': f'=(1/{s}11)/{s}4',
                                   s + '17': f'={s}15/{s}6',
                                   s + '18': f'={s}5*{s}17',
                                   s + '19': f'={s}5*{s}16',
                                   s + '21': f'={s}20/{c0}20',
                                   s + '23': f'={s}22/{c0}22',
                                   s + '25': f'={s}22 + {s}24',
                                   s + '26': f'={s}24/{s}25*100',
                                   s + '27': f'={s}24/({s}25/1000000)',
                                   s + '28': f'={s}27/{c0}27',
                                   s + '29': f'={s}27/{input_dict[s]}27',
                                   s + '30': f'={s}22/{s}29',
                                   s + '31': f'={s}30/{c0}30',
                                   s + '32': f'={s}31/{s}21',
                                   s + '33': f'=(1/{s}28)/{s}21',
                                   s + '34': f'={s}32/{s}23',
                                   s + '35': f'={s}22*{s}34',
                                   s + '36': f'={s}22*{s}33'}

   
    
    #print(formula_dict)

    #test_df = pd.DataFrame.from_dict(formula_dict)
    
    dset = Dataset('ChIPseq_230601', 'hg38', spikein_ref_genome = 'dm6')
    data_p = DataParser(dset)
    smd = data_p.parse_sample_meta()
    col_letter = defaultdict()
    for k in list(smd.keys()):
        col_letter[k] = chr(col_start)
        col_start += 1
    print(col_letter)
    s_list = list(smd.keys())
    selected_keys = ['cpct_full_name', 'group_id', 'group_position', 'group_leader', 'input_sample',
                     'prop_paired_rmdup', 'prop_paired_rmdup_sp', 'prop_paired_sorted', 'prop_paired_sorted_sp',
                     'seq_paired_rmdup', 'seq_paired_sorted', 'column_letter']
    f_d['rows'] = selected_keys

    for s in s_list:
        f_d[smd[s].short_name] = [smd[s].cpct_full_name, smd[s].group_id, smd[s].group_position, 
                               smd[s].group_leader, smd[s].input_sample, smd[s].prop_paired_rmdup, smd[s].prop_paired_rmdup_sp,
                               smd[s].prop_paired_sorted, smd[s].prop_paired_sorted_sp,
                               smd[s].seq_paired_rmdup, smd[s].seq_paired_sorted, chr(col_start)]
        #f_d[chr(col_start)] = [v for k,v in smd[s].__dict__.items()]
        col_start += 1

    #df = pd.DataFrame.from_dict(f_d)
    #df.to_csv('xlsx_test/formula_new.csv', index=False)

    '''
    with pd.ExcelWriter('formula_test/formula_test.xlsx', engine='xlsxwriter') as writer:
        #test_df.to_excel(writer, sheet_name='sample_reads', index=False)
        workbook = writer.book
        worksheet = workbook.add_worksheet('sample_reads')
        for k,v in formula_dict.items():
            for ki,vi in v.items():
                worksheet.write_formula(ki, vi)
    '''

make_formula()
