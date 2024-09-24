import argparse
from data import *
from process import *
from utility import YamlReader
from reconfig import Reconf

parser = argparse.ArgumentParser(description='Bedtool Summary Pipeline')
parser.add_argument('-c', '--config', type=str)
parser.add_argument('--custom', action=argparse.BooleanOptionalAction)
parser.add_argument('-d', '--dataset', type=str)
parser.add_argument('-ref', '--ref_genome', type=str)
parser.add_argument('-ref_sp', '--ref_genome_spikein', type=str)
parser.add_argument('-dset', '--datasets', type=str) # e.g "ChIPseq_230123|ChIPseq_230201" etc. using | as the delimiter 
parser.add_argument('-tss', '--tss_file', type=str)
parser.add_argument('-o', '--output', type=str, required=True)
args = parser.parse_args()

def reorder_dict(sample_dict: Dict = None, new_order: list = None):
    new_sample_dict = defaultdict() 
    try:
        new_sample_dict = {k: sample_dict[k] for k in new_order}
    except KeyError:
        print('Please check if the sample list has problem!')
    return new_sample_dict
    

def run_tool_custom(dataset, ref, ref_sp, tss, output):
    #dset = Dataset('ChIPseq_230209_hs_combine_test', 'hg38', spikein_ref_genome = 'dm6')
    dset = Dataset(dataset, ref, spikein_ref_genome = ref_sp)
    #dset.protocol = 'ChIPseq_custom'
    print(dset)
    p = DataParserCustom(dset)
    #p.parse_dirs('xlsx_test/combined_tss.csv') # args.tss_file
    p.parse_dirs(tss) # args.tss_file
    p.parse_sample_meta()
    df = DF(p.sample_meta, dset)
    print(list(p.sample_meta.keys()))
    #tss = p.parse_tss(list(p.sample_meta.keys()))
    sample_order = []
    for k,v in p.parse_strategy().items():
        sample_order.append(k)
        for v1 in v['group_samples']:
            sample_order.append(v1)
        
    #print(sample_order)
    tss = p.parse_tss(sample_order)
    print(tss)
    #print(tss)
    #print(p.sample_meta)
    raw_stats = df.calculate_sample_stats(p.sample_meta)

    raw_stats_copy = raw_stats.copy()
    raw_stats = reorder_dict(raw_stats_copy, sample_order)
    

    for k,v in raw_stats.items():
        v.update_rs()
        raw_stats[k] = v
    df.calculate_rpk_count_ratio(raw_stats, tss)

    sample_reads_df = defaultdict()
    t = Table(ref_name=dset.ref_genome)
    sample_reads_df['index'] = [v for k,v in t.sample_reads_index.items()]
    for k,v in raw_stats.items():
        print(k)
        n_k = p.sample_meta[k].short_name
        sample_reads_df[n_k] = v.r_np
    #new_sample_reads_df = sample_reorder(sample_reads_df, sample_order)
    #sample_reads_df_out = pd.DataFrame(new_sample_reads_df)
    sample_reads_df_out = pd.DataFrame(sample_reads_df)
    print(raw_stats.keys())

    tss_out = pd.DataFrame(tss)

    #sample_name_for_color = p.sample_meta.keys()
    sample_name_for_color = raw_stats.keys()
    c1_temp = zip(color_cell(len(p.sample_meta) + 1)[1:], list(sample_name_for_color))
    c2_temp = zip(color_cell(len(p.sample_meta) + 5)[5:], list(sample_name_for_color))

    s_m = defaultdict() # sample meta from sample_meta yaml file
    y = YamlReader('metaconfig/sample_meta.yaml')
    for i in y.read_yaml():
        s_m = i[dset.dataset_name]    
    
    c1 = defaultdict()
    c2 = defaultdict()
    c3 = defaultdict()

    # c = {cell:[shortname, color]}
    for c in c1_temp:
        c1[c[0]] = [s_m[c[1]]['short_name'], s_m[c[1]]['color']]
    for c in c2_temp:
        c2[c[0]] = [s_m[c[1]]['short_name'], s_m[c[1]]['color']]
    
    c3_sample_list = []
    if dset.spikein_ref_genome:
        for s in sample_name_for_color:
            if not s_m[s]['is_input'] and not s_m[s]['is_group_leader']:
                c3_sample_list = c3_sample_list + [s] * 9
            elif not s_m[s]['is_input'] and s_m[s]['is_group_leader']:
                c3_sample_list.append(s)
            elif s_m[s]['is_input'] and not s_m[s]['is_group_leader']:
                c3_sample_list.append(s)
            elif s_m[s]['is_input'] and s_m[s]['is_group_leader']:
                c3_sample_list.append(s)
                    
    else:
        for s in sample_name_for_color:
            if not s_m[s]['is_input'] and not s_m[s]['is_group_leader']:
                c3_sample_list = c3_sample_list + [s] * 6
            elif not s_m[s]['is_input'] and s_m[s]['is_group_leader']:
                c3_sample_list.append(s)
            elif s_m[s]['is_input'] and not s_m[s]['is_group_leader']:
                c3_sample_list.append(s)
            elif s_m[s]['is_input'] and s_m[s]['is_group_leader']:
                c3_sample_list.append(s)
        
    c3_temp = zip(color_cell(len(list(df.rpk_count_ratios_df.keys())[6:]) + 6)[6:], list(df.rpk_count_ratios_df.keys())[6:], c3_sample_list)

    for c in c3_temp:
        c3[c[0]] = [c[1], s_m[c[2]]['color']]
    

    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
    #with pd.ExcelWriter('xlsx_test/ChIPseq_230209_hs_combine_test_args_new_order_with_color.xlsx', engine='xlsxwriter') as writer:
        sample_reads_df_out.to_excel(writer, sheet_name='sample_reads', index=False)
        workbook = writer.book
        worksheet = writer.sheets['sample_reads']
        a1_format = workbook.add_format()
        worksheet.write('A1', None, a1_format)
        digit_format = workbook.add_format({'num_format': '#,##0'})
        decimal_format = workbook.add_format({'num_format': '##0.00'})
        for k,v in c1.items():
            color_format = workbook.add_format({'bg_color': v[1]})
            worksheet.write(k, v[0], color_format)
        last_cell = color_cell(len(sample_reads_df_out.columns))[-1][:-1]
        cell_range = 'B2:' + last_cell + str(len(sample_reads_df_out) + 1)
        worksheet.conditional_format(cell_range, {'type': 'cell', 'criteria': '>', 'value': 100, 'format': digit_format})
        worksheet.conditional_format(cell_range, {'type': 'cell', 'criteria': '<=', 'value': 100, 'format': decimal_format})

        tss_out.to_excel(writer, sheet_name='raw_counts', index=False)
        workbook = writer.book
        worksheet = writer.sheets['raw_counts']
        digit_format = workbook.add_format({'num_format': '#,##0'})
        decimal_format = workbook.add_format({'num_format': '##0.00'})
        for k,v in c2.items():
            color_format = workbook.add_format({'bg_color': v[1]})
            worksheet.write(k, v[0], color_format)

        pd.DataFrame(df.rpk_count_ratios_df).to_excel(writer, sheet_name='rpk_count_ratios', index=False)
        workbook = writer.book
        worksheet = writer.sheets['rpk_count_ratios']
        digit_format = workbook.add_format({'num_format': '#,##0'})
        decimal_format = workbook.add_format({'num_format': '##0.00'})
        for k,v in c3.items():
            color_format = workbook.add_format({'bg_color': v[1]})
            worksheet.write(k, v[0], color_format)


        pd.DataFrame(df.rpk_count_ratios_per_gene_df).to_excel(writer, sheet_name='rpk_count_ratios_per_gene', index=False)
        workbook = writer.book
        worksheet = writer.sheets['rpk_count_ratios_per_gene']
        digit_format = workbook.add_format({'num_format': '#,##0'})
        decimal_format = workbook.add_format({'num_format': '##0.00'})
        for k,v in c3.items():
            color_format = workbook.add_format({'bg_color': v[1]})
            worksheet.write(k, v[0], color_format)


def prep_custom(datasets, ref, ref_sp):
    ds = datasets.split('|')
    for d in ds:
        dset = Dataset(d, ref, spikein_ref_genome = ref_sp)
        print(dset)
        datab = Database(dset)
        datab._get_sample_info()
        datab._protocol_check()
        print(f'PROTOCOL is : {datab.dataset.protocol}')
        p = DataParser(dset)
        p.parse_dirs()
        sample_meta_with_group = p.add_group_info(datab.sample_meta)
        p.add_dirs(datab.sample_meta)
        for k,v in datab.sample_meta.items():
            datab.sample_meta[k] = p._parse_flagstat(v, dset.flyreads)
        df = DF(datab.sample_meta, dset)
        dw = DataWriter(datab.sample_meta, dset)
        dw.write_strategy('metaconfig/bedtool_summary_config.yaml')
        dw.write_sample_meta()

def rewrite_sample_meta(dataset, datasets):
    datasets = datasets.split('|')
    reconf = Reconf(dataset, datasets)
    reconf.get_orig_sample_meta()
    reconf.rewrite_sample_meta(reconf.get_new_strategy())
    

def run_tool(dataset, ref, ref_sp, output):
    dset = Dataset(dataset, ref, spikein_ref_genome = ref_sp)
    #dset = Dataset('ChIPseq_221118_JR', 'hg38')
    #dset = Dataset('ChIPseq_221114', 'mm10', spikein_ref_genome = 'dm6')
    #dset = Dataset('ChIPseq_230127_M7', 'hg38')
    #dset = Dataset('ChIPseq_230209_hs', 'hg38', spikein_ref_genome = 'dm6')
    #dset = Dataset('ChIPseq_230209_hs_combine_test', 'hg38', spikein_ref_genome = 'dm6')
    #dset = Dataset('ChIPseq_230127_M7', 'hg38', spikein_ref_genome = 'dm6')
    #dset = Dataset('ChIPseq_230113', 'mm10')
    datab = Database(dset)
    datab._protocol_check()
    if datab.dataset.protocol == 'RNAseq':
        datab._get_sample_info()
        p = DataParser(dset)
        p.parse_dirs()
        sample_meta_with_group = p.add_group_info(datab.sample_meta)
        p.add_dirs(datab.sample_meta)
        for k,v in datab.sample_meta.items():
            datab.sample_meta[k] = p._parse_flagstat(v, dset.flyreads)
        print(datab.sample_meta)
        df = DF(datab.sample_meta, dset)
        sample_name_for_color = datab.sample_meta.keys()
        print(df)
        dw = DataWriter(datab.sample_meta, dset)
        dw.write_strategy('metaconfig/bedtool_summary_config.yaml')
        dw.write_sample_meta()
        raw_stats = df.calculate_sample_stats(p.parse_sample_meta())
        for k,v in raw_stats.items():
            v.update_rs()
            raw_stats[k] = v

        sample_reads_df = defaultdict()
        t = Table(ref_name=dset.ref_genome)
        sample_reads_df['index'] = [v for k,v in t.sample_reads_index.items()]
        for k,v in raw_stats.items():
            n_k = datab.sample_meta[k].short_name
            sample_reads_df[n_k] = v.r_np
        sample_reads_df_out = pd.DataFrame(sample_reads_df)
        print(sample_reads_df_out)

        c1_temp = zip(color_cell(len(datab.sample_meta) + 1)[1:], list(sample_name_for_color))
        y = YamlReader('metaconfig/sample_meta.yaml')
        for i in y.read_yaml():
            s_m = i[dset.dataset_name]    
        
        c1 = defaultdict()
        for c in c1_temp:
            c1[c[0]] = [s_m[c[1]]['short_name'], s_m[c[1]]['color']]
        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
            sample_reads_df_out.to_excel(writer, sheet_name='sample_reads', index=False)
            workbook = writer.book
            worksheet = writer.sheets['sample_reads']
            a1_format = workbook.add_format()
            worksheet.write('A1', None, a1_format)
            digit_format = workbook.add_format({'num_format': '#,##0'})
            decimal_format = workbook.add_format({'num_format': '##0.00'})
            for k,v in c1.items():
                color_format = workbook.add_format({'bg_color': v[1]})
                worksheet.write(k, v[0], color_format)
            last_cell = color_cell(len(sample_reads_df_out.columns))[-1][:-1]
            cell_range = 'B2:' + last_cell + str(len(sample_reads_df_out) + 1)
            print(cell_range)
            worksheet.conditional_format(cell_range, {'type': 'cell', 'criteria': '>', 'value': 100, 'format': digit_format})
            worksheet.conditional_format(cell_range, {'type': 'cell', 'criteria': '<=', 'value': 100, 'format': decimal_format})
            border_line_range = 'A20:' + last_cell + '20'
            border_format = workbook.add_format()
            border_format.set_top(2)
            worksheet.conditional_format(border_line_range, {'type': 'cell', 'criteria': '>', 'value': 0, 'format': border_format})
    else:


        datab._get_sample_info()
        p = DataParser(dset)
        #p.parse_dirs('TSS_ALL_CUSTOMIZED')
        p.parse_dirs()
        print(p.d)
        #p.parse_strategy()
        sample_meta_with_group = p.add_group_info(datab.sample_meta)
        p.add_dirs(datab.sample_meta)
        for k,v in datab.sample_meta.items():
            datab.sample_meta[k] = p._parse_flagstat(v, dset.flyreads)
        tss = p.parse_tss(list(datab.sample_meta.keys()))
        print(datab.sample_meta)
        #print(tss)
        #print(sample_meta_with_group)
        #print(p.strategy)
        df = DF(datab.sample_meta, dset)
        print(dset.flyreads)

        sample_name_for_color = datab.sample_meta.keys()
       
        dw = DataWriter(datab.sample_meta, dset)
        dw.write_strategy('metaconfig/bedtool_summary_config.yaml')
        #dw.write_attribute()
        dw.write_sample_meta()
        
        #raw_stats = df.calculate_sample_stats(datab.sample_meta)
        raw_stats = df.calculate_sample_stats(p.parse_sample_meta())
        for k,v in raw_stats.items():
            v.update_rs()
            raw_stats[k] = v
        #print(raw_stats)
        df.calculate_rpk_count_ratio(raw_stats, tss)
        #rpk = df.calculate_rpk_count_ratio(raw_stats, tss)
        #print(rpk)
        #print(df.rpk_count_ratios_per_gene_df)
        #print(df.rpk_count_ratios_df)

        sample_reads_df = defaultdict()
        t = Table(ref_name=dset.ref_genome)
        sample_reads_df['index'] = [v for k,v in t.sample_reads_index.items()]
        for k,v in raw_stats.items():
            n_k = datab.sample_meta[k].short_name
            sample_reads_df[n_k] = v.r_np
        sample_reads_df_out = pd.DataFrame(sample_reads_df)
        print(sample_reads_df_out)

        tss_out = pd.DataFrame(tss)
        print(tss_out)

        c1_temp = zip(color_cell(len(datab.sample_meta) + 1)[1:], list(sample_name_for_color))
        c2_temp = zip(color_cell(len(datab.sample_meta) + 5)[5:], list(sample_name_for_color))

        s_m = defaultdict() # sample meta from sample_meta yaml file
        y = YamlReader('metaconfig/sample_meta.yaml')
        for i in y.read_yaml():
            s_m = i[dset.dataset_name]    
        
        c1 = defaultdict()
        c2 = defaultdict()
        c3 = defaultdict()

        # c = {cell:[shortname, color]}
        for c in c1_temp:
            c1[c[0]] = [s_m[c[1]]['short_name'], s_m[c[1]]['color']]
        for c in c2_temp:
            c2[c[0]] = [s_m[c[1]]['short_name'], s_m[c[1]]['color']]
        
        c3_sample_list = []
        if dset.spikein_ref_genome:
            for s in sample_name_for_color:
                if not s_m[s]['is_input'] and not s_m[s]['is_group_leader']:
                    c3_sample_list = c3_sample_list + [s] * 9
                elif not s_m[s]['is_input'] and s_m[s]['is_group_leader']:
                    c3_sample_list.append(s)
                elif s_m[s]['is_input'] and not s_m[s]['is_group_leader']:
                    c3_sample_list.append(s)
                elif s_m[s]['is_input'] and s_m[s]['is_group_leader']:
                    c3_sample_list.append(s)
                        
        else:
            for s in sample_name_for_color:
                if not s_m[s]['is_input'] and not s_m[s]['is_group_leader']:
                    c3_sample_list = c3_sample_list + [s] * 6
                elif not s_m[s]['is_input'] and s_m[s]['is_group_leader']:
                    c3_sample_list.append(s)
                elif s_m[s]['is_input'] and not s_m[s]['is_group_leader']:
                    c3_sample_list.append(s)
                elif s_m[s]['is_input'] and s_m[s]['is_group_leader']:
                    c3_sample_list.append(s)
        print(len(c3_sample_list))
            
        c3_temp = zip(color_cell(len(list(df.rpk_count_ratios_df.keys())[6:]) + 6)[6:], list(df.rpk_count_ratios_df.keys())[6:], c3_sample_list)
        print(len(df.rpk_count_ratios_df.keys()))

        for c in c3_temp:
            c3[c[0]] = [c[1], s_m[c[2]]['color']]
        

        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        #with pd.ExcelWriter('xlsx_test/chipseq_221114_test_args_with_color.xlsx', engine='xlsxwriter') as writer:
        #with pd.ExcelWriter('xlsx_test/ChIPseq_230209_hs_combine_test_with_color.xlsx', engine='xlsxwriter') as writer:
            sample_reads_df_out.to_excel(writer, sheet_name='sample_reads', index=False)
            workbook = writer.book
            worksheet = writer.sheets['sample_reads']
            a1_format = workbook.add_format()
            worksheet.write('A1', None, a1_format)
            digit_format = workbook.add_format({'num_format': '#,##0'})
            decimal_format = workbook.add_format({'num_format': '##0.00'})
            for k,v in c1.items():
                color_format = workbook.add_format({'bg_color': v[1]})
                worksheet.write(k, v[0], color_format)
            last_cell = color_cell(len(sample_reads_df_out.columns))[-1][:-1]
            cell_range = 'B2:' + last_cell + str(len(sample_reads_df_out) + 1)
            worksheet.conditional_format(cell_range, {'type': 'cell', 'criteria': '>', 'value': 100, 'format': digit_format})
            worksheet.conditional_format(cell_range, {'type': 'cell', 'criteria': '<=', 'value': 100, 'format': decimal_format})

            tss_out.to_excel(writer, sheet_name='raw_counts', index=False)
            workbook = writer.book
            worksheet = writer.sheets['raw_counts']
            digit_format = workbook.add_format({'num_format': '#,##0'})
            decimal_format = workbook.add_format({'num_format': '##0.00'})
            for k,v in c2.items():
                color_format = workbook.add_format({'bg_color': v[1]})
                worksheet.write(k, v[0], color_format)

            pd.DataFrame(df.rpk_count_ratios_df).to_excel(writer, sheet_name='rpk_count_ratios', index=False)
            workbook = writer.book
            worksheet = writer.sheets['rpk_count_ratios']
            digit_format = workbook.add_format({'num_format': '#,##0'})
            decimal_format = workbook.add_format({'num_format': '##0.00'})
            for k,v in c3.items():
                color_format = workbook.add_format({'bg_color': v[1]})
                worksheet.write(k, v[0], color_format)


            pd.DataFrame(df.rpk_count_ratios_per_gene_df).to_excel(writer, sheet_name='rpk_count_ratios_per_gene', index=False)
            workbook = writer.book
            worksheet = writer.sheets['rpk_count_ratios_per_gene']
            digit_format = workbook.add_format({'num_format': '#,##0'})
            decimal_format = workbook.add_format({'num_format': '##0.00'})
            for k,v in c3.items():
                color_format = workbook.add_format({'bg_color': v[1]})
                worksheet.write(k, v[0], color_format)


#run_tool()
#run_tool_custom()

def run():
    if not args.custom:
        if not args.config:
            run_tool(args.dataset, args.ref_genome, args.ref_genome_spikein, args.output)
    else:
        if not args.config:
            if args.datasets and not args.dataset:
                prep_custom(args.datasets, args.ref_genome, args.ref_genome_spikein)
            if args.datasets and args.dataset:
                rewrite_sample_meta(args.dataset, args.datasets)
            if not args.datasets and args.tss_file and args.dataset:
                run_tool_custom(args.dataset, args.ref_genome, args.ref_genome_spikein, args.tss_file, args.output)

run()
    

