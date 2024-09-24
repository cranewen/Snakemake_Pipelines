# create a spreadsheet with all the colors with its hex code
from utility import YamlDict, YamlReader
import yaml
import sys, os
sys.path.insert(0, '/home/yw900/lab_pipelines')
from collections import defaultdict
import pandas as pd


def form_cell(n):
    A = 65
    Z = 90
    R = 26
    result = None
    # l1, l2, l3 are the letters from right to left => l3l2l1 or l2l1 or l1
    l1 = ''
    l2 = ''
    l3 = ''
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

def write_spreadsheet():
    color_yaml = 'metaconfig/bedtool_summary_color_config.yaml'
    color_dict = defaultdict()
    color = YamlDict(color_yaml)
    a = color.read_yaml() # dict of the entire yaml
    color_rgb = [*a]
    color_dict = defaultdict()
    for c in color_rgb:
        sub = a[c]['sub']
        sub.insert(0, a[c]['base'])
        color_dict[c] = sub
     
    r_c = 1 # row_count
    start = 65
    with pd.ExcelWriter('xlsx_test/color_table.xlsx') as writer:
        workbook = writer.book
        worksheet = workbook.add_worksheet('color_table')
        for k,v in color_dict.items():
            cell = form_cell(start) # cell letter
            worksheet.write(cell + str(r_c), k)
            r_c += 1
            for v_i in v:
                color_format = workbook.add_format({'bg_color': v_i})
                worksheet.write(cell + str(r_c), v_i, color_format)
                r_c += 1
            r_c = 1
            start += 1    
                

write_spreadsheet()
