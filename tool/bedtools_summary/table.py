# preps all the excel sheets related stuffs
from collections import defaultdict
from typing import Dict, List
from pprint import pprint

class Table:
    def __init__(self):
        self.rownames = None
        self.matrix = None

    def __repr__(self):
        print(f'===== Table object members =====')
        pprint(self.__dict__)
        return f'===== Table object members ====='

    # a row name dictionary of {row#:name}
    def gen_rowname_dict(self, ref_name: str) -> Dict:
        rowname_dict = {1: 'rmdup paired reads (1)',
                        2: 'rmdup % total',
                        3: 'rmdup paired read ratio (2)', 
                        4: 'rmdup ' + ref_name + ' aligned reads (1)', 
                        5: 'rmdup ' + ref_name + ' aligned read ratio (4)',
                        6: 'rmdup dm6 aligned reads (1)',
                        7: 'rmdup all aligned reads (3)',
                        8: 'rmdup % aligned',
                        9: 'rmdup spikein reads per million (5)',
                        10: 'rmdup spikein rpm ratio (5b)',
                        11: 'rmdup input normalization factor (6a)',
                        12: 'rmdup input-normalized ' + ref_name + ' reads (7a)',
                        13: 'rmdup input-normalized read ratio (8a)',
                        14: 'rmdup read-adj input-norm read ratio (9a)', 
                        15: 'rmdup read-adj spikein read ratio (9b)',
                        16: 'rmdup final norm factor (input-norm) (10a)',
                        17: 'rmdup expected ' + ref_name + ' aligned reads (input-norm) (11a)',
                        18: 'rmdup expected ' + ref_name + ' aligned reads (spikein rpm) (11b)',
                        19: 'dupseq paired reads (1)',
                        20: 'dupseq paired read ratio (2)',
                        21: 'dupseq ' + ref_name + ' aligned reads (1)',
                        22: 'dupseq ' + ref_name + ' aligned read ratio (4)',
                        23: 'dupseq dm6 aligned reads (1)',
                        24: 'dupseq all aligned reads (3)',
                        25: 'dupseq % aligned',
                        26: 'dupseq spikein reads per million (5)',
                        27: 'dupseq spikein rpm ratio (5b)',
                        28: 'dupseq input normalization factor (6a)',
                        29: 'dupseq input-normalized ' + ref_name + ' reads (7a)',
                        30: 'dupseq input-normalized read ratio (8a)',
                        31: 'dupseq read-adj input-norm read ratio (9a)',
                        32: 'dupseq read-adj spikein read ratio (9b)',
                        33: 'dupseq final norm factor (input-norm) (10a)',
                        34: 'dupseq expected ' + ref_name + ' aligned reads (input-norm) (11a)',
                        35: 'dupseq expected ' + ref_name + ' aligned reads (spikein rpm) (11b)'}

        rowname_dict2 = {'r1': 'rmdup paired reads (1)',
                        'r2': 'rmdup % total',
                        'r3': 'rmdup paired read ratio (2)', 
                        'r4': 'rmdup ' + ref_name + ' aligned reads (1)', 
                        'r5': 'rmdup ' + ref_name + ' aligned read ratio (4)',
                        'r6': 'rmdup dm6 aligned reads (1)',
                        'r7': 'rmdup all aligned reads (3)',
                        'r8': 'rmdup % aligned',
                        'r9': 'rmdup spikein reads per million (5)',
                        'r10': 'rmdup spikein rpm ratio (5b)',
                        'r11': 'rmdup input normalization factor (6a)',
                        'r12': 'rmdup input-normalized ' + ref_name + ' reads (7a)',
                        'r13': 'rmdup input-normalized read ratio (8a)',
                        'r14': 'rmdup read-adj input-norm read ratio (9a)', 
                        'r15': 'rmdup read-adj spikein read ratio (9b)',
                        'r16': 'rmdup final norm factor (input-norm) (10a)',
                        'r17': 'rmdup expected ' + ref_name + ' aligned reads (input-norm) (11a)',
                        'r18': 'rmdup expected ' + ref_name + ' aligned reads (spikein rpm) (11b)',
                        'r19': 'dupseq paired reads (1)',
                        'r20': 'dupseq paired read ratio (2)',
                        'r21': 'dupseq ' + ref_name + ' aligned reads (1)',
                        'r22': 'dupseq ' + ref_name + ' aligned read ratio (4)',
                        'r23': 'dupseq dm6 aligned reads (1)',
                        'r24': 'dupseq all aligned reads (3)',
                        'r25': 'dupseq % aligned',
                        'r26': 'dupseq spikein reads per million (5)',
                        'r27': 'dupseq spikein rpm ratio (5b)',
                        'r28': 'dupseq input normalization factor (6a)',
                        'r29': 'dupseq input-normalized ' + ref_name + ' reads (7a)',
                        'r30': 'dupseq input-normalized read ratio (8a)',
                        'r31': 'dupseq read-adj input-norm read ratio (9a)',
                        'r32': 'dupseq read-adj spikein read ratio (9b)',
                        'r33': 'dupseq final norm factor (input-norm) (10a)',
                        'r34': 'dupseq expected ' + ref_name + ' aligned reads (input-norm) (11a)',
                        'r35': 'dupseq expected ' + ref_name + ' aligned reads (spikein rpm) (11b)'}


        return rowname_dict

    def gen_matrix(self) -> Dict:
        
         

    def excel_cell_num(self, num: int) -> List:
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
    
    def color_picking(self, num: int) -> List:
        if num == 1:
            return 7 // 2
        return [x * (7 // num) for x in range(1, num + 1)]

def main():
    t = Table()    
    print(t.gen_rowname_dict('hg38'))
    print(t.excel_cell_num(100))
    print(t.color_picking(1))

if __name__ == "__main__":
    main()
