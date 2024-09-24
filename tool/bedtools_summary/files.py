# only read files and store the data in a generator
from typing import Dict
from collections import defaultdict
import pandas as pd
from pprint import pprint

class Files:
    def __init__(self, filepath: str):
        self.filepath = filepath
        
    def __repr__(self):
        print(f'===== Files object members =====')
        pprint(self.__dict__)
        return f'===== Files object members ====='     
    
    def read_txt_file(self) -> Dict:
        stats_dict = defaultdict()
        with open(self.filepath) as f:
            for line in f:
                yield line
    
    def read_csv_to_df():
        df = pd.read_csv(self.filepath)
        return df


def main():
    f = Files('/mnt/storage/dept/pedonc/CPCT/projects/Partner_lab_projects/BAM_sorted/ChIPseq_221118_JR_hg38/samtools_flagstat/A549_ChIP_input_10ug_dex_10uM_2h_REP1_221118_TGACCA_S4_hg38_sorted_readgps_samtools_flagstat.txt')
    f.read_txt_file()
    for i in f.read_txt_file():
        print(i)
    print(f)


if __name__ == '__main__':
    main()
                    
