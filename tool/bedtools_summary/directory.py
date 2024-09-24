import sys, os
#sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
#sys.path.append('../../../lab_pipelines')
sys.path.insert(0, '/home/yw900/lab_pipelines')
from metaparser.dirbuilder import DirBuilder
from collections import defaultdict
from pprint import pprint, pformat

class Directory:
    def __init__(self, dataset: str, ref_name: str):
        self.dataset = dataset
        self.ref_name = ref_name
        self.dirs = None
        self.dirs_spikein = None
        self.flagstat_path = defaultdict()

    def __repr__(self):
        print(f'===== Directory object members =====')
        pprint(self.__dict__)
        return f'===== Directory object members ====='

    # set all the required dirs from DirBuilder, no return
    def set_dirs(self) -> None:
        self.dirs = DirBuilder(self.dataset)
        self.dirs.build_chipseq_dirs(self.ref_name) 
        self.tss_all_samples_path = self.dirs.bedtools_coverage_dataset_rmdup_dir + '/tss_all_samples.csv'
        self.flagstat_path['sorted'] = self.dirs.sorted_bam_dataset_dir + 'samtools_flagstat/'
        self.flagstat_path['rmdup'] = self.dirs.rmdup_bam_dataset_dir + 'samtools_flagstat/'

    def set_dirs_spikein(self, ref_genome = 'dm6') -> None:
        self.dirs_spikein = DirBuilder(self.dataset)
        self.dirs_spikein.build_chipseq_spikein_dirs('dm6')
        self.tss_all_samples_path_spikein = self.dirs_spikein.bedtools_coverage_dataset_rmdup_dir + '/tss_all_samples.csv'
        self.flagstat_path['sorted_dm6'] = self.dirs_spikein.sorted_bam_dataset_dir + 'samtools_flagstat/'
        self.flagstat_path['rmdup_dm6'] = self.dirs_spikein.rmdup_bam_dataset_dir + 'samtools_flagstat/'

    '''
    def combine_flagstat_path(self):
        self.set_dirs()
        flagstat_path = self.flagstat_path
        self.set_dirs_spikein(self.ref_name)
        flagstat_path_spikein = self.flagstat_path
        self.flagstat_path = {**flagstat_path, **flagstat_path_spikein}
    '''

def main():
    d = Directory('ChIPseq_221118_JR', 'hg38')
    d.set_dirs()
    d_s = Directory('ChIPseq_221118_JR', 'dm6')
    print(d)
    d_s.set_dirs_spikein()
    print({**d.flagstat_path, **d_s.flagstat_path})
    #print(d.tss_all_samples_path)
    print(d.tss_all_samples_path_spikein)
    print(d)
    

if __name__ == "__main__":
    main()
