# Attributes: dataset_name, ref_genome, spikein, protocal: chipseq, rnaseq,etc. 
from pprint import pprint
from dataclasses import dataclass, is_dataclass, make_dataclass, fields, asdict, astuple

@dataclass
class Dataset:
    dataset_name: str
    ref_genome: str
    reads_norm: bool = True
    flyreads: bool = False
    
    def __repr__(self):
        print(f'===== Dataset object members =====')
        pprint(self.__dict__)
        return f'===== Dataset object members ====='

'''
class Dataset:
    def __init__(self):
        self.dataset_name = None
        self.ref_genome = None
        self.reads_norm: bool = True
        self.flyreads: bool = False
    

    def __repr__(self):
        print(f'===== Dataset object members =====')
        pprint(self.__dict__)
        return f'===== Dataset object members ====='
        
    def set_dataset_name(self, dataset_name: str) -> None:
        self.dataset_name = dataset_name

    def set_ref_genome(self, ref_genome: str) -> None:
        self.ref_genome = ref_genome

    def set_spikein_attr(self, reads_norm: bool, flyreads: bool) -> None:
        self.reads_norm = reads_norm
        self.flyreads = flyreads

'''

def main():
    d = Dataset('ChIPseq_221118_JR', 'hg38', True)
    print(asdict(d))
    print(d.ref_genome)
    '''
    d.set_dataset_name('ChIPseq_221118_JR')
    d.set_ref_genome('hg38')
    d.set_spikein_attr(reads_norm = True, flyreads = True)
    '''
    print(d)


if __name__ == '__main__':
    main()
