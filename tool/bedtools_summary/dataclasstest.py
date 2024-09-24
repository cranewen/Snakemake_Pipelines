from dataclasses import dataclass, is_dataclass, make_dataclass, fields, asdict, astuple

@dataclass
class Dataset:
    dataset_name: str
    ref_genome: str
    
    def justprint(self) -> str:
        return self.dataset_name

def main():
    d = Dataset('ChIPseq_123', 'hg38')
    print(asdict(d))

if __name__ == '__main__':
    main()
        
