# a class deals with sample names properties
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
sys.path.insert(0, '/home/yw900/lab_pipelines')
from collections import defaultdict
import yaml
from db.dbconnection import DBconnector
from db.cpctdb import DatasetSeqRunSample, DatasetSeqRun, Sample, SampleGroup, SampleAttribute, Protocol
from metaparser.sampleparser import SampleParser
from pprint import pprint
from typing import Dict
from dataclasses import dataclass, is_dataclass, make_dataclass, fields, asdict, astuple

@dataclass
class Sample:
    dataset: str  # dataset_name from Dataset
    db_connection: type(DBconnector) = DBconnector
    samples: list = None
    sample_input_dict: Dict = None
    
    
'''
class Sample:
    def __init__(self, dataset):
        self.dataset = dataset
        self.db_connection = DBconnector()
        self.samples = None
        self.sample_input_dict = None
        with self.db_connection.create_session() as session:
            self.protocol = session.query(Protocol.protocolName).filter(DatasetSeqRun.protocolID == Protocol.protocolID).filter(DatasetSeqRun.datasetSeqRunDirectory == self.dataset).first()[0]
        

    def __repr__(self):
        print(f'===== Sample object members =====')
        pprint(self.__dict__)
        return f'===== Sample object members ====='

    # get sample names, sample input dict
    def set_samples(self):
        sample_parser = SampleParser()
        sample_parser.get_sample_list(self.dataset)
        self.samples = sample_parser.sample_list
        if self.protocol == 'ChIPseq':
            sample_parser.get_sample_input_dict(self.dataset)
            self.sample_input_dict = sample_parser.sample_input_dict

'''
    
        
        
def main():
    #s = Sample('ChIPseq_221118_JR')
    s = Sample('RNAseq_220708')
    #s.set_samples()
    print(s)


if __name__ == '__main__':
    main()
