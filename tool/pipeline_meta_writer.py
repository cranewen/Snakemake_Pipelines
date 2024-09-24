import sys, os, re
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from glob import glob
from datetime import datetime
from collections import defaultdict
from db.dbconnection import DBconnector
from db.cpctdb import Dataset, DatasetSeqRun, SequencingRun, Sample, SampleAttribute, SampleAttributeType, DatasetSeqRunSample, Project, Protocol
from metaparser.bclmetaparser import BclMetaParser

class PipelineMetaWriter:
    def __init__(self):
        self.db_connection = DBconnector()
        
        # set key as dataset_name(DatasetSeqRunDirectory in database)
        self.bcl_meta_dict = defaultdict()
        self.sample_meta_dict = defaultdict()
        pedonc_src_path = '/mnt/storage/dept/pedonc/src/'
        pipeline_root_path = 'lab_pipelines/'
        self.bcl_meta_path = pedonc_src_path + pipeline_root_path + 'metaconfig/bcl_meta.yaml'
        self.sample_yaml_path = pedonc_src_path + pipeline_root_path + 'metaconfig/sample.yaml'
        self.rnaseq_meta_yaml_path = pedonc_src_path + pipeline_root_path + 'metaconfig/rnaseq_meta.yaml'
        self.deseq2_meta_yaml_path = pedonc_src_path + pipeline_root_path + 'metaconfig/deseq2_meta.yaml'
        self.chipseq_meta_yaml_path = pedonc_src_path + pipeline_root_path + 'metaconfig/chipseq_meta.yaml'
        self.chipseq_input_sample_yaml_path = pedonc_src_path + pipeline_root_path + 'metaconfig/chipseq_input_sample.yaml'
        self.cuttag_meta_yaml_path = pedonc_src_path + pipeline_root_path + 'metaconfig/cuttag_meta.yaml'
        self.tornado_meta_yaml_path = pedonc_src_path + pipeline_root_path + 'metaconfig/tornado_meta.yaml'
        #self.crispr_meta_yaml_path = pedonc_src_path + pipeline_root_path + 'metaconfig/crispr_meta.yaml'
        self.crispr_meta_yaml_path = 'metaconfig/crispr_meta.yaml'
        self.sgrnaseq_meta_yaml_path = 'metaconfig/sgrnaseq_meta.yaml'

    def write_pipeline_yamls(self):
        with self.db_connection.create_session() as session:
            # select bcl's SequencingRunDate and DatasetSeqRunDirectory data in table "SequencingRun" and "DatasetSeqRun"
            bcl_date_dir_data = session.query(DatasetSeqRun.datasetSeqRunDirectory, SequencingRun.sequencingRunDate, DatasetSeqRun.bclYamlFlag).\
                join(DatasetSeqRun, SequencingRun.sequencingRunID == DatasetSeqRun.sequencingRunID).all()
            bcl_dirs = [a.datasetSeqRunDirectory for a in bcl_date_dir_data]
            # a datasetSeqRunDirectory list of labeled "N" in bclYamlFlag
            dataset_list = []

            # bcl_date_dir_data's bclYamlFlag will be shared for create sample
            for r in bcl_date_dir_data:
                if (r.bclYamlFlag == 'N' or r.bclYamlFlag == None):
                    raw_bcl_path = self.get_raw_bcl_path(str(r.sequencingRunDate.date().strftime("%y%m%d")))
                    self.bcl_meta_dict[r.datasetSeqRunDirectory] = \
                        ['  run_date: ' + str(r.sequencingRunDate.date()), '  raw_bcl_path: \"' + raw_bcl_path + '\"', '  sample_sheet_name: \"' + r.datasetSeqRunDirectory + '_sample_sheet.csv\"']

                    self.sample_meta_dict[r.datasetSeqRunDirectory] = session.query(Sample.sampleCPCTName, DatasetSeqRunSample.seqRunSampleBarcode).\
                        filter(Sample.sampleID == DatasetSeqRunSample.sampleID).\
                        filter(DatasetSeqRunSample.datasetSeqRunID == DatasetSeqRun.datasetSeqRunID).\
                        filter(DatasetSeqRun.datasetSeqRunDirectory == r.datasetSeqRunDirectory).all()

                    dataset_list.append(r.datasetSeqRunDirectory)

            self.write_meta_yaml(session, dataset_list)
            # session.query(DatasetSeqRun).filter(DatasetSeqRun.datasetSeqRunDirectory.in_(bcl_dirs)).update({"bclYamlFlag": "Y"}, synchronize_session = False)

            try:
                with open(self.bcl_meta_path, 'a') as f:
                    f.write('# ' + str(datetime.now().date()) + '\n')
                    for k,v in self.bcl_meta_dict.items():
                        f.write(k + ':\n')
                        for i in v:
                            f.write(i + '\n')
                        f.write('\n')

                with open(self.sample_yaml_path, 'a') as f:
                    f.write('# ' + str(datetime.now().date()) + '\n')
                    for k,v in self.sample_meta_dict.items():
                        self.write_sample_sheet(k)
                        f.write(k + ':\n')
                        for i in v:
                            # print(str(i))
                            f.write(' - ' + i[0] + '\n')
                        f.write('\n')
            except:
                raise
            else:
                print('Change N to Y')
                session.query(DatasetSeqRun).filter(DatasetSeqRun.datasetSeqRunDirectory.in_(bcl_dirs)).update({'bclYamlFlag': 'Y'}, synchronize_session = False)
            
    def write_sample_sheet(self, dataset_seq_run_directory):
        sample_sheet_path = '/mnt/storage/dept/pedonc/CPCT/projects/bcl2fastq_demultiplex/' + dataset_seq_run_directory + '_sample_sheet.csv'
        if not os.path.isfile(sample_sheet_path):
            sample_data = self.sample_meta_dict[dataset_seq_run_directory]
            with open(sample_sheet_path, 'w') as f:
                f.write('[Header],,,,,,,,,\n')
                f.write('IEMFileVersion,4,,,,,,,,\n')
                f.write('Date,' + datetime.now().strftime("%m/%d/%Y") + ',,,,,,,,\n')
                f.write('Workflow,GenerateFASTQ,,,,,,,,\n')
                f.write('Application,NextSeq FASTQ Only,,,,,,,,\n')
                f.write('Assay,TruSeq HT,,,,,,,,\n')
                f.write('Description,,,,,,,,,\n')
                f.write('Chemistry,Amplicon,,,,,,,,\n')
                f.write(',,,,,,,,,\n')
                f.write('[Reads],,,,,,,,,\n')
                f.write('37,,,,,,,,,\n')
                f.write(',,,,,,,,,\n')
                f.write(',,,,,,,,,\n')
                f.write('[Settings],,,,,,,,,\n')
                f.write('Adapter,,,,,,,,,\n')
                f.write('AdapterRead2,,,,,,,,,\n')
                f.write(',,,,,,,,,\n')
                f.write('[Data],,,,,,,,,\n')
                f.write('Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description\n')
                for i in sample_data:
                    f.write(i.sampleCPCTName + ',,,,,' + i.seqRunSampleBarcode + ',,,,\n')



    # writing ChIPseq (ChIPseq_input_sample yaml included), CUTTAGseq and Rnaseq meta yaml
    def write_meta_yaml(self, session, dataset_seq_run_directories):
        meta_dict = defaultdict()
        meta_data = ''
        meta_file_path = ''
        rnaseq_dataset_list = []
        chipseq_dataset_list = []
        cuttagseq_dataset_list = []
        # sgrnaseq is for crispr. crispr is not the downstream analysis, using sgrnaseq is just to navigate to the project
        sgrnaseq_dataset_list = []

        for dsrd in dataset_seq_run_directories:
            meta_data = session.query(SequencingRun.sequencingRunDate, Project.projectDirectory).\
                filter(Project.projectID == Dataset.projectID).\
                filter(Dataset.datasetID == DatasetSeqRun.datasetID).\
                filter(DatasetSeqRun.sequencingRunID == SequencingRun.sequencingRunID).\
                filter(DatasetSeqRun.datasetSeqRunDirectory == dsrd).all()
            
            if self.get_protocol_name(session, dsrd) == 'RNAseq':
                rnaseq_dataset_list.append(dsrd)
            if self.get_protocol_name(session, dsrd) == 'ChIPseq':
                chipseq_dataset_list.append(dsrd)
            if self.get_protocol_name(session, dsrd) == 'ChBIDseq':
                chipseq_dataset_list.append(dsrd)
            if self.get_protocol_name(session, dsrd) == 'CUTTAGseq':
                cuttagseq_dataset_list.append(dsrd)
            if self.get_protocol_name(session, dsrd) == 'sgRNAseq':
                sgrnaseq_dataset_list.append(dsrd)

            for m in meta_data:
                meta_dict[dsrd] = '  run_date: ' + str(m[0].date()) + '\n' + \
                                  '  project: \"' + m[1] + '\"\n'

        
        try:
            if len(rnaseq_dataset_list) > 0:
                meta_file_path = self.rnaseq_meta_yaml_path
                with open(meta_file_path, 'a') as f:
                    f.write('# ' + str(datetime.now().date()) + '\n')
                    for rnaseq in rnaseq_dataset_list:
                        f.write(rnaseq + ":\n")
                        f.write(meta_dict[rnaseq] + '\n')
                with open(self.deseq2_meta_yaml_path, 'a') as f:
                    f.write('# ' + str(datetime.now().date()) + '\n')
                    for rnaseq in rnaseq_dataset_list:
                        f.write(rnaseq + ":\n")
                        deseq2_meta_dict = self.get_deseq2_meta(session, rnaseq)
                        for k,v in deseq2_meta_dict.items():
                            f.write('  ' + k + ':\n')
                            for i in v:
                                f.write('    ' + i.sampleAttributeTypeName + ': \"' + i.sampleAttributeValue + '\"\n')
                        f.write('\n')

            if len(sgrnaseq_dataset_list) > 0:
                sgrnaseq_meta_file_path = self.sgrnaseq_meta_yaml_path
                with open(sgrnaseq_meta_yaml_path, 'a') as f:
                    f.write('# ' + str(datetime.now().date()) + '\n')
                    for sgrnaseq in sgrnaseq_dataset_list:
                        f.write(sgrnaseq + ":\n")
                        f.write(meta_dict[sgrnaseq] + '\n')
                with open(crispr_meta_yaml_path, 'a') as f:
                    f.write('# ' + str(datetime.now().date()) + '\n')
                    for sgrnaseq in sgrnaseq_dataset_list:
                        crispr_meta_dict = self.get_crispr_meta(session, sgrnaseq)
                        for k,v in crispr_meta_dict.items():
                            f.write('  ' + k + ':\n')
                            for i in v:
                                f.write('    ' + i[0] + ': \"' + i[1] + '\"\n')
                        f.write('\n')

                            

            if len(chipseq_dataset_list) > 0:
                meta_file_path = self.chipseq_meta_yaml_path
                chipseq_input_sample_yaml_path = self.chipseq_input_sample_yaml_path
                tornado_meta_yaml_path = self.tornado_meta_yaml_path
                with open(meta_file_path, 'a') as f:
                    f.write('# ' + str(datetime.now().date()) + '\n')
                    for chipseq in chipseq_dataset_list:
                        f.write(chipseq + ":\n")
                        f.write(meta_dict[chipseq] + '\n')

                with open(chipseq_input_sample_yaml_path, 'a') as f:
                    f.write('# ' + str(datetime.now().date()) + '\n')
                    for chipseq in chipseq_dataset_list:
                        f.write(chipseq + ":\n")
                        chipseq_input_sample_data = self.chipseq_input_sample(session, chipseq)
                        for k,v in chipseq_input_sample_data.items():
                            f.write(' ' + k + ':\n')
                            f.write(v + '\n')
                        f.write('\n')
                
                with open(tornado_meta_yaml_path, 'a') as f:
                    f.write('# ' + str(datetime.now().date()) + '\n')
                    for chipseq in chipseq_dataset_list:
                        f.write(chipseq + ':\n')
                        tornado_meta_dict = self.get_tornado_meta(session, chipseq)
                        for k,v in tornado_meta_dict.items():
                            f.write('  ' + k + ':\n')
                            for i in v:
                                f.write('    ' + i.sampleAttributeTypeName + ': \"' + i.sampleAttributeValue + '\"\n')
                        f.write('\n')



            if len(cuttagseq_dataset_list) > 0:
                meta_file_path = self.cuttag_meta_yaml_path
                with open(meta_file_path, 'a') as f:
                    f.write('# ' + str(datetime.now().date()) + '\n')
                    for cuttagseq in cuttagseq_dataset_list:
                        f.write(cuttagseq + ":\n")
                        f.write(meta_dict[cuttagseq] + '\n')
        except:
            raise
        
    # getting chipseq input and samples into a dictionary
    def chipseq_input_sample(self, session, dataset_seq_run_directory):
        chipseq_input_sample_dict = defaultdict()
        chipseq_input_sample_data = session.query(\
            Sample.sampleCPCTName, DatasetSeqRunSample.fastqSampleNumber, \
            DatasetSeqRunSample.chIPInputSampleFlag, \
            DatasetSeqRunSample.sampleID, DatasetSeqRunSample.chIPInputSampleID).\
            filter(DatasetSeqRun.datasetSeqRunID == DatasetSeqRunSample.datasetSeqRunID).\
            filter(DatasetSeqRunSample.sampleID == Sample.sampleID).\
            filter(DatasetSeqRun.datasetSeqRunDirectory == dataset_seq_run_directory).all()

        input_dict = defaultdict()

        for i in chipseq_input_sample_data:
            # CSH 07/31/22 - Changed logic from "if chIPInputSampleID == None"
            # to "if chIPInputSampleFlag == "Y"" -- logic is more direct and
            # custom DB setups are easier with the new logic

            # if (i.chIPInputSampleID == None):
            if (i.chIPInputSampleFlag == "Y"):
                input_dict[i.sampleID] = i.sampleCPCTName + '_S' + str(i.fastqSampleNumber)

        for i in chipseq_input_sample_data:
            if (i.chIPInputSampleID != None):
                chipseq_input_sample_dict[i.sampleCPCTName + '_S' + str(i.fastqSampleNumber)] = '  \'' + input_dict[i.chIPInputSampleID] + '\''

        return(chipseq_input_sample_dict)
        # print(chipseq_input_sample_dict)

    def get_crispr_meta(self, session, dataset_seq_run_directory):
        crispr_meta_dict = defaultdict()
        bcl_meta = BclMetaParser()
        bcl_meta.set_dataset_path(dataset_seq_run_directory)
        sample_ids = session.query(Sample.sampleID).\
                    filter(DatasetSeqRun.datasetSeqRunID == DatasetSeqRunSample.datasetSeqRunID).\
                    filter(DatasetSeqRunSample.sampleID == Sample.sampleID).\
                    filter(DatasetSeqRun.datasetSeqRunDirectory == dataset_seq_run_directory).all()
        for s in sample_ids:
            sample_name = session.query(Sample.sampleCPCTName, DatasetSeqRunSample.fastqSampleNumber).\
                            filter(DatasetSeqRunSample.sampleID == Sample.sampleID).\
                                filter(Sample.sampleID == s.sampleID).first()

            crispr_meta_dict[sample_name.sampleCPCTName + '_S' + str(sample_name.fastqSampleNumber)] = session.query(\
                SampleAttributeType.sampleAttributeTypeName, SampleAttribute.sampleAttributeValue). \
                    filter(SampleAttribute.sampleAttributeTypeID == SampleAttributeType.sampleAttributeTypeID).\
                        filter(DatasetSeqRunSample.sampleID == Sample.sampleID).\
                            filter(SampleAttribute.sampleID == Sample.sampleID).\
                                filter(SampleAttribute.sampleID == s.sampleID).all() 
            sample_path_fullname = bcl_meta.bcl_output_dir + 'output/' + bcl_meta.raw_bcl_path_name + 'merged/' + sample_name.sampleCPCTName + '_S' + str(sample_name.fastqSampleNumber) + '_R1_merged.fastq.gz'
            crispr_meta_dict[sample_name.sampleCPCTName + '_S' + str(sample_name.fastqSampleNumber)].append(('path', sample_path_fullname))
            
        return(crispr_meta_dict)


    def get_deseq2_meta(self, session, dataset_seq_run_directory):
        deseq2_meta_dict = defaultdict()
        sample_ids = session.query(Sample.sampleID).\
                    filter(DatasetSeqRun.datasetSeqRunID == DatasetSeqRunSample.datasetSeqRunID).\
                    filter(DatasetSeqRunSample.sampleID == Sample.sampleID).\
                    filter(DatasetSeqRun.datasetSeqRunDirectory == dataset_seq_run_directory).all()
        for s in sample_ids:
            sample_name = session.query(Sample.sampleCPCTName, DatasetSeqRunSample.fastqSampleNumber).\
                            filter(DatasetSeqRunSample.sampleID == Sample.sampleID).\
                                filter(Sample.sampleID == s.sampleID).first()

            deseq2_meta_dict[sample_name.sampleCPCTName + '_S' + str(sample_name.fastqSampleNumber)] = session.query(\
                SampleAttributeType.sampleAttributeTypeName, SampleAttribute.sampleAttributeValue). \
                    filter(SampleAttribute.sampleAttributeTypeID == SampleAttributeType.sampleAttributeTypeID).\
                        filter(DatasetSeqRunSample.sampleID == Sample.sampleID).\
                            filter(SampleAttribute.sampleID == Sample.sampleID).\
                                filter(SampleAttribute.sampleID == s.sampleID).all()

        return(deseq2_meta_dict)


    def get_tornado_meta(self, session, dataset_seq_run_directory):
        tornado_meta_dict = defaultdict()
        sample_ids = session.query(Sample.sampleID).\
                    filter(DatasetSeqRun.datasetSeqRunID == DatasetSeqRunSample.datasetSeqRunID).\
                    filter(DatasetSeqRunSample.sampleID == Sample.sampleID).\
                    filter(DatasetSeqRun.datasetSeqRunDirectory == dataset_seq_run_directory).all()
        for s in sample_ids:
            sample_name = session.query(Sample.sampleCPCTName, DatasetSeqRunSample.fastqSampleNumber).\
                            filter(DatasetSeqRunSample.sampleID == Sample.sampleID).\
                                filter(Sample.sampleID == s.sampleID).first()

            tornado_meta_dict[sample_name.sampleCPCTName + '_S' + str(sample_name.fastqSampleNumber)] = session.query(\
                SampleAttributeType.sampleAttributeTypeName, SampleAttribute.sampleAttributeValue). \
                    filter(SampleAttribute.sampleAttributeTypeID == SampleAttributeType.sampleAttributeTypeID).\
                        filter(DatasetSeqRunSample.sampleID == Sample.sampleID).\
                            filter(SampleAttribute.sampleID == Sample.sampleID).\
                                filter(SampleAttribute.sampleID == s.sampleID).all()

        return(tornado_meta_dict)



    def get_protocol_name(self, session, dataset_seq_run_directory):
        protocol_name = ''
        protocol_name = session.query(Protocol.protocolName).\
            filter(Protocol.protocolID == Dataset.protocolID).\
            filter(Dataset.datasetID == DatasetSeqRun.datasetID).\
            filter(DatasetSeqRun.datasetSeqRunDirectory == dataset_seq_run_directory).first()

        return(protocol_name[0])

    def get_raw_bcl_path(self, dataset_date):
        file_list = glob('/mnt/storage/dept/pedonc/CPCT/projects/bcl2fastq_demultiplex/*/')
        file_list_dates = [d.split('/')[-2].split('_')[0] for d in file_list]
        raw_bcl_path = ''

        try:
            match_index = file_list_dates.index(dataset_date)
            raw_bcl_path = file_list[match_index].split('/')[-2] + '/'
            return(raw_bcl_path)
        except ValueError:
            return '/'


def main():
    bmw = PipelineMetaWriter()
    with bmw.db_connection.create_session() as session:
        #print(bmw.get_sgrnaseq_meta(session, 'sgRNAseq_220622_NO'))
        bmw.get_crispr_meta(session, 'sgRNAseq_test')
        crispr_meta_file_path = bmw.crispr_meta_yaml_path
        print(crispr_meta_file_path)
        #print(bmw.get_crispr_meta(session, 'sgRNAseq_220622_NO'))
        crispr_meta_dict = bmw.get_crispr_meta(session, 'sgRNAseq_test')
        with open(crispr_meta_file_path, 'a') as f:
            f.write('# ' + str(datetime.now().date()) + '\n')
            f.write('sgRNAseq_test:\n')
            for k,v in crispr_meta_dict.items():
                f.write('  ' + k + ':\n')
                for i in v:
                    f.write('    ' + i[0] + ': \"' + i[1] + '\"\n')
            
       

    #bmw.write_pipeline_yamls()
    #print(bmw.get_raw_bcl_path('220513'))
    #print(bmw.get_tornado_meta('ChIPseq_210820'))

if __name__ == "__main__":
    main()

