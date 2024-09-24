from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, String, DateTime, ForeignKey
from sqlalchemy.orm import relationship
from db.dbconnection import DBconnector

Base = declarative_base()

class AnalysisParameter(Base):
    __tablename__ = 'AnalysisParameter'

    analysisParameterID = Column('AnalysisParameterID', Integer, primary_key = True)
    analysisParameterName = Column(String)
    analysisParameterTypeID = ''
    analysisParameterCode = Column(String)

'''
class AnalysisParameterType(Base):
    __tablename__ = 'AnalysisParameterType'


class AnalysisStep(Base):
    __tablename__ = 'AnalysisStep'



class AnalysisStepParameter(Base):
    __tablename__ = 'AnalysisStepParameter'



class CellLine(Base):
    __tablename__ = 'CellLine'


class ChIPAntibody(Base):
    __tablename__ = 'ChIPAntibody'
    


class ChIPTarget(Base):
    __tablename__ = 'ChIPTarget'


class CohortComparison(Base):
    __tablename__ = 'CohortComparison'



class Compound(Base):
    __tablename__ = 'Compound'
    
'''

class Dataset(Base):
    __tablename__ = 'Dataset'

    datasetID = Column('DatasetID', Integer, primary_key = True)
    datasetName = Column('DatasetName', String)
    datasetStatus = Column('DatasetStatus', String)
    projectID = Column('ProjectID', Integer, ForeignKey("Project.ProjectID"))
    sourceID = Column('SourceID', Integer, ForeignKey("Source.SourceID"))
    protocolID = Column('ProtocolID', Integer, ForeignKey("Protocal.ProtocalID"))
    technologyID = Column('TechnologyID', Integer, ForeignKey("Technology.TechnologyID"))


'''
class DatasetAnalysis(Base):
    __tablename__ = 'DatasetAnalysis'


class DatasetAnalysisStep(Base):
    __tablename__ = 'DatasetAnalysisStep'


class DatasetAnalysisStepParameter(Base):
    __tablename__ = 'DatasetAnalysisStepParameter'

'''

class DatasetSample(Base):
    __tablename__ = 'DatasetSample'

    datasetSampleID = Column('DatasetSampleID', Integer, primary_key = True)
    datasetID = Column('DatasetID', Integer, ForeignKey("Dataset.DatasetID"))
    sampleID = Column('SampleID', Integer, ForeignKey("Sample.SampleID"))

# bcl_meta
class DatasetSeqRun(Base):
    __tablename__ = 'DatasetSeqRun'

    datasetSeqRunID = Column('DatasetSeqRunID', Integer, primary_key = True)
    datasetSeqRunDirectory = Column('DatasetSeqRunDirectory', String)
    datasetID = Column('DatasetID', Integer, ForeignKey("Dataset.DatasetID"))
    sequencingRunID = Column('SequencingRunID', Integer, ForeignKey("SequencingRun.SequencingRunID"))
    scientistID = Column('ScientistID', Integer, ForeignKey("Scientist.ScientistID"))
    bclYamlFlag = Column('BclYamlFlag', String)
    protocolID = Column('ProtocolID', Integer, ForeignKey("Protocol.ProtocolID"))
    


class DatasetSeqRunSample(Base):
    __tablename__ = 'DatasetSeqRunSample'

    datasetSeqRunSampleID = Column('DatasetSeqRunSampleID', Integer, primary_key = True)
    datasetSeqRunSampleIndex = Column('DatasetSeqRunSampleIndex', Integer)
    seqRunSampleBarcode = Column('SeqRunSampleBarcode', String)
    chIPInputSampleFlag = Column('ChIPInputSampleFlag', String)
    fastqSampleNumber = Column('FastqSampleNumber', Integer)
    sampleID = Column('SampleID', Integer, ForeignKey("Sample.SampleID"))
    datasetSeqRunID = Column('DatasetSeqRunID', Integer, ForeignKey("DatasetSeqRun.DatasetSeqRunID"))
    chIPInputSampleID = Column('ChIPInputSampleID', Integer, ForeignKey("Sample.SampleID"))
    sampleGroupID = Column('SampleGroupID', Integer, ForeignKey("SampleGroup.SampleGroupID"))
    sampleGroupPosition = Column('SampleGroupPosition', Integer)
    

'''
class Genome(Base):
    __tablename__ = 'Genome'


class GuideRna(Base):
    __tablename__ = 'GuideRna'


class Lab(Base):
    __tablename__ = 'Lab'


class Modification(Base):
    __tablename__ = 'Modification'


class Pipeline(Base):
    __tablename__ = 'Pipeline'


class PipelineStep(Base):
    __tablename__ = 'PipelineStep'



'''

class Project(Base):
    __tablename__ = 'Project'

    projectID = Column('ProjectID', Integer, primary_key = True)
    projectName = Column('ProjectName', String, unique = True)
    projectStartDate = Column('ProjectStartDate', DateTime)
    projectCode = Column('ProjectCode', String, unique = True)
    projectDirectory = Column('ProjectDirectory', String)
    scientistID = Column('ScientistID', Integer, ForeignKey("Scientist.ScientistID"))



class Protocol(Base):
    __tablename__ = 'Protocol'

    protocolID = Column('ProtocolID', Integer, primary_key = True)    
    protocolName = Column('ProtocolName', String, unique = True)
    protocolPipelineType = Column('ProtocolPipelineType', String)

class Sample(Base):
    __tablename__ = 'Sample'

    sampleID = Column('SampleID', Integer, primary_key = True)
    sampleCPCTName = Column('SampleCPCTName', String, unique = True)
    sampleName = Column('SampleName', String)
    sampleReplicateNumber = Column('SampleReplicateNumber', Integer)
    geneticBackground = Column('GeneticBackground', String)
    spikeinGenomeCode = Column('SpikeinGenomeCode', String)
    sourceID = Column('SourceID', Integer, ForeignKey('Source.SourceID'))
    speciesID = Column('SpeciesID', Integer, ForeignKey('Species.SpeciesID'))
    sampleTypeID = Column('SampleTypeID', Integer, ForeignKey('SampleType.SampleTypeID'))
    sampleTreatmentID = Column('SampleTreatmentID', Integer, ForeignKey('SampleTreatment.SampleTreatmentID'))
    cellLineID = Column('CellLineID', Integer, ForeignKey('CellLine.CellLineID'))
    sampleChIPID = Column('SampleChIPID', Integer, ForeignKey('SampleChIP.SampleChIPID'))
    protocolID = Column('ProtocolID', Integer, ForeignKey('Protocol.ProtocolID'))

class SampleAttribute(Base):
    __tablename__ = "SampleAttribute"

    sampleAttributeID = Column('SampleAttributeID', Integer, primary_key = True)
    sampleAttributeValue = Column('SampleAttributeValue', String)
    sampleID = Column('SampleID', Integer, ForeignKey('Sample.SampleID'))
    sampleAttributeTypeID = Column('SampleAttributeTypeID', Integer, ForeignKey('SampleAttributeType.SampleAttributeTypeID'))


class SampleAttributeType(Base):
    __tablename__ = "SampleAttributeType"

    sampleAttributeTypeID = Column('SampleAttributeTypeID', Integer, primary_key = True)
    sampleAttributeTypeName = Column('SampleAttributeTypeName', String)

class SampleGroup(Base):
    __tablename__ = "SampleGroup"

    sampleGroupID = Column('SampleGroupID', Integer, primary_key = True)
    sampleGroupNumber = Column('SampleGroupNumber', Integer)
    sampleGroupBaseColor = Column('SampleGroupBaseColor', String)

'''
class SampleChIP(Base):
    __tablename__ = 'SampleChIP'


class SampleCohort(Base):
    __tablename__ = 'SampleCohort'


class SampleCohortSample(Base):
    __tablename__ = 'SampleCohortSample'


class SampleModification(Base):
    __tablename__ = 'SampleModification'


class SampleTreatment(Base):
    __tablename__ = 'SampleTreatment'


class SampleType(Base):
    __tablename__ = 'SampleType'




class SeqRunProtocol(Base):
    __tablename__ = 'SeqRunProtocol'


'''


class Scientist(Base):
    __tablename__ = 'Scientist'

    scientistID = Column('ScientistID', Integer, primary_key = True)
    scientistFirstName = Column('ScientistFirstName', String)
    scientistLastName = Column('ScientistLastName', String)
    scientistInitials = Column('ScientistInitials', String)
    labID = Column('LabID', Integer, ForeignKey("Lab.LabID"))
    hashcode = Column('Hashcode', String)
    username = Column('Username', String)
    adminFlag = Column('AdminFlag', String)
    RoleID = Column('RoleID', Integer, ForeignKey("Role.RoleID"))

class SequencingRun(Base):
    __tablename__ = 'SequencingRun'

    sequencingRunID = Column('SequencingRunID', Integer, primary_key = True)
    sequencingRunCode = Column('SequencingRunCode', String, unique = True)
    sequencingRunDate = Column('SequencingRunDate', DateTime)
    sequencingRunFormat = Column('SequencingRunFormat', String)
    sequencingRunReadLength = Column('SequencingRunReadLength', Integer)
    technologyID = Column('TechnologyID', Integer, ForeignKey("Technology.TechnologyID"))


'''
class Source(Base):
    __tablename__ = 'Source'


class Species(Base):
    __tablename__ = 'Species'


class Technology(Base):
    __tablename__ = 'Technology'


class TreatmentCompound(Base):
    __tablename__ = 'TreatmentCompound'


class TreatmentGuideRna(Base):
    __tablename__ = 'TreatmentGuideRna'


class Vendor(Base):
    __tablename__ = 'Vendor'

'''
