# This class configs all the yaml files' path into an Enum
from enum import Enum

class ArgosTools(Enum):
    ARGOS_BIN_PATH = 'argos_bin_path' # Argos conda pipeline envs bin path -- root path for all the bin files
    ARGOS_BIN_PATH_SNAKEMAKE = 'argos_bin_path_snakemake' # Argos conda snakemake envs bin path -- root path for all the bin files
    ARGOS_NGSPLOT_PATH = 'argos_ngsplot_path' # ngsplot tool isn't in pipeline or snakemake env

    COMPUTEMATRIX = 'computeMatrix'
    PLOTHEATMAP = 'plotHeatmap'
    BAMCOVERAGE = 'bamCoverage'
    HTSEQ_COUNT = 'htseq_count'
    BCL2FASTQ = 'bcl2fastq'
    BEDTOOLS = 'bedtools'
    IGVTOOLS = 'igvtools'
    SAMTOOLS = 'samtools'
    MULTIQC = 'multiqc'
    NGSPLOT = 'ngsplot'
    RSCRIPT = 'Rscript'
    PICARD = 'picard'
    FASTQC = 'fastqc'
    MACS2 = 'macs2'
    STAR = 'STAR'

    ARGOS_CONFIG_PATH = 'metaconfig/argos_config.yaml'


class BclMeta(Enum):
    RUN_DATE = 'run_date'
    RAW_BCL_PATH = 'raw_bcl_path'
    FASTQ_PATH = 'fastq_path'
    SAMPLE_SHEET_NAME = 'sample_sheet_name'

    BCL_YAML_PATH = 'metaconfig/bcl_meta.yaml'


class DirConfig(Enum):
    DATASET_BCL_DEMULTIPLEX = 'dataset_bcl_demultiplex'
    DATASET_BCL2FASTQ = 'dataset_bcl2fastq'
    ROOT_PROJECT_PATH = 'root_project_path'
    REFERENCE_PATH = 'reference_path'
    GENE_SPANS_PATH = 'gene_spans_path'
    GENE_LISTS_PATH = 'gene_lists_path'

    DIR_CONFIG_PATH = 'metaconfig/dir_config.yaml'

class BedtoolSummaryConfig(Enum):
    BEDTOOL_SUMMARY_CONFIG_YAML_PATH = 'metaconfig/bedtool_summary_config.yaml'

class BedtoolSummaryColorConfig(Enum):
    BEDTOOL_SUMMARY_COLOR_CONFIG_YAML_PATH = 'metaconfig/bedtool_summary_color_config.yaml'

class SampleMeta(Enum):
    SAMPLE_YAML_PATH = 'metaconfig/sample.yaml'


class RnaseqMeta(Enum):
    RUN_DATE = 'run_date'
    PROJECT = 'project'

    RNA_YAML_PATH = 'metaconfig/rnaseq_meta.yaml'

class DatasetMeta(Enum):
    RUN_DATE = 'run_date'
    PROJECT = 'project'

    DATASET_YAML_PATH = 'metaconfig/dataset_meta.yaml'

class ReferenceMeta(Enum):
    DIR = 'dir'
    STARINDEX = 'starindex'
    STARINDEX_ERCC = 'starindex-ercc'
    GTF = 'gtf'
    FA = 'fa'
    MACS2_REF_GENOME = 'macs2_ref_genome'
    BED = 'bed'
    TXT = 'txt'
    CHRNAMELENGTH = 'chrNameLength'
    HEPI_SGRNA_LIBRARY = 'hEpi_sgRNA_library'

    REFERENCE_YAML_PATH = 'metaconfig/reference_meta.yaml'

class GenomeAnnotationsMeta(Enum):
    GENE_SPANS = 'gene_spans'
    TSS_SPANS = 'tss_spans'
    GENE_LISTS = 'gene_lists'
    BLACKLIST = 'blacklist'
    ENHANCERS = 'enhancers'
    DIR = 'dir'

    YAML_PATH = 'metaconfig/genome_annotations.yaml'

class ChipseqMeta(Enum):
    RUN_DATE = 'run_date'
    PROJECT = 'project'
    FASTQ_PATH = 'fastq_path'

    CHIPSEQ_YAML_PATH = 'metaconfig/chipseq_meta.yaml'


class ChipseqInputSampleMeta(Enum):
    CHIPSEQ_INPUT_SAMPLE_YAML_PATH = 'metaconfig/chipseq_input_sample.yaml'


class ChipseqInputSampleIndexMeta(Enum):
    CHIPSEQ_INPUT_SAMPLE_INDEX_YAML_PATH = 'metaconfig/chipseq_input_sample_index.yaml'

class ChipseqDownsampleMeta(Enum):
    CHIPSEQ_DOWNSAMPLE_META_YAML_PATH = 'metaconfig/chipseq_downsample_meta.yaml'

class CuttagMeta(Enum):
    RUN_DATE = 'run_date'
    PROJECT = 'project'
    CUTTAG_YAML_PATH = 'metaconfig/cuttag_meta.yaml'

class DESeq2Meta(Enum):
    CELLLINE= 'cellline'
    TREATMENT= 'treatment'
    DESEQ2_YAML_PATH = 'metaconfig/deseq2_meta.yaml'
    
class DESeq2StrategyMeta(Enum):
    DESEQ2_STRATEGY_META_YAML_PATH = 'metaconfig/deseq2_strategy.yaml'

class SgRNAseqMeta(Enum):
    RUN_DATE = 'run_date'
    PROJECT = 'project'
    SGRNA_YAML_PATH = 'metaconfig/sgrnaseq_meta.yaml'


class CRISPRMeta(Enum):
    CRISPR_YAML_PATH = 'metaconfig/crispr_meta.yaml'
    
class CRISPRStrategyMeta(Enum):
    CRISPR_STRATEGY_META_YAML_PATH = 'metaconfig/crispr_strategy.yaml'

class TornadoMeta(Enum):
    TORNADO_META_YAML_PATH = 'metaconfig/tornado_meta.yaml'

class TornadoStrategyMeta(Enum):
    TORNADO_STRATEGY_META_YAML_PATH = 'metaconfig/tornado_strategy.yaml'
