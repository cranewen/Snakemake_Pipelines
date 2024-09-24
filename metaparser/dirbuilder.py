# A class for pre-building all the directories for a snakemake run
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from metaparser.yamlhandle import YamlHandle
from metaparser.metaconf import *
from metaparser.sampleparser import SampleParser
from metaparser.genomeannotationsparser import GenomeAnnotationsParser


class DirBuilder:
    def __init__(self, dataset_name):
        # bcl2fastq pipeline varibles
        self.data_dir = '' # where the raw data/demultiplex data at
        self.output_root_dir = '' # the root dir for the outputs, usually its bcl2fastq/
        self.bcl_demultiplex_dir = '' # raw demultiplex directory
        self.dataset_name = dataset_name # Using the dataset_name through out the entire pipeline
        self.dataset_dir = '' # an output folder for each dataset, a subdiretory of output_root_dir
        self.fastq_output_dir = '' # just a output/ subdirectory of dataset_dir
        self.merged_fastq_dir = '' # for merged fastq files in fastq_output_dir
        self.fastqc_dir = '' # a fastqc directory in fastq_output_dir
        self.multiqc_dir = '' # a multiqc directory in fastq_output_dir
        self.sample_sheet = '' # sample sheet name and path

        ### RNAseq and ChIPseq pipeline varibles ###
        # RNAseq pipeline varibles
        self.root_project_path = ''
        self.reference_path = ''
        self.fastq_dir = ''
        self.star_dir = ''
        self.star_dataset_dir = ''
        self.project_path = ''
        self.starindex_ercc = ''
        self.sorted_bam_dir = ''
        self.sorted_bam_dataset_dir = ''
        self.rmdup_bam_dir = ''
        self.rmdup_bam_dataset_dir = ''
        self.counts_rmdup_dir = ''
        self.counts_sorted_dir = ''
        self.tdf_dataset_rmdup_dir = ''
        self.tdf_dataset_sorted_dir = ''
        self.wig_dataset_rmdup_dir = ''
        self.wig_dataset_sorted_dir = ''
        self.gtf = ''

        # ChIPseq pipeline varibles
        # Including all the RNAseq variables with additional varibles
        self.macs2_dir = ''
        self.macs2_dataset_dir = ''
        self.bedtools_coverage_dir = ''
        self.bedtools_coverage_dataset_rmdup_dir = ''

        self.bedtools_chrNameLength = ''
        self.bedtools_tss_1k_2kbuff_bed_ref = ''

        self.bedtools_intersect_dir = ''
        self.bedtools_intersect_input_dir = ''
        self.bedtools_intersect_dataset_rmdup_dir = ''
        
        # deeptools_bamCoverage varibles, it's the first step of deeptools for creating tornado plot
        self.bigwig_dir = ''
        self.bigwig_sorted_dataset_dir = ''
        self.bigwig_rmdup_dataset_dir = ''

        # deeptools_computeMatrix varibles (2nd step of deeptools)
        self.deeptools_dir = ''
        self.matrix_dir = ''
        self.matrix_sorted_dataset_dir = ''
        self.matrix_rmdup_dataset_dir = ''
        self.heatmap_dir = ''
        self.heatmap_sorted_dataset_dir = ''
        self.heatmap_rmdup_dataset_dir = ''
        
        # CUTTAG pipeline varibles
        self.macs2_dataset_rmdup_dir = ''
        self.macs2_dataset_sorted_dir = ''
        #self.bedtools_coverage_dataset_rmdup_dir = ''
        #self.bedtools_coverage_dataset_sorted_dir = ''
        # ChIPseq and CUTTAG shared varibles
        self.macs2_ref_genome = ''

        # DESeq2 pipeline varibles
        self.deseq2_dir = ''
        self.deseq2_dataset = ''
        self.deseq2_rlog_expression = ''
        self.deseq2_diffexp_genes = ''
        self.pca_dir = ''
        self.heatmaps_dir = ''
        self.volcano_plot_dir = ''
        self.deseq2_R = ''

        # tornado plot pipeline varibles, it's based on ChIPseq, so adding all the 
        # variables created in ChIPseq pipeline
        self.tornado_dir = ''
        self.tornado_dataset_dir = ''
        self.tornado_sorted_dataset_dir = ''
        self.tornado_rmdup_dataset_dir = ''
        #self.tornado_config_dir = ''
        self.tornado_sorted_config_dir = ''
        self.tornado_rmdup_config_dir = ''
        self.tornado_bam_sorted_dir = ''
        self.tornado_bam_rmdup_dir = ''

        
        # global multiqc varibles
        # multiqc root path
        self.multiqc_all_root = ''
        # multiqc dataset path
        self.multiqc_all_dir = ''

    def __repr__(self):
        return f'{self.__dict__}'

    # for CRISPR; sgRNA is different than ChIPseq and rnaseq, because there is only crispr analysis step
    # therefore, there is no need to build_crispr_dir()
    def build_sgrnaseq_dirs(self, ref_name):
        conf_handle = YamlHandle(DirConfig.DIR_CONFIG_PATH.value).read_yaml()
        for c in conf_handle:
            self.root_project_path = c[DirConfig.ROOT_PROJECT_PATH.value]

        sgrna_meta_handle = YamlHandle(SgRNAseqMeta.SGRNA_YAML_PATH.value).read_yaml()
        for s in sgrna_meta_handle:
            self.project_path = self.root_project_path + s[self.dataset_name][SgRNAseqMeta.PROJECT.value] + '/'

        self.mageck_vispr_dir = self.project_path + 'mageck-vispr/'
        self.crispr_input_dir = self.mageck_vispr_dir + self.dataset_name + '/input/'  
        self.mvispr_output_dir = self.mageck_vispr_dir + self.dataset_name + '/mvispr_' + self.dataset_name + '/'

    def build_rnaseq_dirs(self, ref_name):
        conf_handle = YamlHandle(DirConfig.DIR_CONFIG_PATH.value).read_yaml()
        for c in conf_handle:
            self.root_project_path = c[DirConfig.ROOT_PROJECT_PATH.value]
            self.reference_path = c[DirConfig.REFERENCE_PATH.value]

        rna_meta_handle = YamlHandle(RnaseqMeta.RNA_YAML_PATH.value).read_yaml()
        for r in rna_meta_handle:
            self.project_path = self.root_project_path + r[self.dataset_name][RnaseqMeta.PROJECT.value] + '/'
        self.fastq_dir = self.project_path + 'FASTQ/'
        self.fastq_dataset_dir = self.fastq_dir + self.dataset_name + '/'
        self.star_dir = self.project_path + 'STAR/'
        self.star_dataset_dir = self.star_dir + self.dataset_name + '_' + ref_name + '/'
        self.sorted_bam_dir = self.project_path + 'BAM_sorted/'
        self.sorted_bam_dataset_dir = self.sorted_bam_dir + self.dataset_name + '_' + ref_name + '/'
        self.rmdup_bam_dir = self.project_path + 'BAM_rmdup/'
        self.rmdup_bam_dataset_dir = self.rmdup_bam_dir + self.dataset_name + '_' + ref_name + '/'
        self.counts_rmdup_dir = self.project_path + 'counts/' + self.dataset_name + '_' + ref_name + '/sorted_rmdup'
        self.counts_sorted_dir = self.project_path + 'counts/' + self.dataset_name + '_' + ref_name + '/sorted'
        self.tdf_dataset_rmdup_dir = self.project_path + 'TDF/' + self.dataset_name + '_' + ref_name + '/sorted_rmdup'
        self.tdf_dataset_sorted_dir = self.project_path + 'TDF/' + self.dataset_name + '_' + ref_name + '/sorted'
        self.wig_dataset_rmdup_dir = self.project_path + 'WIG/' + self.dataset_name + '_' + ref_name + '/sorted_rmdup'
        self.wig_dataset_sorted_dir = self.project_path + 'WIG/' + self.dataset_name + '_' + ref_name + '/sorted'

        self.multiqc_all_root = self.project_path + 'multiqc/'
        self.multiqc_all_dir = self.multiqc_all_root + self.dataset_name + '_' + ref_name

        self.sorted_bam_flagstat_dir = self.sorted_bam_dir + self.dataset_name + '_' + ref_name + '/samtools_flagstat/'
        self.rmdup_bam_flagstat_dir = self.rmdup_bam_dir + self.dataset_name + '_' + ref_name + '/samtools_flagstat/'

        reference_meta_handle = YamlHandle(ReferenceMeta.REFERENCE_YAML_PATH.value).read_yaml()
        for ref in reference_meta_handle:
            ref_dir = ref[ref_name][ReferenceMeta.DIR.value]
            self.starindex_ercc = self.reference_path + ref_dir + ref[ref_name][ReferenceMeta.STARINDEX_ERCC.value]
            self.gtf = self.starindex_ercc + '/' + ref[ref_name][ReferenceMeta.GTF.value]


    def build_cuttag_dirs(self, ref_name):
        conf_handle = YamlHandle(DirConfig.DIR_CONFIG_PATH.value).read_yaml()
        for c in conf_handle:
            self.root_project_path = c[DirConfig.ROOT_PROJECT_PATH.value]
            self.reference_path = c[DirConfig.REFERENCE_PATH.value]

        rna_meta_handle = YamlHandle(CuttagMeta.CUTTAG_YAML_PATH.value).read_yaml()
        for r in rna_meta_handle:
            self.project_path = self.root_project_path + r[self.dataset_name][CuttagMeta.PROJECT.value] + '/'
        self.fastq_dir = self.project_path + 'FASTQ/'
        self.fastq_dataset_dir = self.fastq_dir + self.dataset_name + '/'
        self.star_dir = self.project_path + 'STAR/'
        self.star_dataset_dir = self.star_dir + self.dataset_name + '_' + ref_name + '/'
        self.sorted_bam_dir = self.project_path + 'BAM_sorted/'
        self.sorted_bam_dataset_dir = self.sorted_bam_dir + self.dataset_name + '_' + ref_name + '/'
        self.rmdup_bam_dir = self.project_path + 'BAM_rmdup/'
        self.rmdup_bam_dataset_dir = self.rmdup_bam_dir + self.dataset_name + '_' + ref_name + '/'
        self.counts_rmdup_dir = self.project_path + 'counts/' + self.dataset_name + '_' + ref_name + '/sorted_rmdup'
        self.counts_sorted_dir = self.project_path + 'counts/' + self.dataset_name + '_' + ref_name + '/sorted'
        self.tdf_dataset_rmdup_dir = self.project_path + 'TDF/' + self.dataset_name + '_' + ref_name + '/sorted_rmdup'
        self.tdf_dataset_sorted_dir = self.project_path + 'TDF/' + self.dataset_name + '_' + ref_name + '/sorted'
        self.wig_dataset_rmdup_dir = self.project_path + 'WIG/' + self.dataset_name + '_' + ref_name + '/sorted_rmdup'
        self.wig_dataset_sorted_dir = self.project_path + 'WIG/' + self.dataset_name + '_' + ref_name + '/sorted'
        self.macs2_dir = self.project_path + 'MACS2/'
        self.macs2_dataset_rmdup_dir = self.macs2_dir + self.dataset_name + '_' + ref_name
        self.macs2_dataset_sorted_dir = self.macs2_dir + self.dataset_name + '_' + ref_name
        self.bedtools_coverage_dir = self.project_path + 'bedtoos_coverage/'
        self.bedtools_coverage_dataset_rmdup_dir = self.bedtools_coverage_dir + self.dataset_name + '_' + ref_name + '/sorted_rmdup'
        self.bedtools_coverage_dataset_sorted_dir = self.bedtools_coverage_dir + self.dataset_name + '_' + ref_name + '/sorted'

        self.multiqc_all_root = self.project_path + 'multiqc/'
        self.multiqc_all_dir = self.multiqc_all_root + self.dataset_name + '_' + ref_name

        reference_meta_handle = YamlHandle(ReferenceMeta.REFERENCE_YAML_PATH.value).read_yaml()
        for ref in reference_meta_handle:
            ref_dir = ref[ref_name][ReferenceMeta.DIR.value]
            self.starindex_ercc = self.reference_path + ref_dir + ref[ref_name][ReferenceMeta.STARINDEX_ERCC.value]
            self.gtf = self.starindex_ercc + '/' + ref[ref_name][ReferenceMeta.GTF.value]
            self.macs2_ref_genome = ref[ref_name][ReferenceMeta.MACS2_REF_GENOME.value]


    def build_chipseq_dirs(self, ref_name):
        conf_handle = YamlHandle(DirConfig.DIR_CONFIG_PATH.value).read_yaml()
        for c in conf_handle:
            self.root_project_path = c[DirConfig.ROOT_PROJECT_PATH.value]
            self.reference_path = c[DirConfig.REFERENCE_PATH.value]

        bedtool_summary_config_handle = YamlHandle(BedtoolSummaryConfig.BEDTOOL_SUMMARY_CONFIG_YAML_PATH.value).read_yaml()
        self.bedtool_summary_config = None
        try:
            for b in bedtool_summary_config_handle:
                self.bedtool_summary_config = b[self.dataset_name]
        except KeyError:
            print('check config file if the dataset is there')
            

        chipseq_meta_handle = YamlHandle(ChipseqMeta.CHIPSEQ_YAML_PATH.value).read_yaml()
        for c in chipseq_meta_handle:
            self.project_path = self.root_project_path + c[self.dataset_name][ChipseqMeta.PROJECT.value] + '/'

        self.fastq_dir = self.project_path + 'FASTQ/'
        self.fastq_dataset_dir = self.fastq_dir + self.dataset_name + '/'
        self.star_dir = self.project_path + 'STAR/'
        self.star_dataset_dir = self.star_dir + self.dataset_name + '_' + ref_name + '/'
        self.sorted_bam_dir = self.project_path + 'BAM_sorted/'
        self.sorted_bam_dataset_dir = self.sorted_bam_dir + self.dataset_name + '_' + ref_name + '/'
        self.rmdup_bam_dir = self.project_path + 'BAM_rmdup/'
        self.rmdup_bam_dataset_dir = self.rmdup_bam_dir + self.dataset_name + '_' + ref_name + '/'
        self.counts_rmdup_dir = self.project_path + 'counts/' + self.dataset_name + '_' + ref_name + '/sorted_rmdup'
        self.counts_sorted_dir = self.project_path + 'counts/' + self.dataset_name + '_' + ref_name + '/sorted'
        self.tdf_dataset_rmdup_dir = self.project_path + 'TDF/' + self.dataset_name + '_' + ref_name + '/sorted_rmdup'
        self.tdf_dataset_sorted_dir = self.project_path + 'TDF/' + self.dataset_name + '_' + ref_name + '/sorted'
        self.wig_dataset_rmdup_dir = self.project_path + 'WIG/' + self.dataset_name + '_' + ref_name + '/sorted_rmdup'
        self.wig_dataset_sorted_dir = self.project_path + 'WIG/' + self.dataset_name + '_' + ref_name + '/sorted'
        self.macs2_dir = self.project_path + 'MACS2/'
        self.macs2_dataset_dir = self.macs2_dir + self.dataset_name + '_' + ref_name + '/'
        self.bedtools_coverage_dir = self.project_path + 'bedtools_coverage/'
        self.bedtools_coverage_dataset_rmdup_dir = self.bedtools_coverage_dir + self.dataset_name + '_' + ref_name + '/sorted_rmdup' 
        self.bedtools_coverage_dataset_dir = self.bedtools_coverage_dir + self.dataset_name + '_' + ref_name + '/' 

        self.bedtools_intersect_dir = self.macs2_dir
        # need to add nofloor/ in Snakefile_chipseq
        self.bedtools_intersect_input_dir = self.bedtools_intersect_dir + self.dataset_name + '_' + ref_name + '/output_rmdup/nofloor/'
        self.bedtools_intersect_dataset_rmdup_dir = self.bedtools_intersect_dir + self.dataset_name + '_' + ref_name + '/' + 'peak_overlaps_rmdup/'

        self.bigwig_dir = self.project_path + 'BigWig/'
        self.bigwig_sorted_dataset_dir = self.bigwig_dir + self.dataset_name + '_' + ref_name + '/sorted'
        self.bigwig_rmdup_dataset_dir = self.bigwig_dir + self.dataset_name + '_' + ref_name + '/sorted_rmdup'

        self.sorted_bam_flagstat_dir = self.sorted_bam_dir + self.dataset_name + '_' + ref_name + '/samtools_flagstat/'
        self.rmdup_bam_flagstat_dir = self.rmdup_bam_dir + self.dataset_name + '_' + ref_name + '/samtools_flagstat/'

        # multiqc root path
        self.multiqc_all_root = self.project_path + 'multiqc/'
        # multiqc dataset path
        self.multiqc_all_dir = self.multiqc_all_root + self.dataset_name + '_' + ref_name

        reference_meta_handle = YamlHandle(ReferenceMeta.REFERENCE_YAML_PATH.value).read_yaml()
        gap = GenomeAnnotationsParser(ref_name)
        self.bedtools_tss_1k_2kbuff_bed_ref = gap.get_tss_spans('unique_protcoding_tx_TSS+1k_2kbuffer_coords')
        for ref in reference_meta_handle:
            ref_dir = ref[ref_name][ReferenceMeta.DIR.value]
            self.starindex_ercc = self.reference_path + ref_dir + ref[ref_name][ReferenceMeta.STARINDEX_ERCC.value]
            self.gtf = self.starindex_ercc + '/' + ref[ref_name][ReferenceMeta.GTF.value]
            self.macs2_ref_genome = ref[ref_name][ReferenceMeta.MACS2_REF_GENOME.value]
            self.bedtools_chrNameLength = self.starindex_ercc + '/' + ref[ref_name][ReferenceMeta.CHRNAMELENGTH.value]


    def build_chipseq_spikein_dirs(self, ref_name):
        conf_handle = YamlHandle(DirConfig.DIR_CONFIG_PATH.value).read_yaml()
        for c in conf_handle:
            self.root_project_path = c[DirConfig.ROOT_PROJECT_PATH.value]
            self.reference_path = c[DirConfig.REFERENCE_PATH.value]
            self.data_dir = c[DirConfig.DATASET_BCL_DEMULTIPLEX.value]
            self.output_root_dir = self.data_dir + c[DirConfig.DATASET_BCL2FASTQ.value]

        self.dataset_dir = self.output_root_dir + self.dataset_name + '/'

        
        chipseq_meta_handle = YamlHandle(ChipseqMeta.CHIPSEQ_YAML_PATH.value).read_yaml()
        for c in chipseq_meta_handle:
            if ChipseqMeta.FASTQ_PATH.value in c[self.dataset_name].keys():
                self.dataset_name_original = c[self.dataset_name][ChipseqMeta.FASTQ_PATH.value]
                self.project_path = self.root_project_path + c[self.dataset_name_original][ChipseqMeta.PROJECT.value] + '/'
            else:
                self.dataset_name_original = self.dataset_name
                self.project_path = self.root_project_path + c[self.dataset_name_original][ChipseqMeta.PROJECT.value] + '/'

        self.fastq_dir = self.project_path + 'FASTQ/'
        self.fastq_dataset_dir = self.fastq_dir + self.dataset_name_original + '/'
        self.star_dir = self.project_path + 'STAR/'
        self.star_dataset_dir = self.star_dir + self.dataset_name + '_' + ref_name + '/'
        self.sorted_bam_dir = self.project_path + 'BAM_sorted/'
        self.sorted_bam_dataset_dir = self.sorted_bam_dir + self.dataset_name + '_' + ref_name + '/'
        self.rmdup_bam_dir = self.project_path + 'BAM_rmdup/'
        self.rmdup_bam_dataset_dir = self.rmdup_bam_dir + self.dataset_name + '_' + ref_name + '/'
        self.sorted_bam_flagstat_dir = self.sorted_bam_dir + self.dataset_name + '_' + ref_name + '/samtools_flagstat/'
        self.rmdup_bam_flagstat_dir = self.rmdup_bam_dir + self.dataset_name + '_' + ref_name + '/samtools_flagstat/'
        self.bedtools_coverage_dir = self.project_path + 'bedtools_coverage/'
        self.bedtools_coverage_dataset_rmdup_dir = self.bedtools_coverage_dir + self.dataset_name + '_' + ref_name + '/sorted_rmdup' 
        self.bedtools_coverage_dataset_dir = self.bedtools_coverage_dir + self.dataset_name + '_' + ref_name + '/' 

        # multiqc root path
        self.multiqc_all_root = self.project_path + 'multiqc/'
        # multiqc dataset path
        self.multiqc_all_dir = self.multiqc_all_root + self.dataset_name + '_' + ref_name

        reference_meta_handle = YamlHandle(ReferenceMeta.REFERENCE_YAML_PATH.value).read_yaml()

        for ref in reference_meta_handle:
            ref_dir = ref[ref_name][ReferenceMeta.DIR.value]
            self.starindex_ercc = self.reference_path + ref_dir + ref[ref_name][ReferenceMeta.STARINDEX_ERCC.value]



    def build_bclfastq_dirs(self):
        conf_handle = YamlHandle(DirConfig.DIR_CONFIG_PATH.value).read_yaml()
        for c in conf_handle:
            self.data_dir = c[DirConfig.DATASET_BCL_DEMULTIPLEX.value]
            self.output_root_dir = self.data_dir + c[DirConfig.DATASET_BCL2FASTQ.value]

        self.dataset_dir = self.output_root_dir + self.dataset_name + '/'

        bcl_meta_handle = YamlHandle(BclMeta.BCL_YAML_PATH.value).read_yaml()        
        for b in bcl_meta_handle:
            self.fastq_output_dir = self.dataset_dir + 'output/' + b[self.dataset_name][BclMeta.RAW_BCL_PATH.value]
            self.bcl_demultiplex_dir = self.data_dir + b[self.dataset_name][BclMeta.RAW_BCL_PATH.value]
            self.sample_sheet = self.data_dir + b[self.dataset_name][BclMeta.SAMPLE_SHEET_NAME.value]

        self.merged_fastq_dir = self.fastq_output_dir + 'merged/'
        self.fastqc_dir = self.fastq_output_dir + 'fastqc/'
        self.multiqc_dir = self.fastq_output_dir + 'multiqc/'

    # use all the RNAseq dirs and add some extra ones for deseq2
    def build_deseq2_dirs(self, ref_name):
        conf_handle = YamlHandle(DirConfig.DIR_CONFIG_PATH.value).read_yaml()
        for c in conf_handle:
            self.root_project_path = c[DirConfig.ROOT_PROJECT_PATH.value]
            self.reference_path = c[DirConfig.REFERENCE_PATH.value]

        rna_meta_handle = YamlHandle(RnaseqMeta.RNA_YAML_PATH.value).read_yaml()
        for r in rna_meta_handle:
            self.project_path = self.root_project_path + r[self.dataset_name][RnaseqMeta.PROJECT.value] + '/'
        self.counts_rmdup_dir = self.project_path + 'counts/' + self.dataset_name + '_' + ref_name + '/sorted_rmdup'
        self.counts_sorted_dir = self.project_path + 'counts/' + self.dataset_name + '_' + ref_name + '/sorted'
        self.deseq2_dir = self.project_path + 'DESeq2/'
        self.deseq2_dataset_sorted = self.deseq2_dir + self.dataset_name + '_' + ref_name + '/' + 'sorted/'
        self.deseq2_dataset_sorted_rmdup = self.deseq2_dir + self.dataset_name + '_' + ref_name + '/' + 'sorted_rmdup/'
        self.deseq2_rlog_expression = self.deseq2_dataset + 'rlog_expression'
        self.deseq2_diffexp_genes = self.deseq2_dataset + 'diffexp_genes'
        self.pca_dir = self.deseq2_dataset + 'PCA/'
        self.heatmaps_dir = self.deseq2_dataset + 'heatmaps/'
        self.volcano_plot_dir = self.deseq2_dataset + 'volcano_plot/'
        self.deseq2_R = '/mnt/storage/dept/pedonc/src/lab_pipelines/deseq2_r.R'

        reference_meta_handle = YamlHandle(ReferenceMeta.REFERENCE_YAML_PATH.value).read_yaml()
        for ref in reference_meta_handle:
            ref_dir = ref[ref_name][ReferenceMeta.DIR.value]
            self.starindex_ercc = self.reference_path + ref_dir + ref[ref_name][ReferenceMeta.STARINDEX_ERCC.value]
            self.gtf = self.starindex_ercc + '/' + ref[ref_name][ReferenceMeta.GTF.value]

    # use all the ChIPseq dirs and add some extra ones for tornado
    def build_tornado_dirs(self, ref_name):
        conf_handle = YamlHandle(DirConfig.DIR_CONFIG_PATH.value).read_yaml()
        for c in conf_handle:
            self.root_project_path = c[DirConfig.ROOT_PROJECT_PATH.value]
            self.reference_path = c[DirConfig.REFERENCE_PATH.value]
        
        chipseq_meta_handle = YamlHandle(ChipseqMeta.CHIPSEQ_YAML_PATH.value).read_yaml()
        for c in chipseq_meta_handle:
            self.project_path = self.root_project_path + c[self.dataset_name][ChipseqMeta.PROJECT.value] + '/'

        self.tornado_dir = self.project_path + 'ngsplot/'
        self.tornado_bam_sorted_dir = self.project_path + 'BAM_sorted/' + self.dataset_name + '_' + ref_name + '/'
        self.tornado_bam_rmdup_dir = self.project_path + 'BAM_rmdup/' + self.dataset_name + '_' + ref_name + '/'
        self.tornado_sorted_dataset_dir = self.tornado_dir + self.dataset_name + '_' + ref_name + '/sorted/'
        self.tornado_rmdup_dataset_dir = self.tornado_dir + self.dataset_name + '_' + ref_name + '/rmdup/'
        self.tornado_sorted_config_dir = self.tornado_sorted_dataset_dir + 'config/'
        self.tornado_rmdup_config_dir = self.tornado_rmdup_dataset_dir + 'config/'
    
    # deeptools_computeMatrix (this is the 2nd step for producing the tornado plot)
    def build_tornado_deeptools_dirs(self, ref_name):
        conf_handle = YamlHandle(DirConfig.DIR_CONFIG_PATH.value).read_yaml()
        for c in conf_handle:
            self.root_project_path = c[DirConfig.ROOT_PROJECT_PATH.value]
            self.reference_path = c[DirConfig.REFERENCE_PATH.value]

        chipseq_meta_handle = YamlHandle(ChipseqMeta.CHIPSEQ_YAML_PATH.value).read_yaml()
        for c in chipseq_meta_handle:
            self.project_path = self.root_project_path + c[self.dataset_name][ChipseqMeta.PROJECT.value] + '/'

        self.deeptools_dir = self.project_path + 'deeptools/'
        self.matrix_dir = self.deeptools_dir + self.dataset_name + '_' + ref_name + '/matrix/'
        self.matrix_sorted_dataset_dir = self.matrix_dir + 'sorted/'
        self.matrix_rmdup_dataset_dir = self.matrix_dir + 'sorted_rmdup/'
        self.heatmap_dir = self.deeptools_dir + self.dataset_name + '_' + ref_name + '/heatmap/'
        self.heatmap_sorted_dataset_dir = self.heatmap_dir + 'sorted/'
        self.heatmap_rmdup_dataset_dir = self.heatmap_dir + 'sorted_rmdup/'

        

        
        
# Tests
def main():
    dirobj = DirBuilder('ChIPseq_test')
    dirobj.build_chipseq_dirs('hg38')


if __name__ == "__main__":
    main()

