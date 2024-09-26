import subprocess
import argparse
import sys, os
from metaparser.bclmetaparser import BclMetaParser
from metaparser.dirbuilder import DirBuilder


parser = argparse.ArgumentParser(description='User info for snakemake run')
parser.add_argument('-d', '--dataset', type=str, required=True)
parser.add_argument('-f', '--ref_genome', type=str)
parser.add_argument('-s', '--snakefile', type=str)
parser.add_argument('-r', '--dryrun', type=str)
parser.add_argument('-g', '--dag', type=str)
parser.add_argument('-l', '--localrun', type=str)
parser.add_argument('-t', '--datatype', type=str) # datatypes are: rnaseq, chipseq, etc
parser.add_argument('-res', '--resources', type=str)
parser.add_argument('-nl', '--nolock', type=str)
parser.add_argument('-k', '--keep_going', type=str)
parser.add_argument('-ambiguity', '--allow_ambiguity', type=str)
parser.add_argument('-e', '--enhancers', type=str)
parser.add_argument('-gs', '--gene_spans', type=str)
parser.add_argument('-ts', '--tss_spans', type=str)
parser.add_argument('-bl', '--blacklist', type=str)
parser.add_argument('-ds', '--downsample_method', type=str)
parser.add_argument('-p', '--print_shell', type=str)
args = parser.parse_args()


# If non-default snakefile is specified, run the specified snakefile
snakefile = "Snakefile"
if args.snakefile is not None:
    snakefile = args.snakefile

# data_type are R: rnaseq, Ch: chipseq, CT: cut&tag, CR: cut&run
data_type = False
if args.datatype is not None:
    data_type = args.datatype

# If reference genome is specified, capture value to pass as snakemake param
ref_genome = False
if args.ref_genome is not None:
    ref_genome = args.ref_genome

# If local run (i.e., no cluster/qsub call) is specified, set localrun flag
localrun = False
if args.localrun is not None:
    localrun = args.localrun

resources = False
if args.resources is not None:
    resources = args.resources

# If dry run is specified, set dry run flag (superceded by DAG param)
dryrun = False
if args.dryrun is not None:
    dryrun = args.dryrun

# If DAG computation is specified, set DAG flag (supercedes dry run param)
dag = False
if args.dag is not None:
    dag = args.dag

# Set defaults for number of cores (for local run) or jobs (for cluster run);
#  we may later make these parameters that can override the defaults
cores = 4
jobs = 60
# Using nodes to control the max number of a batch jobs
nodes = 60

lock = False
if args.nolock is not None:
    lock = args.nolock

keep_going = False
if args.keep_going is not None:
    keep_going = args.keep_going

# set ambiguity for cuttage
ambiguity = False
if args.allow_ambiguity is not None:
    ambiguity = args.allow_ambiguity

# ChIPseq bedtools_intersect arguments
enhancers = False
if args.enhancers is not None:
    enhancers = args.enhancers

gene_spans = False
if args.gene_spans is not None:
    gene_spans = args.gene_spans

tss_spans = False
if args.tss_spans is not None:
    tss_spans = args.tss_spans

blacklist = False
if args.blacklist is not None:
    blacklist = args.blacklist

ds_method = False
if args.downsample_method is not None:
    ds_method = args.downsample_method

print_shell = False
if args.print_shell is not None:
    print_shell = args.print_shell

workingdir = ''

def runsnakemake(dataset_name):
    ###
    # Prepare and submit the snakemake command to run the pipeline
    ###

    # Setting log_dir 
    #   If no ref genome specified, set dataset dir in bcl2fastq project;
    #     otherwise, set dataset dir in project associated with dataset

    # Notes on deseq2 snakemake pipeline command:
    # e.g. snakemake --snakefile Snakefile_deseq2 --config dataset=RNAseq_210201 ref_genome=hg38 --cluster qsub -j 60
    dirbuilder = DirBuilder(dataset_name)

    if ref_genome:

        # NOTE: Using build_rnaseq_dirs(), but lookup should work for any 
        #   project-based (i.e., non-bcl2fastq) dataset -- should we do this
        #   in a more specific / explicit way?
        if data_type == 'R':
            dirbuilder.build_rnaseq_dirs(ref_genome)
            project_path = dirbuilder.project_path
            workingdir = project_path + 'logs' + '/' + dataset_name
            os.makedirs(workingdir, mode=0o777, exist_ok=True)
        if data_type == 'Ch':
            dirbuilder.build_chipseq_dirs(ref_genome)
            project_path = dirbuilder.project_path
            workingdir = project_path + 'logs' + '/' + dataset_name
            os.makedirs(workingdir, mode=0o777, exist_ok=True)
        if data_type == 'Chds':
            dirbuilder.build_chipseq_dirs(ref_genome)
            project_path = dirbuilder.project_path
            workingdir = project_path + 'logs' + '/' + dataset_name + '/downsample'
            os.makedirs(workingdir, mode=0o777, exist_ok=True)
        if data_type == 'ChSI':
            dirbuilder.build_chipseq_spikein_dirs(ref_genome)
            project_path = dirbuilder.project_path
            workingdir = project_path + 'logs' + '/' + dataset_name
            if not 'dm6' in dataset_name:
                workingdir = project_path + 'logs' + '/' + dataset_name + '_dm6'
            else:
                workingdir = workingdir
            os.makedirs(workingdir, mode=0o777, exist_ok=True)
        if data_type == 'CT':
            dirbuilder.build_cuttag_dirs(ref_genome)
            project_path = dirbuilder.project_path
            workingdir = project_path + 'logs' + '/' + dataset_name
            os.makedirs(workingdir, mode=0o777, exist_ok=True)
        if data_type == 'CR':
            dirbuilder.build_cutrun_dirs(ref_genome)
            project_path = dirbuilder.project_path
            workingdir = project_path + 'logs' + '/' + dataset_name
            os.makedirs(workingdir, mode=0o777, exist_ok=True)
        if data_type == 'DE':
            dirbuilder.build_deseq2_dirs(ref_genome)
            project_path = dirbuilder.project_path
            workingdir = project_path + 'logs' + '/' + dataset_name
            os.makedirs(workingdir, mode=0o777, exist_ok=True)
        if data_type == 'TO':
            dirbuilder.build_tornado_dirs(ref_genome)
            project_path = dirbuilder.project_path
            workingdir = project_path + 'logs' + '/' + dataset_name
            os.makedirs(workingdir, mode=0o777, exist_ok=True)

    else:
        dirbuilder.build_bclfastq_dirs()
        dataset_dir = dirbuilder.dataset_dir
        workingdir = dataset_dir + 'logs' + '/'
        os.makedirs(workingdir, mode=0o777, exist_ok=True)

    # Create directory for log files, if necessary
    #if not os.path.exists(log_dir):
    #    os.makedirs(log_dir)

    # Set custom names for output and error files for better readability
    output_log = dataset_name + "_snakemake.out"
    error_log = dataset_name + "_snakemake.err"

    ###
    # Build arguments to snakemake based on snakefile and specified parameters
    ###

    # If specified, perform a dry run or calculate the DAG only
    args = ["snakemake"]
    if dryrun:
        args.extend(["-n"])
    if dag:
        args.extend(["--dag"])

    # Run on specified snakefile (or default "Snakefile", if none specified)  
    args.extend(["--snakefile", snakefile])

    # Pass dataset to snakefile as config "dataset" parameter
    args.extend(["--config", "dataset=" + dataset_name])

    # If reference genome is specified, pass "ref_genome" as second parameter
    if ref_genome:
        args.extend(["ref_genome=" + ref_genome])

    # If ChIPseq bedtools_intersect arguments are specified, pass the parameters
    if enhancers:
        args.extend(["enhancers=" + enhancers])
    if gene_spans:
        args.extend(["gene_spans=" + gene_spans])
    if tss_spans:
        args.extend(["tss_spans=" + tss_spans])
    if blacklist:
        args.extend(["blacklist=" + blacklist])

    if ds_method:
        args.extend(["downsample_method=" + ds_method]) 

    if print_shell:
        args.extend(["-p"])

    # Cluster args -- exclude these if "localrun" is specified (with any value)
    #  (but set number of cores instead)
    if localrun:
        args.extend(["--cores", str(cores)])
    else:
        if snakefile == 'Snakefile_tornado':
            cluster_params = "qsub -V -wd " + workingdir
        else:
            cluster_params = "qsub -wd " + workingdir
        args.extend(["--cluster", cluster_params])
        args.extend(["-j", str(jobs)])

    if resources:
        args.extend(["--resources", "nodes="+str(nodes)])

    if lock:
        args.extend(["--nolock"])

    if keep_going:
        args.extend(["-k"])

    if ambiguity:
        args.extend(["--allow-ambiguity"])

    # Run the final snakemake command
    print('Running: ' + ' '.join(args))
    subprocess.Popen(args)

runsnakemake(args.dataset)
