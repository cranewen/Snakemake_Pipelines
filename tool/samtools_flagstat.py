from collections import defaultdict
import argparse
import sys, os 
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

parser = argparse.ArgumentParser(description = "samtools flagstat")
parser.add_argument('-fin', '--filein', type = str, required = True)
parser.add_argument('-fout', '--fileout', type = str, required = True)
parser.add_argument('-s', '--sample', nargs='+', required = True)
args = parser.parse_args()

# filein is a file path for each samtools_flagstat.txt. The format: {rmdup_bam_dataset_dir}/samtools_flagstat/sample_{ref_genome}_sorted_readgps_rmdup_samtools_flagstat.txt
# here sample is intented to be left as a string. So we won't get error in Snakemake due to not allowing to use single element from a list without modification
# So, use filein.split('sample') and re-assemble the whole path with sample argument
def parse_paired(filein, fileout, sample):
    for s in sample:
        filepath = filein.split('sample')[0] + s + filein.split('sample')[1]
        with open(filepath) as f:
            sample_stats_dict = defaultdict()
            sample_stats_dict[s] = []
            for line in f:
                if 'properly paired' in line or 'paired in sequencing' in line:
                    sample_stats_dict[s].append(line.split(' ')[0])
        with open(fileout, 'a') as fo:
            for k,v in sample_stats_dict.items():
                fo.write(k + '\t' + v[0] + '\t' + v[1] + '\n')
                
parse_paired(args.filein, args.fileout, args.sample)
