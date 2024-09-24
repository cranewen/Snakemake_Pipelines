import os, sys, glob, re

class FileUtil:
    def __init__(self, filepath):
        self.filepath = filepath
 
    @staticmethod
    def get_fastq_sample_list(filepath, sample_list):
        file_list = ''
        fastq_sample_name_list = []

        if os.path.exists(filepath):

            # Look for _R1 files; should always exist for single- or paired-end
            file_list = glob.glob(filepath + '*_R1_*.fastq.gz') 

            # Compile all requested sample names into a single regex
            sample_regex = re.compile('|'.join(sample_list))

            # Cross-check all FASTQ (_R1) files against all sample names
            filtered_fastq_list = list(filter(sample_regex.search, file_list))

            # Strip off path and trailing info we don't need
            fastq_sample_name_list = [re.sub(r'_R1_.*','',fastq_file).split('/')[-1] for fastq_file in filtered_fastq_list]

        return fastq_sample_name_list

    @staticmethod
    def get_gz_file_list(filepath):
        if os.path.exists(filepath):
            return glob.glob(filepath + '*.fastq.gz') 

    @staticmethod
    # a fuction gets all the gz files with '.fastq.gz' stripped
    def get_gz_sample_name_list(filepath):
        file_list = ''
        if os.path.exists(filepath):
            file_list = glob.glob(filepath + '*.fastq.gz') 
        sample_name_list = [s.strip('.fastq.gz').split('/')[-1] for s in file_list]
        return sample_name_list 



def main():
    print('some tests here')

if __name__ == "__main__":
    main()
