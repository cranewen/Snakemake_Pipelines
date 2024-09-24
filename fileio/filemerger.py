from collections import defaultdict
import subprocess
import argparse
import os, errno
import logging
import glob
import sys
import atexit


class FileMerger():
    def __init__(self):
        self.input_path = ''
        self.merged_file_path = '' # create a new directory for merged files
        self.file_list = [] # A array stores all the file names from one directory
        self.file_groups = [] # An array stores all the groups of samples, e.g. S1, S2...
        self.merged_file_list = [] # An array stores merged files' name with full path 
        self.file_dict = defaultdict() # A dictionary stores all the samples

    # read file path, and get a self.file_groups
    def set_input_path(self, file_path):
        if file_path[-1] != 0:
            self.input_path = file_path + '/'

    def set_merged_file_path(self, file_path):
        self.merged_file_path = file_path
    
    # read input path and store all the gz files into a list
    def read_path(self):
        if os.path.exists(self.input_path):
            self.file_list = sorted(glob.glob(self.input_path + '*.fastq.gz'))
        else:
            log_filename = 'logs/bcl2fastq_logs/' + self.input_path.split('/')[-1] + '.log'
            logging.basicConfig(filename=log_filename, level=logging.DEBUG,\
                 format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
            logging.debug(str(self.input_path) + " can not be found.")
            logging.info(str(self.input_path) + " can not be found.")

    def create_path(self):
        if not os.path.exists(self.merged_file_path):
            os.makedirs(self.merged_file_path)

    def set_merged_file_name(self, file_list):
        merged_file_name = file_list[0].split('/')[-1][0:-21] + '_' + file_list[0].split('/')[-1].split('_')[-2]
        return(merged_file_name)

    # Loopping through a folder and creating a dictionary which has key as sample name, 
    # value as multiple parts of the gz files
    def get_file_dict(self):
        temp_key = ''
        temp_s_num = ''

        if len(self.file_list) != 0:
            for f in self.file_list:
                if temp_key == '':
                    temp_key = f.split('/')[-1][0:-21]
                    self.file_dict[temp_key] = []
                    self.file_dict[temp_key].append(f)
                elif temp_key != f.split('/')[-1][0:-21]:
                    temp_key = f.split('/')[-1][0:-21]
                    self.file_dict[temp_key] = []
                    self.file_dict[temp_key].append(f)
                else:
                    self.file_dict[temp_key].append(f)
    
    # This is a function associated with get_file_dict(), to get the merged file name list, 
    # iterating the dictionary and get the keys as the file names    
    def get_merged_file_name_list(self):
        merged_file_name_list = []
        for k, v in self.file_dict.items():
            if not 'Undetermined' in k:
                merged_file_name_list.append(k)
        return merged_file_name_list


    def file_merge(self):
        temp_list = []
        temp_list_r1 = []
        temp_list_r2 = []
        temp_s_num = '' # tracking the S# for grouping purpose

        if len(self.file_list) != 0:
            for f in self.file_list:
                if temp_s_num == '':
                    temp_s_num = f.split('/')[-1].split('_')[-4]
                    if f.split('/')[-1].split('_')[-2] == 'R1':
                        temp_list_r1.append(f)
                    if f.split('/')[-1].split('_')[-2] == 'R2':
                        temp_list_r2.append(f)
                if temp_s_num != f.split('/')[-1].split('_')[-4]:
                    self.file_groups.append(temp_list_r1)
                    self.file_groups.append(temp_list_r2)
                    temp_list_r1 = []
                    temp_list_r2 = []
                    if f.split('/')[-1].split('_')[-2] == 'R1':
                        temp_list_r1.append(f)
                    if f.split('/')[-1].split('_')[-2] == 'R2':
                        temp_list_r2.append(f)
                    temp_s_num = f.split('/')[-1].split('_')[-4]
                else:
					
                    if f.split('/')[-1].split('_')[-2] == 'R1':
                        temp_list_r1.append(f)
                        temp_s_num = f.split('/')[-1].split('_')[-4]
                    if f.split('/')[-1].split('_')[-2] == 'R2':
                        temp_list_r2.append(f)
                        temp_s_num = f.split('/')[-1].split('_')[-4]


            for f in self.file_groups:
                merged_file_name = self.merged_file_path + self.set_merged_file_name(f) + '_merged.fastq.gz'
                self.merged_file_list.append(merged_file_name)
                file_list = ' '.join(f)
                cat_cmd = "cat " + file_list + " > " + merged_file_name
                process = subprocess.Popen(cat_cmd, stdout=subprocess.PIPE, shell=True)
                process.communicate()


        else:
            log_filename = 'logs/bcl2fastq_logs/' + self.input_path.split('/')[-1] + '.log'
            logging.basicConfig(filename=log_filename , level=logging.DEBUG, \
                format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
            logging.debug("please check the bcl file folder.")
            logging.info("There is no file in the file_list, the bcl2fastq folder needs to be checked.")


def main():
    parser = argparse.ArgumentParser(description='Merge fastq gz files')
    parser.add_argument('--inpath', type=str)
    parser.add_argument('--outpath', type=str)
    args = parser.parse_args()
    fm = FileMerger()
    fm.set_input_path('path_for_input/')
    fm.set_merged_file_path('path_for_merged_file/')
    fm.create_path()
    fm.read_path()
    fm.file_merge()


if __name__ == "__main__":
    main()





    
