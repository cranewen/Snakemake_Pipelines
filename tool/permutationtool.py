import sys, os, re
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from metaparser.yamlhandle import YamlHandle
from metaparser.metaconf import SampleMeta
from metaparser.metaconf import ChipseqInputSampleMeta
from fileio.fileutil import FileUtil
from collections import defaultdict
from itertools import permutations, chain

class PermutationTool:
    def __init__(self, dataset_name, ref_genome):
        self.dataset_name = dataset_name
        self.ref_genome = ref_genome


    def permutation_samples(self, sample_list, num_of_files):
        return(permutations(sample_list, num_of_files))

    # permutation of 2 files in a list
    def perm2_samples(self, sample_list):
        perm_list = []
        len_of_list = len(sample_list)
        for i in range(len_of_list - 1):
            for j in range(i + 1, len_of_list):
               perm_list.append((sample_list[i], sample_list[j])) 
               perm_list.append((sample_list[j], sample_list[i])) 
        return(perm_list)

    # permutation of 2 lists to each other
    def perm2_sample_lists(self, a, b):
        perm_list = []
        for i in a:
            for j in b:
                perm_list.append((i, j))
                perm_list.append((j, i))
        return(perm_list)


    # permutation of 2 lists to each other include permutation of a itself
    def perm2_sample_lists_with_a(self, a, b):
        perm_list = []
        for i in a:
            for j in b:
                perm_list.append((i, j))
                perm_list.append((j, i))
        a_list = self.perm2_samples(a)
        [perm_list.append(a) for a in a_list]
        return(perm_list)

    # get a dictionary of {sample:perm_samples}
    def perm2_sample_dict(self, sample_list):
        perm_list = self.perm2_samples(sample_list)
        perm_dict = defaultdict()
        for s in sample_list:
            perm_dict[s] = []
        
        for p in perm_list:
            perm_dict[p[0]].append(p[1])
            
        return(perm_dict)

    # using gene_spans and TSS reference to compare other files 
    #print(enhancers)
    def ref_compare(self, ref_list, b):
        result_list = []
        for i in ref_list:
            for j in b:
                result_list.append((i, j))
        return(result_list)

        


def main():
    btp = PermutationTool('a', 'b')
    '''
    for i in btp.permutation_samples(['sample_s1', 'sample_s5','sample_s4','sample_s3','sample_s2'], 2):
        print(i)
    print(btp.perm2_samples(['sample_s1', 'sample_s5']))
    '''
    #print(btp.perm2_sample_lists_with_a(['sample_s1', 'sample_s5'] ,['sample_s4','sample_s3','sample_s2']))
    #print(btp.ref_compare(['sample_s1', 'sample_s5'] ,['sample_s4','sample_s3','sample_s2'], 2))
    #print(btp.perm2_samples(['sample_s1', 'sample_s5','sample_s4','sample_s3','sample_s2']))
    print(btp.ref_compare(['sample_s1', 'sample_s5','sample_s4','sample_s3','sample_s2'], ['enhancer1', 'enhancer2']))

    
    
    

if __name__ == "__main__":
    main()
        
