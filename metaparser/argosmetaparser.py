import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
from metaparser.yamlhandle import YamlHandle
from metaparser.metaconf import ArgosTools 

class ArgosMetaParser:
    def __init__(self, argos_conf_path):
        self.argos_bin_path = ''
        self.argos_tool_path = ''
        self.argos_ngsplot_path = ''
        self.abs_tool_path = ''
    
    # Input argument tool as ArgosTools.BCL2FASTQ format
    def set_tool(self, tool):
        yaml_handle = YamlHandle(ArgosTools.ARGOS_CONFIG_PATH.value).read_yaml()
        for s in yaml_handle:
            self.argos_bin_path = s[ArgosTools.ARGOS_BIN_PATH.value]
            self.argos_tool_path = s[tool.value]
        
        self.abs_tool_path = self.argos_bin_path + self.argos_tool_path
        return(self.abs_tool_path)

    def set_tool_snakemake(self, tool):
        yaml_handle = YamlHandle(ArgosTools.ARGOS_CONFIG_PATH.value).read_yaml()
        for s in yaml_handle:
            self.argos_bin_path = s[ArgosTools.ARGOS_BIN_PATH_SNAKEMAKE.value]
            self.argos_tool_path = s[tool.value]
        
        self.abs_tool_path = self.argos_bin_path + self.argos_tool_path
        return(self.abs_tool_path)
    
    # ngsplot tool for tornado plot, it's separated from pipeline and snakemake envs
    def set_tool_ngsplot(self, tool):
        argos_ngsplot_path = ''
        argos_tool_path = ''
        yaml_handle = YamlHandle(ArgosTools.ARGOS_CONFIG_PATH.value).read_yaml()
        for s in yaml_handle:
            argos_ngsplot_path = s[ArgosTools.ARGOS_NGSPLOT_PATH.value] 
            argos_tool_path = s[tool.value]
        
        abs_tool_path = argos_ngsplot_path + argos_tool_path
        return(abs_tool_path)




        


def main():
    amp = ArgosMetaParser(ArgosTools.ARGOS_CONFIG_PATH.value)
    print(amp.set_tool(ArgosTools.PLOTHEATMAP))
    


if __name__ == "__main__":
    main()
