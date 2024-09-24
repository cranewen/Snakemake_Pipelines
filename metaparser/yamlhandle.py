import yaml
import os, errno

class YamlHandle:
    def __init__(self, yaml_path):
        self.yaml_path = yaml_path


    def read_yaml(self):
        try:
            return(yaml.safe_load_all(open(self.yaml_path)))
        except FileNotFoundError:
            print(str(self.yaml_path) + ' can not be found')