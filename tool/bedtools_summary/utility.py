import yaml

class DoublelyLinkedList:
    def __init__(self):
        print('doublely linked list')

class LinkedList:
    def __init__(self):
        self.next = None
        self.previous = None
        print('linked list')

class YamlReader:
    def __init__(self, path):
        self.path = path
    # yaml.safe_load_all
    def read_yaml(self):
        try:
            return(yaml.safe_load_all(open(self.path)))
        except FileNotFoundError:
            print(str(self.path) + ' can not be found')
    # yaml.load_all
    def read_yaml_load_all(self):
        try:
            return(yaml.load_all(open(self.path), Loader=yaml.FullLoader))
        except FileNotFoundError:
            print(str(self.path) + ' can not be found')
        
class YamlDict:
    def __init__(self, path):
        self.path = path
    def read_yaml(self):
        try:
            with open(self.path) as f:
                return(yaml.full_load(f))
        except FileNotFoundError:
            print(str(self.path) + ' can not be found')
