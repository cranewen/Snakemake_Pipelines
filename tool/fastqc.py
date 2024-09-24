import glob
import os
import sys
import subprocess
sys.path.append(os.path.abspath(os.path.join(
    os.path.dirname(__file__), os.path.pardir)))
from metaparser.argosmetaparser import ArgosMetaParser
from metaparser.metaconf import ArgosTools

class Fastqc:
    def __init__(self):
        self.input = ''
        self.output = ''
        self.argos_fastqc = ArgosMetaParser()
        self.argos_fastqc.set_tool(ArgosTools.FASTQC)

    def set_input(self, input_path):
        if os.path.exists(input_path):
            self.input = (', ').join(sorted(glob.glob(input_path + '*.fastq.gz')))

    def set_output(self, output_path):
        if not os.path.exists(output_path):
            os.mkdir(output_path)
            self.output = output_path
        self.output = output_path

    def run_fastqc(self):
        cat_cmd = self.argos_fastqc.abs_tool_path + \
            " -t 32 " + self.input + " -o " + self.output
        subprocess.run([self.argos_fastqc.abs_tool_path, "-t", "32", self.input, "-o", self.output])


def main():
    fq = Fastqc()
    fq.set_input(
        'test_input_path/')
    fq.set_output(
        'test_output_path/')
    print(fq.input)


if __name__ == "__main__":
    main()
