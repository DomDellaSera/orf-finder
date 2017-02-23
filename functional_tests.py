import sys
import os
import unittest
from re_comp import *
from fasta_io import * #THIS IS BAD FORM


#bashCommand = "grep '>' test_fasta/brugia_malayi_reviewed.fasta | wc -l"
class bash_Command(object):
    def __init__(self, bash_Command):
        self.bash_Command = bash_Command
    def io(self):

        self.file_input = ""
        self.file_output = ""
        

class number_of_seq_io_equivalantTestCase(unittest.TestCase):
    def setUp(self):
        bash_Command = 'python fasta_io.py -i "test_fasta/brugia_malayi_reviewed.fasta" -o "test_fasta/test_output.fasta" --debug' 
        os.system(bash_Command)


    def tearDown(self):
        bash_Command = "echo '' > test_fasta/test_output.fasta"
        os.system(bash_Command)

    def greaterthan_number_eqiv(self):
        #Based on grep command, both sequences have similar lengths
        bash_Command = ["grep '>' test_fasta/brugia_malayi_reviewed.fasta | wc -l","grep '>' test_fasta/test_output.fasta | wc -l"]
        assertEqual(os.system(bash_Command[0]), os.system(bash_Command[1]))


#if __name__ == '__main__':
 #   unittest.main(warnings='ignore')
#os.system(bashCommand)

