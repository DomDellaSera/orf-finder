

#Input/output functions
import os

from re_comp import *

#C:\Users\DDell\Documents\Programming\Python\ORF_Finder\test_fasta

logger = logging.getLogger()
def convert_bytes(num):
    """
    this function will convert bytes to MB.... GB... etc
    """
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if num < 1024.0:
            return("%3.1f %s" % (num, x))
        num /= 1024.0

def file_size(file_path):
    """
    this function will return the file size
    """
    if os.path.isfile(file_path):
        file_info = os.stat(file_path)
        
        if file_info.st_size > 300*(1024**2):
        
            raise Exception("File is too large for current implementation >(300MB)")
            
        return((convert_bytes(file_info.st_size)))
    
    else:
        pass
        #raise Exception("File not found")
        
class fasta_sequence(object): 
    """
    Formatting for Fasta Sequences
    """
    
    def __init__(self, header, sequence = ""):
        self.header = header.rstrip("\n")
        self.sequence = sequence
        
        
    def __repr__(self): # This forces the sequence to be printed, abstracting the attributes
        return('{}'.format(self.header)) 
    
def basic_fasta_parser(fasta_file):
    #Input: fasta
    #Output: A list of sequences
    
    

    logging.debug(fasta_file)
    logging.info(file_size(fasta_file))
    fasta_sequences = []
    
    with open(fasta_file) as f:
        fasta_lines = f.readlines()
        fasta_size = file_size(fasta_file)
        
        logging.info(
        """
    
                ##############
                # File Stats #
                ##############
                
                File:  %s
                Lines: %s
                Size:  %s
                """ % (fasta_file,
                len(fasta_lines),
                fasta_size))
        
        
        header = None # This will prevent appending 
        #
        #Fasta Generator - Should be a function
        # 
        for i in fasta_lines:
            if i.startswith(">"):
                if header is not None: 
                    logger.debug(fasta_sequence_tmp)
                    fasta_sequences.append(fasta_sequence_tmp)
                    
                header = i[1:]
                fasta_sequence_tmp = fasta_sequence(header)
                
                
                
            elif i.startswith("\n"):
                continue
            
            elif type(i) is str and len(i) > 0:
                stripped = i.rstrip("\n")
                fasta_sequence_tmp.sequence += stripped
        
      
    return([x for x in fasta_sequences])


    
                
                
            
                
        
    







    



