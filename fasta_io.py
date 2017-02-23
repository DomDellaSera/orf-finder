#!/usr/bin/python3.5

#Input/output functions
import logging
import sys
import sys, getopt
import os




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
        raise Exception("File not found")
        
class fasta_sequence(object): 
    def __init__(self, header, sequence = ""):
        #self.sequence = ""
        #self.length = len(sequence)
        self.protein = False
        self.header = header
        self.sequence = sequence
        
        
    def __repr__(self): # This forces the sequence to be printed, abstracting the attributes
        return('{}'.format(self.header)) 
    
def basic_fasta_parser(fasta_file):
    
    

    logging.debug(fasta_file)
    print(file_size(fasta_file))
    fasta_sequences = []
    with open(fasta_file) as f:
        fasta_lines = f.readlines()
        
        fasta_size = file_size(fasta_file)
        
        
        logging.debug(
        """
    
                ##############
                # File Stats #
                ##############
                
                File:  %s
                Lines: %s
                Size:  %s
                """ % (fasta_file,
                len(fasta_lines),
                fasta_size)
                )
        #print(fasta_lines)
        seq_data = False
        
        header = None
        
        
        for i in fasta_lines:
            #logger.debug(i)
            
            if i.startswith(">"):
            
                if header is not None:#Append PREVIOUS fasta_sequence to list so we can reinitialize
                    logger.debug(fasta_sequence_tmp)
                    fasta_sequences.append(fasta_sequence_tmp)
                    
                header = i[1:]
                fasta_sequence_tmp = fasta_sequence(header)
                
                
                
            elif i.startswith("\n"):
                continue
            
            elif type(i) is str and len(i) > 0:
                stripped = i.rstrip("\n")
                fasta_sequence_tmp.sequence += stripped
        
        #fasta_sequences.append(fasta_sequence_tmp) #Add the last sequence
        #if fasta_sequences == 0:
        #    raise Exception("Something Went Wrong")
        logger.debug("Final Fasta Sequence")
        #logger.debug("Length of final fasta list" +len(fasta_sequences))
        #for i in fasta_sequences:
            #logger.debug(i)
            #logger.debug(i.sequence)
    logger.debug([x.sequence for x in fasta_sequences])
    
                
                
            
                
        
    








def main(argv):
    # Function for command line arguements
    debug_mode = False


    inputfile = ''
    outputfile = ''
    #try:
    opts, args = getopt.getopt(argv,"dhi:o:",["ifile=","ofile=","debug"])
    #except getopt.GetoptError:
     #   print('fasta_io.py -i <inputfile> -o <outputfile> [-debug]')
     #   sys.exit(2)
        
    for opt, arg in opts:
        if opt == '-h':
            print('test.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt == ("--debug"):
            logging.basicConfig(level=logging.DEBUG)
            
            debug_mode = True
            
            
            
    if debug_mode is False:
        fasta_file_relative_path = input("Enter the Fasta file name:")
        
            
         
    basic_fasta_parser(inputfile)
    #return()

if __name__ == "__main__":
   main(sys.argv[1:])
   
   

    



