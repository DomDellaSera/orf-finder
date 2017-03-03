#!/usr/bin/python3
import logging
import sys
import getopt
import os
import traceback

import pdb

import textwrap as _textwrap
help_Message = """
Version: 0.1.1
Contact: Dominic Della Sera <dellasera@health.usf.edu>



This program attempts to find the most likely Open Reading frames for each DNA
sequence by finding the largest possible continuous sequences. ORF-Finder allows
for tolerance to missing nucleotides encoded as 'N', and assumes unidentifiable
codons are not stop codons.

"""

logger = logging.getLogger()
logging.basicConfig(level=logging.DEBUG)
logging.debug('logging loaded')




try:
    assert sys.version_info.major == 3
except:
    raise Exception('Python 3.x required')
    raise Exception('Python 3.x required')


def main():


    import re_comp
    import fasta_io
    logging.debug('main execution begins')
    inputfile = ''
    outputfile = 'output.fasta'
    #try:
    
    try:
        opts,args = getopt.getopt(sys.argv[1:], "i:ho")
        #if len(sys.argv)== 1:

           
    except:
        print('orf-finder.py -i <inputfile> -o <outputfile> [--debug] [--set_trace]')
        sys.exit()
        
    pdb.set_trace()
    
     
     
    stuff = []
    #if opts == []:
     #   opts = None
    for i in opts_args:
        if 'fasta' in i[1] & i[0] == []:
            opt = '-i'
        
        stuff.append((opt,arg))


    for opt,arg in stuff:
        print(opt,arg)
        if opt in ("-h", "--help"):
            print(help_Message)
            sys.exit()
        elif opt in ("-o", "--output"):
            outputfile = arg
        elif opt in ("-i", "--input"):
            inputfile = arg
        elif "fasta" in arg is True:
            inputfile = arg 
        else:
            sys.exit()
            
            
            
            
    # Function for command line arguements
    debug_mode = True
    pdb_debug = False



    #os.mkdir("output_data")
    # Ideally the output should be put in some directory/database
    



            
    if pdb_debug is True:
        pdb.set_trace()
    #if sys.version[0] == 2:
    #inputfile = input("Enter the Fasta file name:")
    #elif sys.version[0] == 3:
    #   inputfile = input("Enter the Fasta file name:")
    
        
    #inputfile = str(inputfile)
  
    logging.info(os.name)
    
    
    if os.name == "nt":
       
        #Windows doesn't do relative paths, so I'm going to grab the user directory from here
        script_dir = os.path.dirname(__file__)
        rel_path = inputfile
        inputfile = os.path.join(script_dir,rel_path)
        
            
    try:
        seq_objs=fasta_io.basic_fasta_parser(inputfile)
        #Now we have our sequences from this file but they haven't been converted    
    except:
        print(inputfile)
        tb = traceback.format_exc()
        print(tb)
        raise Exception("Basic Fasta Parser not working")
        
    #THis is outputting my class object
    
    #dna_reading_frames =
    
    for seq_obj in seq_objs:
        reading_frames=reading_frame_generator(seq_obj.sequence)
        try:
            reading_frames=reading_frame_generator(seq_obj.sequence)
            final_orfs = peptide_frame_chooser(reading_frames)
            seq_obj.sequence = final_orfs[0] 
        except:
            
            logger.debug('not found')


    formatted_fasta = [">"+x.header+"\n"+x.sequence+"\n\n" for x in seq_objs]
    with open(outputfile, "a") as out:
        for i in formatted_fasta:
            out.write(i)
    print("\n\noutput written to %s" % outputfile)
            
            
if __name__ == "__main__":
    main()
