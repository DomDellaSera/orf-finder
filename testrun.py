from re_comp import *
from fasta_io import *
import logging
import os

inputfile = "test_fasta/random_seqs_dna.fasta"
outputfile = "testrun_output_ORFs.txt"

if os.name == "nt": #Windows doesn't do relative paths
    script_dir = os.path.dirname(__file__)
    rel_path = inputfile
    inputfile = os.path.join(script_dir,rel_path)

logger = logging.getLogger()
logger.setLevel(logging.INFO)
#logging.basicConfig(level=logging.INFO)



testdna = "TTATAGCGCATGTTAGCGCATGCATGCTATGCGCGATGTGTATAGTGACTGATATACAATGCCATGCATGCAAAATATAGCGCATGTTAGCGCATGCATGCTATGCGCGATGTGTATAGTGACTGATATACAATGCCATGCATGCAAAAATAGCGCATGTTAGCGCATGCATGCTATGCGCGATGTGTATAGTGACTGATATACAATGCCATGCATGCAAAATATAGCGCATGTTAGCGCATGCATGCTATGCGCGATGTGTATAGTGACTGATATACAATGCCATGCATGCAAAA"
td = reading_frame_generator(testdna)
print(peptide_frame_chooser(td))


seq_objs=basic_fasta_parser(inputfile) #Now we have our sequences from this file but they haven't been converted    
    #THis is outputting my class object
    
    #dna_reading_frames =
for seq_obj in seq_objs:
    reading_frames=reading_frame_generator(seq_obj.sequence)
    final_orfs = peptide_frame_chooser(reading_frames)
    seq_obj.sequence = final_orfs[0]
    
    logger.debug(final_orfs)


formatted_fasta = [">"+x.header+"\n"+x.sequence+"\n\n" for x in seq_objs]
with open(outputfile, "w") as out:
    for i in formatted_fasta:
        out.write(i)