import unittest
from re_comp import translation, complement, reverse_complement,reading_frame_generator, reading_frame_translator,peptide_frame_chooser

class translationTestCase(unittest.TestCase):
	
	def test_AUG_start(self):
		self.assertTrue(translation("AUG") == "M")
			
	def test_Stop_codon(self):
		self.assertTrue(translation('UAG') == "*")
			
	def test_DNA_Start(self):
		self.assertTrue(translation('ATG') == "M")
			
	def test_UUA_and_CUU_Return_Leu(self):
		self.assertEqual(translation("UUA"), translation('CUU'))
			
			
class complementTestCase(unittest.TestCase):
	def test_A_T(self):
		self.assertEqual(complement("A"),"T")
	def test_G_C(self):
		self.assertTrue(complement("G") == "C")


class reverse_complementTestCase(unittest.TestCase):
	def test_ATGC(self):
		self.assertEqual(reverse_complement("ATGC"), "GCAT")
		
class reading_frame_generatorTestCase(unittest.TestCase):
	def test_ATGC(self):
		self.assertEqual(reading_frame_generator("ATGC"), ["ATGC", "TGC", "GC", "GCAT", "CAT", "AT"])

class reading_frame_translatorTestCase(unittest.TestCase):
	def test_ATGCGC(self):
		self.assertEqual(reading_frame_translator("ATGCGC"), "MR")

dna="TTATAGCGCATGTTAGCGCATGCATGCTATGCGCGATGTGTATAGTGACTGATATACAATGCCATGCATGCAAAATATAGCGCATGTTAGCGCATGCATGCTATGCGCGATGTGTATAGTGACTGATATACAATGCCATGCATGCAAAAATAGCGCATGTTAGCGCATGCATGCTATGCGCGATGTGTATAGTGACTGATATACAATGCCATGCATGCAAAATATAGCGCATGTTAGCGCATGCATGCTATGCGCGATGTGTATAGTGACTGATATACAATGCCATGCATGCAAAA"  

class peptide_frame_chooserTestCase(unittest.TestCase):

    def test_Sequence(self):
        self.assertEqual(peptide_frame_chooser("TTATAGCGCATGTTAGCGCATGCATGCTATGCGCGATGTGTATAGTGACTGATATACAATGCCATGCATGCAAAATATAGCGCATGTTAGCGCATGCATGCTATGCGCGATGTGTATAGTGACTGATATACAATGCCATGCATGCAAAAATAGCGCATGTTAGCGCATGCATGCTATGCGCGATGTGTATAGTGACTGATATACAATGCCATGCATGCAAAATATAGCGCATGTTAGCGCATGCATGCTATGCGCGATGTGTATAGTGACTGATATACAATGCCATGCATGCAAAA"  ),['HACANMRYILHAWHCISVTIHIAHSMHALTCAIFACMALY','HALFLHAWHCISVTIHIAHSMHALTCAIFCMHGIVYQSLY'])

        
        
        
if __name__ == '__main__':
	unittest.main()