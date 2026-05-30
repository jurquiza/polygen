import unittest

from engine import PEG_PBS_LENGTHS, pegbldr

seq = 'atggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgacccacggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccaagctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaagtaa'


class TestPegbldr(unittest.TestCase):
    def assertPegResult(self, result, expected):
        self.assertEqual(result[1], expected[1])
        self.assertEqual(len(result[0]), len(expected[0]))
        for actual_part,expected_part in zip(result[0], expected[0]):
            self.assertEqual(actual_part[0], expected_part[0])
            self.assertEqual(actual_part[1], expected_part[1])
            self.assertEqual(actual_part[3], expected_part[3])
            if actual_part[1] == 'pegRNA':
                pbs_len = PEG_PBS_LENGTHS[actual_part[2].upper()]
                self.assertTrue(8 <= pbs_len <= 17)
                self.assertEqual(actual_part[2][:-pbs_len], expected_part[2][:-13])
            else:
                self.assertEqual(actual_part[2], expected_part[2])

    def test_peg_mut_single_forw(self):
        '''
        test if one pegRNA with a single mutation is designed with forward spacer
        '''
        edts = [['201', 'T', 'mut']]
        result = pegbldr(seq, edts, 'PE3')
        expected = ([['pegRNA0', 'pegRNA', 'CTCGTGACCACCCTGACCCAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCACTGCACGCAGTGGGTCAGGGTGGTCA', 'f'], ['gRNA0', 'gRNA', 'AGAAGTCGTGCTGCTTCATG', 'f']], None)
        self.assertPegResult(result, expected)
    
    def test_peg_mut_single_rev(self):
        '''
        test if one pegRNA with a single mutation is designed with reverse spacer
        '''
        edts = [['635', 'A', 'mut']]
        result = pegbldr(seq, edts, 'PE3')
        expected = ([['pegRNA0', 'pegRNA', 'CATGTGATCGCGCTTCTCGTGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCCAAAGACCCAAACGAGAAGCGCGATCA', 'r'], ['gRNA0', 'gRNA', 'CAGAACACCCCCATCGGCGA', 'r']], None)
        self.assertPegResult(result, expected)
    
    def test_peg_mut_multiple(self):
        '''
        test if one pegRNA with multiple mutations is designed reliably
        '''
        edts = [['280,283', 'C,T', 'mut']]
        result = pegbldr(seq, edts, 'PE3')
        expected = ([['pegRNA0', 'pegRNA', 'TGAAGAAGATGGTGCGCTCCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCAGGCTACGCCCTGGAGCGCACCATCTTC', 'r'], ['gRNA0', 'gRNA', 'TTCAAGTCCGCCATGCCCGA', 'r']], None)
        self.assertPegResult(result, expected)
        
    def test_peg_del_single_forw(self):
        '''
        test if one pegRNA with a single deletion is designed with forward spacer
        '''
        edts = [['342', '344', 'del']]
        result = pegbldr(seq, edts, 'PE3')
        expected = ([['pegRNA0', 'pegRNA', 'CAACTACAAGACCCGCGCCGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTGTCGCCCTCCTTCACCTCGGCGCGGGTCTTGTA', 'f'], ['gRNA0', 'gRNA', 'CGATGCCCTTCAGCTCGATG', 'f']], None)
        self.assertPegResult(result, expected)
        
    def test_peg_del_single_rev(self):
        '''
        test if one pegRNA with a single deletion is designed with reverse spacer
        '''
        edts = [['504', '508', 'del']]
        result = pegbldr(seq, edts, 'PE3')
        expected = ([['pegRNA0', 'pegRNA', 'CGCTGCCGTCCTCGATGTTGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCACTTCAAGATCCAACATCGAGGACGGC', 'r'], ['gRNA0', 'gRNA', 'CAGCCACAACGTCTATATCA', 'r']], None)
        self.assertPegResult(result, expected)
    
    def test_peg_del_multiple(self):
        '''
        test if one pegRNA with multiple deletions is designed reliably
        '''
        edts = [['425,430', '427,435', 'del']]
        result = pegbldr(seq, edts, 'PE3')
        expected = ([['pegRNA0', 'pegRNA', 'CAACATCCTGGGGCACAAGCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCACGTTGTGGCTGTTGTACAGCTTGTGCCCCAGGAT', 'f'], ['gRNA0', 'gRNA', 'GATGCCGTTCTTCTGCTTGT', 'f']], None)
        self.assertPegResult(result, expected)
    
    def test_peg_ins_single_forw(self):
        '''
        test if one pegRNA with a single insertion is designed with forward spacer
        '''
        edts = [['50', 'ACGT', 'ins']]
        result = pegbldr(seq, edts, 'PE3')
        expected = ([['pegRNA0', 'pegRNA', 'CACCGGGGTGGTGCCCATCCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCCAGCTCGACGTACCAGGATGGGCACCACCCC', 'f'], ['gRNA0', 'gRNA', 'CGCCGGACACGCTGAACTTG', 'f']], None)
        self.assertPegResult(result, expected)
        
    def test_peg_ins_single_rev(self):
        '''
        test if one pegRNA with a single insertion is designed with reverse spacer
        '''
        edts = [['603','CCC','ins']]
        result = pegbldr(seq, edts, 'PE3')
        expected = ([['pegRNA0', 'pegRNA', 'TCAGCTTGGACTGGGTGCTCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCAACCACTACCCCCTGAGCACCCAGTCCAAG', 'r'], ['gRNA0', 'gRNA', 'TACCAGCAGAACACCCCCAT', 'r']], None)
        self.assertPegResult(result, expected)
    
    def test_peg_ins_multiple(self):
        '''
        test if one pegRNA with multiple insertions is designed reliably
        '''
        edts = [['120,130', 'AAAA,GGGG', 'ins']]
        result = pegbldr(seq, edts, 'PE3')
        expected = ([['pegRNA0', 'pegRNA', 'GGTGGTGCAGATGAACTTCAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCCACCTACAAAAGGCAAGCTGAGGGGCCCTGAAGTTCATCTGCAC', 'r'], ['gRNA0', 'gRNA', 'AAGTTCAGCGTGTCCGGCGA', 'r']], None)
        self.assertPegResult(result, expected)
    
    def test_peg_multiple(self):
        '''
        test if multiple independent mutations are designed reliably on forward and reverse strand
        '''
        edts = [['201', 'T', 'mut'], ['504', '508', 'del'], ['120,130', 'AAAA,GGGG', 'ins']]
        result = pegbldr(seq, edts, 'PE3')
        expected = ([['pegRNA0', 'pegRNA', 'CTCGTGACCACCCTGACCCAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCACTGCACGCAGTGGGTCAGGGTGGTCA', 'f'], ['gRNA0', 'gRNA', 'AGAAGTCGTGCTGCTTCATG', 'f'], ['pegRNA1', 'pegRNA', 'CGCTGCCGTCCTCGATGTTGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCACTTCAAGATCCAACATCGAGGACGGC', 'r'], ['gRNA1', 'gRNA', 'CAGCCACAACGTCTATATCA', 'r'], ['pegRNA2', 'pegRNA', 'GGTGGTGCAGATGAACTTCAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCCACCTACAAAAGGCAAGCTGAGGGGCCCTGAAGTTCATCTGCAC', 'r'], ['gRNA2', 'gRNA', 'AAGTTCAGCGTGTCCGGCGA', 'r']], None)
        self.assertPegResult(result, expected)
