import unittest

from engine import (
    PTGbldr,
    oligo_fragment_mismatches,
    pegbldr,
    scarless_gg,
    tas_guide_design,
)


PE_TEST_SEQUENCE = (
    'atggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgacccacggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccaagctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaagtaa'
)


class TestOligoIntegrity(unittest.TestCase):
    def assert_no_oligo_mismatches(self, parts, poltype, enzm='bbsi'):
        polycistron = scarless_gg(parts, poltype=poltype, enzm=enzm)
        mismatches = oligo_fragment_mismatches(polycistron)
        self.assertEqual(mismatches, [])

    def test_ptg_cas9_oligos_match_fragments(self):
        parts = PTGbldr('ptg_check',
                        [['gRNA', 'ATAGGCTCCGATCGTGAACC'],
                         ['gRNA', 'GGAATTTGGCGAGGCCTATT'],
                         ['smRNA', 'CGTATGGCGTATCGGATTCTAT']],
                        poltype='ptg')
        self.assert_no_oligo_mismatches(parts, 'ptg')

    def test_prime_editing_oligos_match_fragments(self):
        inserts,msg = pegbldr(PE_TEST_SEQUENCE, [['201', 'T', 'mut']], 'PE3')
        self.assertIsNone(msg)
        parts = PTGbldr('prime_check',
                        [[insert[1], insert[2]] for insert in inserts],
                        poltype='ptg')
        self.assert_no_oligo_mismatches(parts, 'ptg')

    def test_cas12a_oligos_match_fragments(self):
        parts = PTGbldr('ca_check',
                        [['crRNA', 'ATAGGCTCCGATCGTGAACC'],
                         ['crRNA', 'GGAATTTGGCGAGGCCTATT'],
                         ['crRNA', 'CGTATGGCGTATCGGATTCT']],
                        poltype='ca')
        self.assert_no_oligo_mismatches(parts, 'ca')

    def test_tasa_oligos_match_fragments(self):
        targets = ['ATCATACTTAAAAATGAA',
                   'CGGATTATTAGGTGTCAT',
                   'GCGGAGTCGGAGCGCTGG']
        guides,msg = tas_guide_design('', system='TasA', exact_spacer='\n'.join(targets))
        self.assertIsNone(msg)
        parts = PTGbldr('tasa_check',
                        [['tigRNA', guide['tigRNA_dna']] for guide in guides],
                        poltype='tigRNA')
        self.assert_no_oligo_mismatches(parts, 'tigRNA')

    def test_tash_oligos_match_fragments(self):
        targets = ['ATCATACTTAAAATGA',
                   'CGGATTATTAGGTGTC',
                   'GCGGAGTCGGAGCGCT']
        guides,msg = tas_guide_design('', system='TasH', exact_spacer='\n'.join(targets))
        self.assertIsNone(msg)
        parts = PTGbldr('tash_check',
                        [['tigRNA', guide['tigRNA_dna']] for guide in guides],
                        poltype='tigRNA')
        self.assert_no_oligo_mismatches(parts, 'tigRNA')


if __name__ == '__main__':
    unittest.main()
