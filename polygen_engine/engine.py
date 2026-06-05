## Import all necessary packages
import numpy as np
import csv
import itertools
import re
import pandas as pd
import json
import copy
import os
from Bio.SeqUtils import MeltingTemp as mt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


class InvalidUsage(Exception):
    '''
    Handler for user-interpretable errors
    
    :param message: A user-interpretable error message
    :type message: str
    :param status_code: The status code of the server, defaults to 400
    :type status_code: int, optional
    :param payload: Additional information on the error location, defaults to None
    :type payload: dict, optional
    '''

    def __init__(self, message, status_code=400, payload=None):
        '''inits InvalidUsage'''
        
        Exception.__init__(self)
        self.message = message
        if status_code is not None:
            self.status_code = status_code
        self.payload = payload

    def to_dict(self):
        '''
        Turns error object into dictionary
        
        :return: A dictionary containing the error message under key 'message' and the payload under the respective keys
        :rtype: dict
        '''
        
        rv = dict(self.payload or ())
        rv['message'] = self.message
        return rv

#Copyright (c) 2019 Scott Weisberg
class Part:
    '''
    Object used for defining sequence characteristics
    
    :param name: Sequence name
    :type name: str
    :param type: Sequence type
    :type type: type
    :param sequence: A nucleotide sequence consisting of [ACGTacgt]
    :type sequence: str
    '''
    
    def __init__(self,name,type,sequence):
        '''inits Part'''
        
        self.name = name
        self.sequence = sequence
        self.type = type
    
    def to_json(self):
        '''transforms data to json format'''
        
        return {
            "name": self.name,
            "sequence": self.sequence,
            "type": self.type,
            "primer_forward": self.primer_forward,
            "primer_reverse": self.primer_reverse,
            "primer_forward_tm": self.primer_forward_tm,
            "primer_reverse_tm": self.primer_reverse_tm,
            "localisation": self.localisation
        }
    
    @classmethod
    def from_json(cls, json_data):
        return cls(json_data['name'],
                   json_data['type'],
                   json_data['sequence']
                   )
        
    primer_forward = ""
    primer_reverse = ""
    primer_forward_tm = 0
    primer_reverse_tm = 0
    
    
class Polycistron:
    
    def __init__(self):
        '''inits Polycistron'''
    
    sequence = ''
    parts = []
    oligos = []
    features = []
    warning = None


def add_feature(polycistron, start, end, feature_type, label=None, note=None, strand=1):
    '''Adds a GenBank feature with readable label and note qualifiers.'''

    qualifiers = {}
    if label is not None:
        qualifiers['label'] = label
    if note is not None:
        qualifiers['note'] = note
    polycistron.features.append(SeqFeature(FeatureLocation(start, end, strand=strand),
                                           type=feature_type,
                                           qualifiers=qualifiers))


# Define functions to ease the main computation
def polyToJson(poly):
    '''transforms data to json format''' 
    
    p = copy.deepcopy(poly)
    
    ftrs = [vars(feat) for feat in p.features]
    for feat in ftrs:
        feat['location'] = vars(feat['location'])
    
    return {
        "sequence": p.sequence,
        "parts": [p.to_json() for p in p.parts],
        "oligos": p.oligos,
        "features": ftrs,
        "warning": p.warning,
        "overhang_selection_mode": getattr(p, 'overhang_selection_mode', 'optimal'),
        "overhang_warning": getattr(p, 'overhang_warning', ''),
        "selected_overhangs": getattr(p, 'selected_overhangs', []),
        "overhang_statuses": getattr(p, 'overhang_statuses', []),
        "fill_in_overlap_tms": getattr(p, 'fill_in_overlap_tms', [])
    }

def reverse_complement(sequence):
    '''
    Finds the reverse complement of a provided sequence
    
    :param sequence: The nucleic acid sequence consisting of [ACGTacgt]
    :type sequence: str
    
    :return: The reverse complement of the provided sequence
    :rtype: str
    '''
    
    rev_comp = ""
    Watson_Crick = {"A":"T","C":"G","T":"A","G":"C","a":"t","t":"a","c":"g","g":"c"}
    for base in sequence:
        rev_comp = Watson_Crick[base] + rev_comp
    return rev_comp

def complement(sequence):
    '''
    Finds the complement of a provided sequence
    
    :param sequence: A nucleic acid sequence consisting of [ACGTacgt]
    :type sequence: str
    
    :return: The complement of the provided sequence
    :rtype: str
    '''
    
    comp = ''
    Watson_Crick = {"A":"T","C":"G","T":"A","G":"C","a":"t","t":"a","c":"g","g":"c"}
    for base in sequence:
        comp += Watson_Crick[base]
    return comp

def dna_to_rna(sequence):
    '''Converts DNA bases to RNA bases while preserving orientation.'''

    return sequence.upper().replace('T', 'U')

def gc_content(sequence):
    '''Calculates GC content as a percentage.'''

    if len(sequence) == 0:
        return 0
    sequence = sequence.upper()
    return 100 * float(sequence.count('G') + sequence.count('C')) / len(sequence)

def longest_homopolymer(sequence):
    '''Finds the longest single-base run in a sequence.'''

    sequence = sequence.upper()
    if not sequence:
        return 0
    longest = current = 1
    for i in range(1, len(sequence)):
        if sequence[i] == sequence[i-1]:
            current += 1
        else:
            current = 1
        longest = max(longest, current)
    return longest

def reverse(sequence):
    '''
    Finds the reverse of a provided sequence
    
    :param sequence: A nucleic acid sequence consisting of [ACGTacgt]
    :type sequence: str
    
    :return: The reverse of the provided sequence
    :rtype: str
    '''
    
    rev = ''
    for base in sequence:
        rev = base + rev
    return rev

def flattn(list_2d):
    '''
    Flattens a 2D list of lists into one 1D list
    
    :param list_2d: 2D list
    :type list_2d: list
    
    :return: 1D list with all sublists of the input concatenated
    :rtype: list
    '''
    
    flattnd = []
    for subl in list_2d:
        if type(subl) is list:
            flattnd += [s for s in subl]
        else:
            flattnd.append(subl)
    return flattnd

def Diff(li1, li2):
    '''
    Finds all elements of two lists, which are unique to one of the two
    
    :param li1: First list
    :type li1: list
    :param li2: Second list
    :type li2: list
    
    :return: Differential list
    :rtype: list
    '''
    
    return (list(list(set(li1)-set(li2)) + list(set(li2)-set(li1))))
    
    
## Set template sequences tRNA and gRNA scaffold
tRNA = 'aacaaagcaccagtggtctagtggtagaatagtaccctgccacggtacagacccgggttcgattcccggctggtgca'
scaffld = 'gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgc'
DR = 'aatttctactgttgtagat'
OVERHANGSET_DIR = os.path.join(os.path.dirname(__file__), 'overhangsets')
PBS_TARGET_TM = 30
PBS_MIN_LEN = 8
PBS_MAX_LEN = 17
PEG_PBS_LENGTHS = {}
TIGRNA_MIN_SPLIT_LEN = 12
MIN_FILL_IN_OVERLAP_TM = 45
TAS_SYSTEMS = {
    'TasA': {
        'spacer_len': 9,
        'edge_5': 'AACCG',
        'loop': 'AGTAACCCC',
        'edge_3': 'AGTG',
        'spacer_a_transform': 'reverse_complement',
        'spacer_b_transform': 'provided',
        'note': 'TasA uses two 9 nt programmable spacers between AACCG and AGTG edge repeats.'
    },
    'TasH': {
        'spacer_len': 8,
        'edge_5': 'GAGCGAGTTA',
        'loop': 'AAAACAATCA',
        'edge_3': 'AAGCGAGCCA',
        'spacer_a_transform': 'provided',
        'spacer_b_transform': 'reverse_complement',
        'note': 'TasH uses two 8 nt programmable spacers in an asymmetric multiplexing scaffold; spacer B is reverse-complemented from the target window.'
    }
}

def make_tas_candidate(target, system, edge_5, loop, edge_3, strand='provided', start=0, end=None, index=1, source=None):
    '''
    Builds one mature tigRNA candidate from one target window.

    Target windows are split into two programmable spacers. TasA keeps the
    previous target-window convention where spacer A is reverse-complemented.
    TasH uses multiplexing inputs where spacer A is embedded as provided and
    spacer B is reverse-complemented from the second half of the target window.
    
    :param target: Exact target window sequence
    :type target: str
    :param system: Tas system name, currently 'TasA' or 'TasH'
    :type system: str
    :param edge_5: 5 prime tigRNA edge-repeat sequence
    :type edge_5: str
    :param loop: tigRNA loop-repeat sequence between spacer A and spacer B
    :type loop: str
    :param edge_3: 3 prime tigRNA edge-repeat sequence
    :type edge_3: str
    :returns: Candidate metadata dictionary containing spacers, mature tigRNA,
        GC content, score and PTG transfer string
    :rtype: dict
    '''

    config = TAS_SYSTEMS[system]
    spacer_len = config['spacer_len']
    target = target.upper()
    if end is None:
        end = len(target)
    spacer_a_input = target[:spacer_len]
    spacer_a_dna = reverse_complement(spacer_a_input) if config.get('spacer_a_transform') == 'reverse_complement' else spacer_a_input
    spacer_b_input = target[spacer_len:]
    spacer_b_dna = reverse_complement(spacer_b_input) if config.get('spacer_b_transform') == 'reverse_complement' else spacer_b_input
    tigRNA_dna = edge_5.replace('U', 'T') + spacer_a_dna + loop.replace('U', 'T') + spacer_b_dna + edge_3.replace('U', 'T')
    gc = gc_content(target)
    penalty = abs(50 - gc) + max(0, longest_homopolymer(tigRNA_dna) - 3) * 10
    return {
        'name': system + '_tigRNA_' + str(index),
        'system': system,
        'strand': strand,
        'start': start,
        'end': end,
        'source': source if source is not None else str(start) + '-' + str(end),
        'target': target,
        'spacer_a': dna_to_rna(spacer_a_dna),
        'spacer_b': dna_to_rna(spacer_b_dna),
        'tigRNA': dna_to_rna(tigRNA_dna),
        'tigRNA_dna': tigRNA_dna,
        'gc': round(gc, 1),
        'score': round(penalty, 2),
        'ptg_part': 'tigRNA;' + tigRNA_dna
    }

def annotate_tigRNA_features(polycistron, start, sequence, shared_edge_5_start=None, shared_edge_5=None, shared_edge_3_start=None, shared_edge_3=None):
    '''
    Annotates scaffold and spacer regions within one mature tigRNA.

    Recognized TasA/TasH tigRNAs receive separate GenBank ``misc_feature``
    entries for edge repeat 5 prime, spacer A, loop repeat, spacer B and edge
    repeat 3 prime. Unknown tigRNA scaffolds receive two broad fallback
    spacer/scaffold annotations.
    
    :param polycistron: Polycistron object receiving GenBank features
    :type polycistron: Polycistron object
    :param start: Start coordinate of the tigRNA within the assembled sequence
    :type start: int
    :param sequence: Mature tigRNA DNA sequence
    :type sequence: str
    :returns: Matched Tas system name, or None when the scaffold is unknown
    :rtype: str or None
    '''

    def append_tig_feature(rel_start, rel_end, label):
        add_feature(polycistron, start+rel_start, start+rel_end, 'misc_feature',
                    label=label, note=label)

    sequence = sequence.lower()
    shared_edge_5 = shared_edge_5.replace('U', 'T').lower() if shared_edge_5 else None
    shared_edge_3 = shared_edge_3.replace('U', 'T').lower() if shared_edge_3 else None
    for system, config in TAS_SYSTEMS.items():
        edge_5 = config['edge_5'].replace('U', 'T').lower()
        loop = config['loop'].replace('U', 'T').lower()
        edge_3 = config['edge_3'].replace('U', 'T').lower()
        spacer_len = config['spacer_len']
        expected_len = len(edge_5) + len(edge_3) + len(loop) + 2 * spacer_len
        if len(sequence) == expected_len and sequence.startswith(edge_5) and sequence.endswith(edge_3):
            spacer_a_start = len(edge_5)
            spacer_a_end = spacer_a_start + spacer_len
            loop_start = spacer_a_end
            loop_end = loop_start + len(loop)
            spacer_b_start = loop_end
            spacer_b_end = spacer_b_start + spacer_len
            if sequence[loop_start:loop_end] == loop:
                append_tig_feature(0, len(edge_5), system + ' tigRNA edge repeat 5 prime')
                append_tig_feature(spacer_a_start, spacer_a_end, system + ' tigRNA spacer A')
                append_tig_feature(loop_start, loop_end, system + ' tigRNA loop repeat')
                append_tig_feature(spacer_b_start, spacer_b_end, system + ' tigRNA spacer B')
                append_tig_feature(spacer_b_end, len(sequence), system + ' tigRNA edge repeat 3 prime')
                if system == 'TasH' and edge_5 == edge_3 and len(edge_5) % 2 == 0:
                    cut_offset = int(len(edge_5) / 2)
                    processed_start = cut_offset
                    processed_end = len(sequence) - cut_offset
                    processed_seq = sequence[processed_start:processed_end].upper()
                    append_tig_feature(0, cut_offset, system + ' discarded 5 prime ER half')
                    append_tig_feature(cut_offset, len(edge_5), system + ' retained 5 prime ER half')
                    append_tig_feature(spacer_b_end, spacer_b_end + cut_offset, system + ' retained 3 prime ER half')
                    append_tig_feature(spacer_b_end + cut_offset, len(sequence), system + ' discarded 3 prime ER half')
                    add_feature(polycistron, start+processed_start, start+processed_end, 'misc_feature',
                                label=system + ' processed mature tigRNA',
                                note=system + ' mature tigRNA after ER processing; sequence=' + processed_seq)
                    add_feature(polycistron, start+cut_offset-1, start+cut_offset+1, 'misc_feature',
                                label=system + ' 5 prime ER processing cut',
                                note='Processing cut occurs between the two halves of the 5 prime edge repeat')
                    add_feature(polycistron, start+processed_end-1, start+processed_end+1, 'misc_feature',
                                label=system + ' 3 prime ER processing cut',
                                note='Processing cut occurs between the two halves of the 3 prime edge repeat')
                return system
        shared_expected_len = 2 * spacer_len + len(loop) + len(edge_3)
        if (shared_edge_5_start is not None and shared_edge_5 == edge_5 and
                len(sequence) == shared_expected_len and sequence.endswith(edge_3)):
            spacer_a_start = 0
            spacer_a_end = spacer_len
            loop_start = spacer_a_end
            loop_end = loop_start + len(loop)
            spacer_b_start = loop_end
            spacer_b_end = spacer_b_start + spacer_len
            if sequence[loop_start:loop_end] == loop:
                append_tig_feature(spacer_a_start, spacer_a_end, system + ' tigRNA spacer A')
                append_tig_feature(loop_start, loop_end, system + ' tigRNA loop repeat')
                append_tig_feature(spacer_b_start, spacer_b_end, system + ' tigRNA spacer B')
                append_tig_feature(spacer_b_end, len(sequence), system + ' tigRNA edge repeat 3 prime')
                return system

    midpoint = int(np.floor(len(sequence) / 2))
    append_tig_feature(0, midpoint, 'tigRNA spacer/scaffold region 1')
    append_tig_feature(midpoint, len(sequence), 'tigRNA spacer/scaffold region 2')
    return None

def tigRNA_spacer_ranges(sequence):
    '''
    Returns programmable spacer coordinate ranges for a mature tigRNA.

    For recognized TasA/TasH scaffolds, the returned list contains spacer A and
    spacer B ranges. Golden Gate optimization can use either programmable
    spacer as the source for internal tigRNA overhangs. For unknown scaffolds,
    the whole sequence is returned as one fallback range.
    
    :param sequence: Mature tigRNA DNA sequence
    :type sequence: str
    :returns: List of ``(start, end)`` coordinate tuples
    :rtype: list
    '''

    sequence = sequence.lower()
    for system, config in TAS_SYSTEMS.items():
        edge_5 = config['edge_5'].replace('U', 'T').lower()
        loop = config['loop'].replace('U', 'T').lower()
        edge_3 = config['edge_3'].replace('U', 'T').lower()
        spacer_len = config['spacer_len']
        spacer_a_start = len(edge_5)
        spacer_a_end = spacer_a_start + spacer_len
        loop_start = spacer_a_end
        loop_end = loop_start + len(loop)
        spacer_b_start = loop_end
        spacer_b_end = spacer_b_start + spacer_len
        expected_len = len(edge_5) + len(edge_3) + len(loop) + 2 * spacer_len
        if (len(sequence) == expected_len and sequence.startswith(edge_5) and
                sequence[loop_start:loop_end] == loop and sequence.endswith(edge_3)):
            return [(spacer_a_start, spacer_a_end), (spacer_b_start, spacer_b_end)]
        shared_expected_len = 2 * spacer_len + len(loop) + len(edge_3)
        shared_loop_start = spacer_len
        shared_loop_end = shared_loop_start + len(loop)
        shared_spacer_b_start = shared_loop_end
        shared_spacer_b_end = shared_spacer_b_start + spacer_len
        if (len(sequence) == shared_expected_len and
                sequence[shared_loop_start:shared_loop_end] == loop and sequence.endswith(edge_3)):
            return [(0, spacer_len), (shared_spacer_b_start, shared_spacer_b_end)]
    return [(0, len(sequence))]

def all_curated_overhangs():
    '''Returns all unique 4 bp overhangs present in the curated overhang tables.'''

    overhangs = set()
    for filename in os.listdir(OVERHANGSET_DIR):
        if not filename.startswith('setsof') or not filename.endswith('.csv'):
            continue
        with open(os.path.join(OVERHANGSET_DIR, filename), 'r') as f:
            reader = csv.reader(f, delimiter=',')
            next(reader, None)
            for row in reader:
                for overhang in row:
                    overhang = overhang.strip().lower()
                    if len(overhang) == 4 and re.search(r'^[acgt]*$', overhang) is not None:
                        overhangs.add(overhang)
    return sorted(overhangs)

def primer_annealing_tm(primer, max_ann_len=30):
    '''
    Calculates Tm for the 3 prime annealing portion of an oligo.

    This helper is available for template-annealing designs. tigRNA mode uses
    :func:`primer_pair_overlap_tm` instead, because tigRNA fragments are formed
    from overlapping oligo pairs.
    '''

    annealing = primer[-min(max_ann_len, len(primer)):]
    return mt.Tm_NN(annealing, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)

def pbs_tm(sequence):
    '''
    Calculates PBS melting temperature using the nearest-neighbor model used
    elsewhere in PolyGEN.
    '''

    return mt.Tm_NN(sequence, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)

def choose_pbs(sequence, pegPAM, target_tm=PBS_TARGET_TM, min_len=PBS_MIN_LEN, max_len=PBS_MAX_LEN):
    '''
    Chooses a prime-editing PBS with Tm closest to the target temperature.

    Candidate PBS lengths default to 8-17 nt. The default target Tm is 30 C,
    matching the experimentally favored PBS Tm region.
    '''

    options = []
    nick = pegPAM - 3
    for pbs_len in range(min_len, max_len + 1):
        if nick - pbs_len < 0:
            continue
        pbs = reverse_complement(sequence[nick - pbs_len:nick])
        tm = pbs_tm(pbs)
        options.append((abs(tm - target_tm), pbs_len, tm, pbs))

    if options == []:
        raise InvalidUsage("The provided sequence does not cover enough area around the edit", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'sequence'})

    options = sorted(options, key=lambda x: (x[0], x[1]))
    return options[0][3], options[0][1], options[0][2]

def pegRNA_pbs_length(pegRNA_sequence):
    '''
    Infers PBS length from the pegRNA suffix when explicit metadata is missing.
    '''

    pegRNA_sequence = pegRNA_sequence.upper()
    if pegRNA_sequence in PEG_PBS_LENGTHS:
        return PEG_PBS_LENGTHS[pegRNA_sequence]

    options = []
    for pbs_len in range(PBS_MIN_LEN, PBS_MAX_LEN + 1):
        if len(pegRNA_sequence) < pbs_len:
            continue
        pbs = pegRNA_sequence[-pbs_len:]
        tm = pbs_tm(pbs)
        options.append((abs(tm - PBS_TARGET_TM), pbs_len))

    if options == []:
        return 13
    return sorted(options, key=lambda x: (x[0], x[1]))[0][1]

def longest_common_substring(seq_a, seq_b):
    '''
    Finds the longest common substring between two sequences.

    :param seq_a: First sequence
    :type seq_a: str
    :param seq_b: Second sequence
    :type seq_b: str
    :returns: Longest exact substring shared by both inputs
    :rtype: str
    '''

    best = ''
    for a in range(len(seq_a)):
        for b in range(len(seq_b)):
            n = 0
            while a+n < len(seq_a) and b+n < len(seq_b) and seq_a[a+n] == seq_b[b+n]:
                n += 1
            if n > len(best):
                best = seq_a[a:a+n]
    return best

def primer_pair_overlap(forward_primer, reverse_primer):
    '''
    Finds the polymerase-fill overlap shared by a forward/reverse oligo pair.

    The reverse primer is reverse-complemented before comparison, so the result
    is the sequence that can anneal between both oligos and be filled in by
    polymerase.
    
    :param forward_primer: Forward oligo sequence
    :type forward_primer: str
    :param reverse_primer: Reverse oligo sequence
    :type reverse_primer: str
    :returns: Shared fill-in overlap sequence
    :rtype: str
    '''

    return longest_common_substring(forward_primer.lower(), reverse_complement(reverse_primer).lower())

def primer_pair_overlap_tm(forward_primer, reverse_primer):
    '''
    Calculates Tm for the shared overlap of a forward/reverse oligo pair.

    This is the Tm reported for tigRNA oligo-extension fragments. It uses only
    the overlap returned by :func:`primer_pair_overlap`, not the full forward or
    reverse oligo sequence.
    
    :param forward_primer: Forward oligo sequence
    :type forward_primer: str
    :param reverse_primer: Reverse oligo sequence
    :type reverse_primer: str
    :returns: Nearest-neighbor Tm of the shared overlap, or 0 if no overlap is found
    :rtype: float
    '''

    overlap = primer_pair_overlap(forward_primer, reverse_primer)
    if overlap == '':
        return 0
    return mt.Tm_NN(overlap, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)


def sequence_overlap_tm(seq_a, seq_b):
    '''Calculates Tm for the longest shared sequence between two same-strand oligo cores.'''

    overlap = longest_common_substring(seq_a.lower(), seq_b.lower())
    if overlap == '':
        return 0
    return mt.Tm_NN(overlap, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)


def overhang_combination_overlap_tms(parts_list, comb, poltype_opt='ptg', max_ann_len=30, enzm='bsai', bb_linkers=['tgcc','gttt']):
    '''
    Estimates fill-in overlap Tm values produced by an internal overhang choice.

    CA and tigRNA oligo-extension fragments are short enough that the exact split
    position strongly affects the overlap shared by adjacent oligos. This helper
    scores a candidate combination before primers are committed.
    '''

    if comb is None:
        return []

    tms = []
    enzms={'bsai': ['gaggtctcg', 'cgagacctc'], 'bsmbi': ['tgcgtctca', 'tgagacgca'], 'btgzi': ['ctgcgatggagtatgtta', 'taacatactccatcgcag'], 'bbsi': ['ttgaagactt', 'aagtcttcaa']}
    if poltype_opt == 'ca':
        current_seq = enzms[enzm][0] + reverse_complement(bb_linkers[0]) + parts_list[0].sequence.lower()
        for x in range(len(parts_list)-1):
            cut = current_seq.find(comb[x]) + len(comb[x])
            if cut < len(comb[x]):
                continue
            reverse_primer = reverse_complement(current_seq[current_seq.find(DR):cut] + enzms[enzm][1])
            next_forward_primer = enzms[enzm][0] + current_seq[cut-4:] + DR
            tms.append(primer_pair_overlap_tm(next_forward_primer, reverse_primer))
            current_seq = enzms[enzm][0] + current_seq[cut-4:] + parts_list[x+1].sequence.lower()
    elif poltype_opt == 'tigRNA':
        current_seq = enzms[enzm][0] + reverse_complement(bb_linkers[0]) + parts_list[0].sequence.lower()
        forward_primers = [None for i in parts_list]
        reverse_primers = [None for i in parts_list]
        forward_primers[0] = enzms[enzm][0] + reverse_complement(bb_linkers[0]) + parts_list[0].sequence.lower()[:max_ann_len]
        for x in range(len(parts_list)-1):
            cut = current_seq.find(comb[x]) + len(comb[x])
            if cut < len(comb[x]):
                continue
            current_tail = current_seq[cut-4:]
            if (getattr(parts_list[x+1], 'full_tigRNA_unit', None) and
                    getattr(parts_list[x+1], 'shared_edge_5', None) and
                    parts_list[x+1].shared_edge_5.lower() not in current_tail):
                return []
            if (getattr(parts_list[x+1], 'shared_edge_junction', None) and
                    parts_list[x+1].shared_edge_junction.lower() not in current_tail + parts_list[x+1].sequence.lower()[:len(parts_list[x+1].shared_edge_5)]):
                return []
            current_start = current_seq.find(parts_list[x].template_sequence)
            reverse_primers[x] = reverse_complement(current_seq[current_start:cut] + enzms[enzm][1])
            forward_primers[x+1] = enzms[enzm][0] + current_tail + parts_list[x+1].sequence.lower()[:max_ann_len]
            current_seq = enzms[enzm][0] + current_tail + parts_list[x+1].sequence.lower()
        reverse_primers[-1] = reverse_complement(current_seq[-max_ann_len:] + bb_linkers[-1] + enzms[enzm][1])
        for forward_primer,reverse_primer in zip(forward_primers, reverse_primers):
            if forward_primer is None or reverse_primer is None:
                return []
            tms.append(primer_pair_overlap_tm(forward_primer, reverse_primer))
    return tms


def overhang_combination_score(parts_list, comb, poltype_opt='ptg', enzm='bsai', bb_linkers=['tgcc','gttt']):
    '''Returns a sortable score for a candidate overhang combination.'''

    tms = overhang_combination_overlap_tms(parts_list, comb, poltype_opt, enzm=enzm, bb_linkers=bb_linkers)
    if tms == []:
        return (-1, 0, 0, 0, 0)
    min_tm = min(tms)
    avg_tm = sum(tms) / len(tms)
    spread = max(tms) - min_tm
    passes_floor = 1 if min_tm >= MIN_FILL_IN_OVERLAP_TM else 0
    return (passes_floor, min_tm, avg_tm, -spread, -sum([len(c) for c in comb]))


def candidate_overhang_prefixes(part, poltype_opt, allowed_overhangs):
    '''Returns possible overhang prefixes for one part.'''

    seq = part.sequence.lower()
    candidates = []
    if poltype_opt == 'tigRNA':
        for r_start,r_end in tigRNA_spacer_ranges(seq):
            for pos in range(r_start, r_end-3):
                overhang = seq[pos:pos+4]
                if overhang in allowed_overhangs and pos + 4 >= TIGRNA_MIN_SPLIT_LEN:
                    candidates.append(seq[:pos+4])
    elif poltype_opt == 'ca':
        spacer = seq[seq.find(DR)+len(DR):]
        offset = seq.find(DR)+len(DR)
        for pos in range(0, len(spacer)-3):
            overhang = spacer[pos:pos+4]
            if overhang in allowed_overhangs:
                candidates.append(seq[:offset+pos+4])
    return candidates


def improve_overhang_combination_tms(parts_list, comb, poltype_opt='ptg', enzm='bsai', bb_linkers=['tgcc','gttt'], allowed_overhangs=None):
    '''
    Locally improves a valid overhang combination by swapping weak junctions.

    This keeps runtime bounded for large Tas arrays while still allowing later
    spacer A/B choices to rescue low-Tm oligo-extension fragments.
    '''

    if poltype_opt not in ['ca', 'tigRNA'] or comb is None:
        return comb, []
    if allowed_overhangs is None:
        allowed_overhangs = set(all_curated_overhangs())
        allowed_overhangs.update([reverse_complement(o) for o in list(allowed_overhangs)])

    best_comb = tuple(comb)
    best_score = overhang_combination_score(parts_list, best_comb, poltype_opt, enzm=enzm, bb_linkers=bb_linkers)
    best_tms = overhang_combination_overlap_tms(parts_list, best_comb, poltype_opt, enzm=enzm, bb_linkers=bb_linkers)
    weak_indices = set()
    for idx,tm in enumerate(best_tms):
        if tm < MIN_FILL_IN_OVERLAP_TM:
            if idx > 0:
                weak_indices.add(idx - 1)
            if idx < len(best_comb):
                weak_indices.add(idx)
    if not weak_indices:
        return best_comb, []
    changed_indices = set()

    improved = True
    while improved and min(best_tms) < MIN_FILL_IN_OVERLAP_TM:
        improved = False
        current_overhangs = [c[-4:] for c in best_comb]
        for idx in sorted(weak_indices):
            for candidate in candidate_overhang_prefixes(parts_list[idx], poltype_opt, allowed_overhangs):
                candidate_overhang = candidate[-4:]
                if candidate_overhang in current_overhangs and candidate_overhang != current_overhangs[idx]:
                    continue
                trial = list(best_comb)
                trial[idx] = candidate
                trial = tuple(trial)
                if len([c[-4:] for c in trial]) != len(set([c[-4:] for c in trial])):
                    continue
                score = overhang_combination_score(parts_list, trial, poltype_opt, enzm=enzm, bb_linkers=bb_linkers)
                trial_tms = overhang_combination_overlap_tms(parts_list, trial, poltype_opt, enzm=enzm, bb_linkers=bb_linkers)
                if score[0] >= best_score[0] and trial_tms != [] and min(trial_tms) > min(best_tms):
                    if trial[idx] != best_comb[idx]:
                        changed_indices.add(idx)
                    best_comb = trial
                    best_score = score
                    best_tms = trial_tms
                    improved = True
                    break
            if improved:
                break

    return best_comb, sorted(changed_indices)


def tigRNA_seq_matches_for_overhangs(parts_list, golden_gate_overhangs, cov):
    '''Returns possible tigRNA overhang prefixes for each internal junction.'''

    oh_list = [x.sequence.lower() for x in parts_list]
    tigRNA_ranges = [tigRNA_spacer_ranges(x.sequence) for x in parts_list]
    seq_matches = []
    for x in range(len(parts_list)-1):
        matches = []
        for overhang in golden_gate_overhangs:
            for r_start,r_end in tigRNA_ranges[x]:
                spacer_seq = oh_list[x][r_start:r_end]
                if 2*cov >= len(spacer_seq):
                    spacer_cov = spacer_seq
                else:
                    spacer_mid = int(np.ceil(np.true_divide(len(spacer_seq),2)))
                    spacer_cov = spacer_seq[spacer_mid-cov:spacer_mid+cov]
                if overhang in spacer_cov:
                    overhang_start = oh_list[x].find(overhang, r_start, r_end)
                    if overhang_start + 4 >= TIGRNA_MIN_SPLIT_LEN:
                        matches.append(oh_list[x][:overhang_start+4])
        seq_matches.append(sorted(set(matches), key=lambda item: (len(item), item)))
    return seq_matches


def ptg_subset_rescue_allowed(parts_list):
    '''Returns whether PTG subset rescue can be used without prime editing.'''

    return all(part.type in ['tRNA', 'gRNA'] for part in parts_list)


def ptg_seq_matches_for_overhangs(parts_list, golden_gate_overhangs, cov):
    '''Returns possible Cas9/gRNA PTG overhang prefixes for each junction.'''

    oh_list = []
    for part in parts_list:
        if part.type == 'gRNA':
            oh_list.append(part.sequence[:part.sequence.find(scaffld)])
        else:
            oh_list.append('plc')

    cov_list = []
    for oh in oh_list:
        if 2*cov >= len(oh):
            cov_list.append(oh)
        else:
            midpoint = int(np.ceil(np.true_divide(len(oh),2)))
            cov_list.append(oh[midpoint-cov:midpoint+cov])

    seq_matches = []
    for x in range(len(parts_list)-1):
        matches = []
        if parts_list[x+1].type == 'gRNA':
            for overhang in golden_gate_overhangs:
                if overhang in cov_list[x+1]:
                    matches.append(oh_list[x+1][:oh_list[x+1].find(overhang)+4])
        seq_matches.append(sorted(set(matches), key=lambda item: (len(item), item)))
    return seq_matches


def ca_seq_matches_for_overhangs(parts_list, golden_gate_overhangs, cov):
    '''Returns possible Cas12a/crRNA overhang prefixes for each junction.'''

    oh_list = [part.sequence[part.sequence.find(DR)+len(DR):] for part in parts_list]
    cov_list = []
    for oh in oh_list:
        if 2*cov >= len(oh):
            cov_list.append(oh)
        else:
            midpoint = int(np.ceil(np.true_divide(len(oh),2)))
            cov_list.append(oh[midpoint-cov:midpoint+cov])

    seq_matches = []
    for x in range(len(parts_list)-1):
        matches = []
        for overhang in golden_gate_overhangs:
            if overhang in cov_list[x]:
                matches.append(oh_list[x][:oh_list[x].find(overhang)+4])
        seq_matches.append(sorted(set(matches), key=lambda item: (len(item), item)))
    return seq_matches


def junction_candidate_prefixes(parts_list, idx, poltype_opt, allowed_overhangs):
    '''Returns rescue overhang prefixes for one internal junction.'''

    if poltype_opt == 'tigRNA':
        return candidate_overhang_prefixes(parts_list[idx], 'tigRNA', allowed_overhangs)
    if poltype_opt == 'ca':
        return candidate_overhang_prefixes(parts_list[idx], 'ca', allowed_overhangs)
    if poltype_opt == 'ptg':
        if not ptg_subset_rescue_allowed(parts_list) or parts_list[idx+1].type != 'gRNA':
            return []
        spacer = parts_list[idx+1].sequence[:parts_list[idx+1].sequence.find(scaffld)].lower()
        candidates = []
        for pos in range(0, len(spacer)-3):
            overhang = spacer[pos:pos+4]
            if overhang in allowed_overhangs:
                candidates.append(spacer[:pos+4])
        return candidates
    return []


def rescue_dropped_overhangs(parts_list, partial_comb, dropped_indices, used_overhangs, allowed_overhangs, poltype_opt, enzm='bsai', bb_linkers=['tgcc','gttt']):
    '''Fills dropped junctions with curated rescue overhangs.'''

    rescue_options = []
    for idx in dropped_indices:
        options = []
        for candidate in junction_candidate_prefixes(parts_list, idx, poltype_opt, allowed_overhangs):
            candidate_overhang = candidate[-4:]
            if candidate_overhang not in used_overhangs:
                options.append(candidate)
        options = sorted(set(options), key=lambda item: (len(item), item))
        if options == []:
            return None, None
        rescue_options.append(options)

    best_comb = None
    best_score = None
    for rescue_tuple in itertools.product(*rescue_options):
        trial = list(partial_comb)
        for idx,candidate in zip(dropped_indices, rescue_tuple):
            trial[idx] = candidate
        trial = tuple(trial)
        overhangs = [c[-4:] for c in trial]
        if len(overhangs) != len(set(overhangs)):
            continue
        if poltype_opt in ['ca', 'tigRNA']:
            score = overhang_combination_score(parts_list, trial, poltype_opt, enzm=enzm, bb_linkers=bb_linkers)
            if score[0] < 0:
                continue
            if best_score is None or score > best_score:
                best_comb = trial
                best_score = score
            if score[0] == 1:
                return trial, score
        else:
            return trial, (1, 0, 0, 0, 0)
    return best_comb, best_score


def subset_rescue_optimization(parts_list, free_overhangsets, rescue_overhangs, poltype_opt, enzm='bsai', bb_linkers=['tgcc','gttt']):
    '''
    Finds a mostly optimal overhang set by dropping difficult junctions.

    The kept junctions must come from one validated optimal overhang collection.
    Dropped junctions are filled from the curated rescue pool and can be flagged
    individually in the output.
    '''

    n_junctions = len(parts_list) - 1
    if n_junctions <= 1:
        return None, []
    if poltype_opt == 'ptg' and not ptg_subset_rescue_allowed(parts_list):
        return None, []
    max_cov = max([int(np.ceil(np.true_divide(len(prt.sequence),2))) for prt in parts_list])
    fallback_comb = None
    fallback_score = None
    fallback_dropped = []
    max_partial_trials = 2000

    for drop_count in range(1, n_junctions + 1):
        for dropped_indices in itertools.combinations(range(n_junctions), drop_count):
            dropped_indices = tuple(dropped_indices)
            kept_indices = [idx for idx in range(n_junctions) if idx not in dropped_indices]
            for cov in range(2, max_cov):
                for golden_gate_overhangs in free_overhangsets:
                    if poltype_opt == 'tigRNA':
                        seq_matches = tigRNA_seq_matches_for_overhangs(parts_list, golden_gate_overhangs, cov)
                    elif poltype_opt == 'ca':
                        seq_matches = ca_seq_matches_for_overhangs(parts_list, golden_gate_overhangs, cov)
                    elif poltype_opt == 'ptg':
                        seq_matches = ptg_seq_matches_for_overhangs(parts_list, golden_gate_overhangs, cov)
                    else:
                        return None, []
                    if any(seq_matches[idx] == [] for idx in kept_indices):
                        continue
                    partial_trials = 0
                    for kept_tuple in itertools.product(*[seq_matches[idx] for idx in kept_indices]):
                        partial_trials += 1
                        if partial_trials > max_partial_trials:
                            break
                        kept_overhangs = [candidate[-4:] for candidate in kept_tuple]
                        if len(kept_overhangs) != len(set(kept_overhangs)):
                            continue
                        partial_comb = [None for i in range(n_junctions)]
                        for idx,candidate in zip(kept_indices, kept_tuple):
                            partial_comb[idx] = candidate
                        full_comb,score = rescue_dropped_overhangs(parts_list,
                                                                   partial_comb,
                                                                   dropped_indices,
                                                                   set(kept_overhangs),
                                                                   rescue_overhangs,
                                                                   poltype_opt,
                                                                   enzm=enzm,
                                                                   bb_linkers=bb_linkers)
                        if full_comb is None:
                            continue
                        if score[0] == 1:
                            return full_comb, list(dropped_indices)
                        if fallback_score is None or score > fallback_score:
                            fallback_comb = full_comb
                            fallback_score = score
                            fallback_dropped = list(dropped_indices)
    if fallback_comb is not None:
        return fallback_comb, fallback_dropped
    return None, []


def tigRNA_subset_rescue_optimization(parts_list, free_overhangsets, rescue_overhangs, enzm='bsai', bb_linkers=['tgcc','gttt']):
    '''Compatibility wrapper for tigRNA subset rescue.'''

    return subset_rescue_optimization(parts_list, free_overhangsets, rescue_overhangs, 'tigRNA', enzm=enzm, bb_linkers=bb_linkers)

def primer_construct_match(primer, polycistron_sequence, enzm_seq, strand=1):
    '''
    Returns the primer segment that matches the final assembled polycistron.

    The Type IIS recognition site is excluded before matching. The 4 bp Golden
    Gate overhang is retained because it matches the final GenBank construct.
    '''

    primer = primer.lower()
    polycistron_sequence = polycistron_sequence.lower()
    trimmed = primer[len(enzm_seq):]
    if strand == -1:
        trimmed = reverse_complement(trimmed)

    for length in range(len(trimmed), 3, -1):
        for start in range(0, len(trimmed) - length + 1):
            match = trimmed[start:start + length]
            strt = polycistron_sequence.find(match)
            if strt != -1:
                return strt, strt + length
    return None

def tas_guide_design(sequence, system='TasH', edge_5=None, loop=None, edge_3=None, max_guides=20, min_gc=25, max_gc=75, exact_spacer=None):
    '''
    Designs programmable split-spacer tigRNAs for TIGR-Tas systems.

    The TIGR-Tas papers describe mature tigRNAs as edge-repeat fragment,
    spacer A, loop repeat, spacer B, edge-repeat fragment. TasA/FpTIGR uses
    two 9 nt spacers, while TasH multiplexing uses two 8 nt spacers. No
    PAM/TAM is required, so this scans all possible target windows on both
    orientations.

    When ``exact_spacer`` is provided, it is interpreted as one or more exact
    target windows instead of scanning the longer target sequence. Entries may
    be separated by whitespace, commas, semicolons or vertical bars. Exact
    target windows must be 18 bp for TasA and 16 bp for TasH.
    
    :param sequence: Target sequence to scan on both strands
    :type sequence: str
    :param system: Tas system name, currently 'TasA' or 'TasH'
    :type system: str, optional
    :param edge_5: Optional custom 5 prime edge-repeat sequence
    :type edge_5: str, optional
    :param loop: Optional custom loop-repeat sequence
    :type loop: str, optional
    :param edge_3: Optional custom 3 prime edge-repeat sequence
    :type edge_3: str, optional
    :param max_guides: Maximum number of candidates returned from scanning mode
    :type max_guides: int or str, optional
    :param min_gc: Minimum target-window GC percentage for scanning mode
    :type min_gc: float or str, optional
    :param max_gc: Maximum target-window GC percentage for scanning mode
    :type max_gc: float or str, optional
    :param exact_spacer: Optional exact target window or windows for multiplexing
    :type exact_spacer: str, optional
    :returns: Candidate list and optional warning message
    :rtype: tuple
    '''

    if sequence == '' and (exact_spacer is None or exact_spacer == ''):
        raise InvalidUsage("No exact target window input", status_code=400, payload={'pge': 'tas_generation.html', 'box': 'exact_spacer'})
    sequence = sequence.replace(' ', '').replace('\r\n', '').replace('\n', '').upper()
    if sequence != '' and re.search(r'^[ACGT]*$', sequence) is None:
        raise InvalidUsage("Invalid sequence input", status_code=400, payload={'pge': 'tas_generation.html', 'box': 'sequence'})
    if system not in TAS_SYSTEMS:
        raise InvalidUsage("Invalid Tas system", status_code=400, payload={'pge': 'tas_generation.html', 'box': 'system'})

    config = TAS_SYSTEMS[system]
    spacer_len = config['spacer_len']
    edge_5 = (edge_5 if edge_5 is not None and edge_5 != '' else config['edge_5']).replace(' ', '').upper()
    loop = (loop if loop is not None and loop != '' else config['loop']).replace(' ', '').upper()
    edge_3 = (edge_3 if edge_3 is not None and edge_3 != '' else config['edge_3']).replace(' ', '').upper()
    scaffold = edge_5 + loop + edge_3
    if re.search(r'^[ACGTU]*$', scaffold) is None:
        raise InvalidUsage("Invalid scaffold input", status_code=400, payload={'pge': 'tas_generation.html', 'box': 'scaffold'})
    expected_len = len(scaffold) + 2 * spacer_len
    if system == 'TasA' and expected_len != 36:
        raise InvalidUsage("TasA scaffold and spacer lengths must produce a 36 nt mature tigRNA", status_code=400, payload={'pge': 'tas_generation.html', 'box': 'scaffold'})
    if system == 'TasH' and expected_len != 46:
        raise InvalidUsage("TasH scaffold and spacer lengths must produce a 46 nt mature tigRNA", status_code=400, payload={'pge': 'tas_generation.html', 'box': 'scaffold'})

    max_guides = int(max_guides)
    min_gc = float(min_gc)
    max_gc = float(max_gc)
    if max_guides < 1:
        raise InvalidUsage("Maximum guides must be at least 1", status_code=400, payload={'pge': 'tas_generation.html', 'box': 'max_guides'})

    window_len = 2 * spacer_len
    exact_spacer = '' if exact_spacer is None else exact_spacer.upper()
    if exact_spacer.strip() != '':
        exact_spacers = [s for s in re.split(r'[\s,;|]+', exact_spacer.strip()) if s != '']
        candidates = []
        for c,spacer in enumerate(exact_spacers):
            if re.search(r'^[ACGT]*$', spacer) is None:
                raise InvalidUsage("Invalid exact spacer input in entry " + str(c + 1), status_code=400, payload={'pge': 'tas_generation.html', 'box': 'exact_spacer'})
            if len(spacer) != window_len:
                raise InvalidUsage(system + " exact target window " + str(c + 1) + " must be " + str(window_len) + " bp", status_code=400, payload={'pge': 'tas_generation.html', 'box': 'exact_spacer'})
            candidates.append(make_tas_candidate(spacer, system, edge_5, loop, edge_3, index=c + 1, source='exact ' + str(c + 1)))
        return candidates, None

    candidates = []
    for strand, scan_seq, rev in [('forward', sequence, False), ('reverse', reverse_complement(sequence), True)]:
        for start in range(0, len(scan_seq) - window_len + 1):
            target = scan_seq[start:start + window_len]
            gc = gc_content(target)
            if gc < min_gc or gc > max_gc:
                continue

            if rev:
                genomic_start = len(sequence) - start - window_len
                genomic_end = len(sequence) - start
            else:
                genomic_start = start
                genomic_end = start + window_len
            candidates.append(make_tas_candidate(target, system, edge_5, loop, edge_3, strand, genomic_start, genomic_end, len(candidates) + 1))

    candidates = sorted(candidates, key=lambda x: (x['score'], x['start'], x['strand']))[:max_guides]
    msg = None
    if candidates == []:
        msg = 'No Tas guide candidates passed the current GC filters.'
    return candidates, msg


## Optimize the overhangs used in Golden Gate assembly
def golden_gate_optimization(parts_list, free_overhangsets, poltype_opt='ptg', enzm='bsai', bb_linkers=['tgcc','gttt']):
    '''
    Finds linkers from optimal linker sets in a provided parts list

    In ``tigRNA`` mode, the optimizer uses spacer A or spacer B from recognized
    TasA/TasH tigRNAs as the source region for internal Golden Gate overhangs.
    This avoids placing variable overhangs in fixed edge or loop repeat
    scaffold sequence while giving the linker search more flexibility.
    
    :param parts_list: An array of objects of class part
    :type parts_list: list
    :param free_overhangsets: An array of arrays containing 4 bp linkers
    :type free_overhangsets: list
    :param poltype_opt: The type of polycistronic architecture to use. Must be
        one of 'ptg', 'ca' or 'tigRNA', defaults to 'ptg'
    :type poltype_opt: str, optional
    
    :return: A tuple of sequences indicating the linker between two adjacent parts. The sequence spans from the start of the respective variable sequence to the end of the determined 4 bp linker. 
    :rtype: tuple
    '''
    
    for p in parts_list:
        p.sequence = p.sequence.lower()
        
    ## Write all variable sequences in same order in a new list
    oh_list = []
    if poltype_opt=='ptg':
        for x in parts_list:
            if x.type == 'pegRNA':
                oh_list.append([x.sequence[:x.sequence.find(scaffld)],x.sequence[x.sequence.find(scaffld)+len(scaffld):]])
            elif x.type == 'gRNA':
                oh_list.append(x.sequence[:x.sequence.find(scaffld)])
            elif x.type == 'smRNA':
                oh_list.append(x.sequence[:-len(tRNA)])
            else: # if part is tRNA
                oh_list.append('plc')
    elif poltype_opt=='ca':
        for x in parts_list:
            oh_list.append(x.sequence[x.sequence.find(DR)+len(DR):])
    elif poltype_opt=='tigRNA':
        tigRNA_ranges = []
        for x in parts_list:
            oh_list.append(x.sequence)
            ranges = tigRNA_spacer_ranges(x.sequence)
            tigRNA_ranges.append(ranges)
    
    ## Starting in the middle of the variable sequences and moving outwards, find overhang combinations
    for cov in range(2,max([int(np.ceil(np.true_divide(len(prt),2))) for prt in flattn(oh_list)])):
        
        cov_list = []
        for oh in oh_list:
            if type(oh) is list:
                o_list = []
                for o in oh:
                    if 2*cov >= len(o):
                        o_list.append(o)
                    else:
                        o_list.append(o[int(np.ceil(np.true_divide(len(o),2)))-cov:int(np.ceil(np.true_divide(len(o),2)))+cov])
                cov_list.append(o_list)
            else:
                if 2*cov >= len(oh):
                    cov_list.append(oh)
                else:
                    cov_list.append(oh[int(np.ceil(np.true_divide(len(oh),2)))-cov:int(np.ceil(np.true_divide(len(oh),2)))+cov])
                    
        ## Find possible overhang combinations
        if poltype_opt=='ptg':
            for golden_gate_overhangs in free_overhangsets:
                seq_matches = []
                for x in range(len(parts_list)-1):
                    seq_matches.append([])
                    for overhang in golden_gate_overhangs:
                        if parts_list[x].type == 'pegRNA':
                            if overhang in cov_list[x][1]:
                                seq_matches[x].append(oh_list[x][1][:oh_list[x][1].find(overhang)+4])
                        elif parts_list[x].type == 'gRNA':
                            if parts_list[x+1].type == 'gRNA':
                                if overhang in cov_list[x+1]:
                                    seq_matches[x].append(oh_list[x+1][:oh_list[x+1].find(overhang)+4])
                            elif parts_list[x+1].type == 'pegRNA':
                                if overhang in cov_list[x+1][0]:
                                    seq_matches[x].append(oh_list[x+1][0][:oh_list[x+1][0].find(overhang)+4])
                            elif parts_list[x+1].type == 'smRNA':
                                if overhang in cov_list[x+1]:
                                    seq_matches[x].append(oh_list[x+1][:oh_list[x+1].find(overhang)+4])
                        elif parts_list[x].type == 'smRNA':
                            if parts_list[x+1].type == 'gRNA':
                                if overhang in cov_list[x+1]:
                                    seq_matches[x].append(oh_list[x+1][:oh_list[x+1].find(overhang)+4])
                            elif parts_list[x+1].type == 'pegRNA':
                                if overhang in cov_list[x+1][0]:
                                    seq_matches[x].append(oh_list[x+1][0][:oh_list[x+1][0].find(overhang)+4])
                            elif parts_list[x+1].type == 'smRNA':
                                if overhang in cov_list[x+1]:
                                    seq_matches[x].append(oh_list[x+1][:oh_list[x+1].find(overhang)+4])
                        else: # if part is tRNA
                            if parts_list[x+1].type == 'gRNA':
                                if overhang in cov_list[x+1]:
                                    seq_matches[x].append(oh_list[x+1][:oh_list[x+1].find(overhang)+4])
                            elif parts_list[x+1].type == 'pegRNA':
                                if overhang in cov_list[x+1][0]:
                                    seq_matches[x].append(oh_list[x+1][0][:oh_list[x+1][0].find(overhang)+4])
                            elif parts_list[x+1].type == 'smRNA':
                                if overhang in cov_list[x+1]:
                                    seq_matches[x].append(oh_list[x+1][:oh_list[x+1].find(overhang)+4])
        elif poltype_opt=='ca':
            for s in free_overhangsets:
                golden_gate_overhangs=s
                seq_matches = []
                for x in range(len(parts_list)-1):
                    seq_matches.append([])
                    for overhang in golden_gate_overhangs:
                        if overhang in cov_list[x]:
                            seq_matches[x].append(oh_list[x][:oh_list[x].find(overhang)+4])
        elif poltype_opt=='tigRNA':
            for s in free_overhangsets:
                golden_gate_overhangs=s
                seq_matches = []
                for x in range(len(parts_list)-1):
                    seq_matches.append([])
                    for overhang in golden_gate_overhangs:
                        for r_start,r_end in tigRNA_ranges[x]:
                            spacer_seq = oh_list[x][r_start:r_end]
                            if 2*cov >= len(spacer_seq):
                                spacer_cov = spacer_seq
                            else:
                                spacer_mid = int(np.ceil(np.true_divide(len(spacer_seq),2)))
                                spacer_cov = spacer_seq[spacer_mid-cov:spacer_mid+cov]
                            if overhang in spacer_cov:
                                overhang_start = oh_list[x].find(overhang, r_start, r_end)
                                if overhang_start + 4 >= TIGRNA_MIN_SPLIT_LEN:
                                    seq_matches[x].append(oh_list[x][:overhang_start+4])
                            
        # Copyright (c) 2019 Scott Weisberg
        combs = []
        for x in itertools.product(*seq_matches):
            combs.append(x)
        best_comb = None
        best_score = None
        for comb in combs:
            if len([c[-4:] for c in comb]) == len(set([c[-4:] for c in comb])):
                if poltype_opt in ['ca', 'tigRNA']:
                    score = overhang_combination_score(parts_list, comb, poltype_opt, enzm=enzm, bb_linkers=bb_linkers)
                    if score[0] < 0:
                        continue
                    if best_score is None or score > best_score:
                        best_comb = comb
                        best_score = score
                else:
                    return comb
        if best_comb is not None and best_score[0] == 1:
            return best_comb
        if best_comb is not None and ('fallback_comb' not in locals() or best_score > fallback_score):
            fallback_comb = best_comb
            fallback_score = best_score
    # If there are no possible combinations
    if 'fallback_comb' in locals():
        return fallback_comb
    return None


# Copyright (c) 2019 Scott Weisberg
## Perform scarless Golden Gate assembly computation with provided parts
def scarless_gg(parts_list, tm_range=[55,65], max_ann_len=30, bb_linkers=['tgcc','gttt'], ad_linkers=[], poltype='ptg', enzm='bsai'):
    '''
    Uses a list of desired parts and additional arguments to compute a corresponsing PTG. Returns a list of newly computed parts and their primers which can be used to generate the PTG.

    For ``tigRNA`` assemblies, the input parts are mature tigRNAs. Oligos are
    designed as scarless extension pairs around Golden Gate overhangs selected
    from spacer B. Primer Tm values reported on the output fragments are the
    Tm of the shared forward/reverse oligo overlap used for polymerase fill-in.
    
    :param parts_list: An array of objects of class part, to be included in a PTG
    :type parts_list: list
    :param tm_range: An array of integers of length 2. Temperature range to aim for during primer optimization, defaults to [55,65]
    :type tm_range: list, optional
    :param max_ann_len: The maximal annealing length of the static part of the primer, defaults to 30
    :type max_ann_len: int, optional
    :param bb_linkers: Linkers of the destination plasmid flanking the final PTG, defaults to ['tgcc','gttt']
    :type bb_linkers: list, optional
    :param ad_linkers: Additional linkers in the destination plasmid, defaults to []
    :type ad_linkers: list, optional
    :param poltype: Type of polycistronic architecture to use. Must be one of
        'ptg', 'ca' or 'tigRNA', defaults to 'ptg'
    :type poltype: str, optional
    :param enzm: Type II restriction enzyme to use for the Golden Gate assembly. Defaults to 'bsai'
    :type enzm: str, optional
    
    :return: Returns three objects. (1) A list of newly computed parts of class Part, (2) A list of features of class SeqRecord, (3) An error or warning message
    :rtype: list, list, str
    '''
    
    ## Catching errors
    for lnk in bb_linkers+ad_linkers:
        if len(lnk) != 4 or re.search(r'^[ACGTacgt]*$', lnk) is None:
            raise InvalidUsage("Invalid linker input", status_code=400, payload={'pge': 'sequence.html', 'box': 'link'})
    
    
    ## Setting up variables
    bb_linkers = [i.lower() for i in bb_linkers]
    ad_linkers = [i.lower() for i in ad_linkers]
    enzms={'bsai': ['gaggtctcg', 'cgagacctc'], 'bsmbi': ['tgcgtctca', 'tgagacgca'], 'btgzi': ['ctgcgatggagtatgtta', 'taacatactccatcgcag'], 'bbsi': ['ttgaagactt', 'aagtcttcaa']} #templates found in pUU080 (bsai), pUPD2 (bsmbi), Ortega-Escalante et al. 2018 (btgzi), Kun (bsi)
    
    
    ## Initiating the polycistron object
    polycistron = Polycistron()
    # must overwrite with empty list because features list would accumulate across runs in same session through append command
    polycistron.features = [] 
    polycistron.oligos = []
    polycistron.overhang_selection_mode = 'optimal'
    polycistron.overhang_warning = ''
    polycistron.selected_overhangs = []
    polycistron.overhang_statuses = []

    
    ## Go through parts and write all known annotations into list   
    mmry = len(enzms[enzm][0])+4
    polycistron.sequence = enzms[enzm][0] + reverse_complement(bb_linkers[0])
    
    for part in parts_list:
        part.sequence = part.sequence.lower()
        part.template_sequence = part.sequence
        part_label = part.name + ' ' + part.type
        unit_note = 'PolyGEN array unit ' + part.name + '; RNA type=' + part.type
        polycistron.sequence += part.sequence
        add_feature(polycistron, mmry, mmry+len(part.sequence), part.type,
                    label=part_label,
                    note=unit_note)
        if part.type == 'pegRNA':
            pbs_len = getattr(part, 'pbs_len', pegRNA_pbs_length(part.sequence))
            spacer_end = mmry+part.sequence.find(scaffld)
            rt_start = mmry+part.sequence.find(scaffld)+len(scaffld)
            add_feature(polycistron, mmry, spacer_end, 'spacer',
                        label=part.name + ' spacer',
                        note=part.name + ' pegRNA spacer; sequence=' + part.sequence[:part.sequence.find(scaffld)].upper())
            add_feature(polycistron, rt_start, mmry+len(part.sequence)-pbs_len, 'RT template',
                        label=part.name + ' RT template',
                        note=part.name + ' pegRNA reverse-transcription template')
            add_feature(polycistron, mmry+len(part.sequence)-pbs_len, mmry+len(part.sequence), 'PBS',
                        label=part.name + ' PBS',
                        note=part.name + ' pegRNA primer binding site; length=' + str(pbs_len) + ' nt')
        elif part.type == 'gRNA':
            spacer = part.sequence[:part.sequence.find(scaffld)]
            add_feature(polycistron, mmry, mmry+part.sequence.find(scaffld), 'spacer',
                        label=part.name + ' spacer',
                        note=part.name + ' sgRNA spacer; sequence=' + spacer.upper())
            add_feature(polycistron, mmry+part.sequence.find(tRNA), mmry+len(part.sequence), 'tRNA',
                        label=part.name + ' tRNA',
                        note='tRNA processing element for ' + part.name)
        elif part.type == 'smRNA':
            smrna = part.sequence[:-len(tRNA)]
            add_feature(polycistron, mmry, mmry+len(part.sequence)-len(tRNA), 'smRNA',
                        label=part.name + ' smRNA cargo',
                        note=part.name + ' small RNA cargo; sequence=' + smrna.upper())
            add_feature(polycistron, mmry+len(part.sequence)-len(tRNA), mmry+len(part.sequence), 'tRNA',
                        label=part.name + ' tRNA',
                        note='tRNA processing element for ' + part.name)
        elif part.type == 'tigRNA':
            shared_edge = getattr(part, 'shared_edge_5', None)
            shared_edge_start = mmry - len(shared_edge) if shared_edge and getattr(part, 'full_tigRNA_unit', None) else None
            shared_edge_3 = getattr(part, 'shared_edge_3', None)
            shared_edge_3_start = mmry - len(shared_edge_3) if shared_edge_3 else None
            annotate_tigRNA_features(polycistron, mmry, part.sequence,
                                     shared_edge_5_start=shared_edge_start,
                                     shared_edge_5=shared_edge,
                                     shared_edge_3_start=shared_edge_3_start,
                                     shared_edge_3=shared_edge_3)
        elif part.type == 'DR':
            add_feature(polycistron, mmry, mmry+len(part.sequence), 'direct_repeat',
                        label=part.name + ' terminal direct repeat',
                        note='Terminal Cpf1/Cas12a direct repeat')
        if part.type == 'crRNA':
            spacer = part.sequence[part.sequence.find(DR)+len(DR):]
            add_feature(polycistron, mmry, mmry+part.sequence.find(DR)+len(DR), 'direct_repeat',
                        label=part.name + ' direct repeat',
                        note='Cpf1/Cas12a direct repeat for ' + part.name)
            add_feature(polycistron, mmry+part.sequence.find(DR)+len(DR), mmry+len(part.sequence), 'spacer',
                        label=part.name + ' spacer',
                        note=part.name + ' crRNA spacer; sequence=' + spacer.upper())
        mmry += len(part.sequence)
    polycistron.sequence += bb_linkers[1] + enzms[enzm][1]

    ## Iterate through overhang sets with increasing size until fitting one is found
    gg_opt = None
    breakit = False
    collection_rescue_indices = set()
    tm_rescue_indices = []
    existing_overhangs = ad_linkers+bb_linkers
    compatible_free_overhangsets = []
    for p in range(10,51):
        free_overhangsets = []
        for q in range(5):
            with open(os.path.join(OVERHANGSET_DIR, 'setsof%s.csv'%p),'r') as f:
                reader = csv.reader(f, delimiter=",")
                sets = list(reader)[1:]
                
            temp = []
            for s in sets:
                if len(s) != 0:
                    temp.append(s)

            # Only grab sets that include all existing overhangs and delete the existing from the set
            st = [x.lower() for x in temp[q]]
            if all(i in st or reverse_complement(i) in st for i in existing_overhangs):
                free_overhangsets.append([i for i in [x.lower() for x in temp[q]] if i not in ad_linkers+bb_linkers])
            else:
                continue
        if free_overhangsets:
            compatible_free_overhangsets += free_overhangsets
            gg_opt = golden_gate_optimization(parts_list, free_overhangsets, poltype, enzm=enzm, bb_linkers=bb_linkers)
        if gg_opt is not None:
            breakit = True
        if breakit:
            break
        
    # No sets include all existing overhangs
    if gg_opt is None:
            
        overhangs_sublist = []
        for l in range(len(existing_overhangs)-1,0,-1):
            overhangs_subsublist = list(itertools.combinations(existing_overhangs,l))
            overhangs_sublist += overhangs_subsublist
            
        for sublist in overhangs_sublist:
            for p in range(10,51):
                free_overhangsets = []
                for q in range(5):
                    with open(os.path.join(OVERHANGSET_DIR, 'setsof%s.csv'%p),'r') as f:
                        reader = csv.reader(f, delimiter=",")
                        sets = list(reader)[1:]

                    temp = []
                    for s in sets:
                        if len(s) != 0:
                            temp.append(s)

                    # Only grab sets that include all linkers in the current sublist and delete the existing from the set
                    if all(i in [x.lower() for x in temp[q]] for i in sublist):
                        free_overhangsets.append([i for i in [x.lower() for x in temp[q]] if i not in ad_linkers+bb_linkers])
                    else:
                        continue

                if free_overhangsets:
                    gg_opt = golden_gate_optimization(parts_list, free_overhangsets, poltype, enzm=enzm, bb_linkers=bb_linkers)
                if gg_opt is not None:
                    breakit = True
                    
                    exist = ','.join(sublist)
                    nexist = ','.join(Diff(sublist, existing_overhangs))
                    polycistron.warning = 'The given combination of existing overhangs is not compatible with an optimal overhang set. Found the set including the largest possible fraction of existing overhangs ('+exist+'). The following overhangs were disregarded: '+nexist+'. There might be interference between overhangs.'
                if breakit:
                    break
            if breakit:
                break
        
    subset_rescue_allowed = (poltype in ['tigRNA', 'ca'] or
                             (poltype == 'ptg' and ptg_subset_rescue_allowed(parts_list)))
    if gg_opt is None and subset_rescue_allowed and compatible_free_overhangsets:
        disallowed = set(ad_linkers + bb_linkers + [reverse_complement(i) for i in ad_linkers + bb_linkers])
        curated = set(all_curated_overhangs())
        curated.update([reverse_complement(o) for o in list(curated)])
        rescue_overhangs = sorted([o for o in curated if o not in disallowed])
        subset_opt,subset_rescue_indices = subset_rescue_optimization(parts_list,
                                                                      compatible_free_overhangsets,
                                                                      rescue_overhangs,
                                                                      poltype,
                                                                      enzm=enzm,
                                                                      bb_linkers=bb_linkers)
        if subset_opt is not None:
            gg_opt = subset_opt
            collection_rescue_indices = set(subset_rescue_indices)
            polycistron.overhang_selection_mode = 'rescue'
            polycistron.overhang_warning = ('No validated optimal overhang set covered every internal junction. '
                                            'PolyGEN selected a validated optimal set for the largest compatible subset '
                                            'and used curated rescue overhangs for junction(s) ' +
                                            ','.join([str(i + 1) for i in subset_rescue_indices]) +
                                            '. Assembly should be experimentally validated.')
            polycistron.warning = polycistron.overhang_warning

    if gg_opt is None:
        disallowed = set(ad_linkers + bb_linkers + [reverse_complement(i) for i in ad_linkers + bb_linkers])
        curated = set(all_curated_overhangs())
        curated.update([reverse_complement(o) for o in list(curated)])
        rescue_overhangs = sorted([o for o in curated if o not in disallowed])
        rescue_opt = golden_gate_optimization(parts_list, [rescue_overhangs], poltype, enzm=enzm, bb_linkers=bb_linkers)
        if rescue_opt is None:
            raise InvalidUsage("No combination of optimal or rescue linkers possible for the provided existing linkers", status_code=400, payload={'pge': 'sequence.html', 'box': 'link'})
        gg_opt = rescue_opt
        collection_rescue_indices = set(range(len(gg_opt)))
        polycistron.overhang_selection_mode = 'rescue'
        polycistron.overhang_warning = 'Rescue overhang mode was used. The selected overhangs are individually present in curated overhang tables, but this exact overhang collection was not found as a validated optimal set. Assembly may fail and should be experimentally validated.'
        polycistron.warning = polycistron.overhang_warning
    if poltype in ['ca', 'tigRNA']:
        initial_overlap_tms = overhang_combination_overlap_tms(parts_list, gg_opt, poltype, enzm=enzm, bb_linkers=bb_linkers)
        if initial_overlap_tms != [] and min(initial_overlap_tms) < MIN_FILL_IN_OVERLAP_TM:
            disallowed = set(ad_linkers + bb_linkers + [reverse_complement(i) for i in ad_linkers + bb_linkers])
            curated = set(all_curated_overhangs())
            curated.update([reverse_complement(o) for o in list(curated)])
            allowed_tm_rescue = set([o for o in curated if o not in disallowed])
            improved_gg_opt, tm_rescue_indices = improve_overhang_combination_tms(parts_list, gg_opt, poltype, enzm=enzm, bb_linkers=bb_linkers, allowed_overhangs=allowed_tm_rescue)
            if tm_rescue_indices:
                gg_opt = improved_gg_opt
                if polycistron.overhang_selection_mode == 'optimal':
                    polycistron.overhang_selection_mode = 'rescue'
                tm_rescue_note = 'Tm-aware overhang rescue adjusted one or more split positions because the initial fill-in overlap Tm was below target.'
                if polycistron.overhang_warning:
                    polycistron.overhang_warning += ' ' + tm_rescue_note
                else:
                    polycistron.overhang_warning = tm_rescue_note
                polycistron.warning = (polycistron.warning + ' ' + tm_rescue_note) if polycistron.warning else tm_rescue_note

    polycistron.selected_overhangs = [g[-4:] for g in gg_opt]
    polycistron.overhang_statuses = []
    for idx,overhang in enumerate(polycistron.selected_overhangs):
        if idx in tm_rescue_indices:
            status = 'tm_rescue'
            note = 'This overhang split was adjusted by Tm-aware rescue to improve the oligo-pair overlap Tm.'
        elif idx in collection_rescue_indices:
            status = 'collection_rescue'
            note = 'This overhang is present in curated tables and was used as rescue because the full tigRNA collection did not fit one validated optimal set.'
        else:
            status = 'optimal'
            note = 'This overhang was selected from a validated optimal overhang set.'
        polycistron.overhang_statuses.append({'index': idx + 1,
                                              'overhang': overhang,
                                              'status': status,
                                              'note': note})
    polycistron.fill_in_overlap_tms = overhang_combination_overlap_tms(parts_list, gg_opt, poltype, enzm=enzm, bb_linkers=bb_linkers) if poltype == 'tigRNA' else []
    if polycistron.fill_in_overlap_tms != []:
        min_overlap_tm = min(polycistron.fill_in_overlap_tms)
        if min_overlap_tm < MIN_FILL_IN_OVERLAP_TM:
            tm_warning = ('Low fill-in overlap Tm warning. The best available overhang choice has a minimum '
                          'predicted oligo-pair overlap Tm of ' + str(round(min_overlap_tm, 1)) +
                          ' C, below the ' + str(MIN_FILL_IN_OVERLAP_TM) +
                          ' C target. Assembly may be less efficient for short fragments.')
            if polycistron.overhang_warning:
                polycistron.overhang_warning += ' ' + tm_warning
            else:
                polycistron.overhang_warning = tm_warning
            polycistron.warning = (polycistron.warning + ' ' + tm_warning) if polycistron.warning else tm_warning
    
    
    ## Modify sequences and design primers
    if poltype=='ptg':
        for i in range(len(parts_list)):
            
            # If current part is first part, define forward primer with left backbone overlap and reverse primer ordinarily
            if i == 0:
                parts_list[i].primer_forward = enzms[enzm][0] + reverse_complement(bb_linkers[0]) + parts_list[i].sequence[:max_ann_len]
                parts_list[i].sequence = enzms[enzm][0] + reverse_complement(bb_linkers[0]) + parts_list[i].sequence
                if parts_list[i].type == 'pegRNA':
                    parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[parts_list[i].sequence.find(scaffld.lower())+76-max_ann_len:parts_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])] + enzms[enzm][1])
                    parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i].sequence[parts_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])-4:] + parts_list[i+1].sequence[:max_ann_len]
                    parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i].sequence[parts_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])-4:] + parts_list[i+1].sequence
                    parts_list[i].sequence = parts_list[i].sequence[:parts_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])] + enzms[enzm][1]
                elif parts_list[i].type == 'gRNA':
                    parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[-max_ann_len:] + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1])
                    parts_list[i].sequence = parts_list[i].sequence + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1]
                    if parts_list[i+1].type == 'smRNA':
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(tRNA.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(scaffld.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                elif parts_list[i].type == 'smRNA':
                    parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[-max_ann_len:] + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1])
                    parts_list[i].sequence = parts_list[i].sequence + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1]
                    if parts_list[i+1].type == 'smRNA':
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(tRNA.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(scaffld.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                else: #part is tRNA
                    parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[-max_ann_len:] + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1])
                    parts_list[i].sequence = parts_list[i].sequence + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1]
                    if parts_list[i+1].type == 'smRNA':
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(tRNA.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(scaffld.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
            
        # If current part is last part, define reverse primer with right backbone overlap and forward primer ordinarily
            elif i == len(parts_list)-1:
                parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[-max_ann_len:] + bb_linkers[-1] + enzms[enzm][1])
                parts_list[i].sequence = parts_list[i].sequence + bb_linkers[-1] + enzms[enzm][1]
            
        # If current part is not first or last part, do ordinary computation
            else:
                if parts_list[i].type == 'pegRNA':
                    parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[parts_list[i].sequence.find(scaffld.lower())+76-max_ann_len:parts_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])] + enzms[enzm][1])
                    parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i].sequence[parts_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])-4:] + parts_list[i+1].sequence[:max_ann_len]
                    parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i].sequence[parts_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])-4:] + parts_list[i+1].sequence
                    parts_list[i].sequence = parts_list[i].sequence[:parts_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])] + enzms[enzm][1]
                elif parts_list[i].type == 'gRNA':
                    parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[-max_ann_len:] + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1])
                    parts_list[i].sequence = parts_list[i].sequence + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1]
                    if parts_list[i+1].type == 'smRNA':
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(tRNA.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(scaffld.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                elif parts_list[i].type == 'smRNA':
                    parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[-max_ann_len:] + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1])
                    parts_list[i].sequence = parts_list[i].sequence + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1]
                    if parts_list[i+1].type == 'smRNA':
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(tRNA.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(scaffld.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                else: #part is tRNA
                    parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[-max_ann_len:] + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1])
                    parts_list[i].sequence = parts_list[i].sequence + parts_list[i+1].sequence[:parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1]
                    if parts_list[i+1].type == 'smRNA':
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(tRNA.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:parts_list[i+1].sequence.find(scaffld.lower())+max_ann_len]
                        parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i+1].sequence[parts_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                            
    elif poltype=='ca':
        for i in range(len(parts_list)):
            
            # If current part is first part, define forward primer with left backbone overlap and reverse primer ordinarily
            if i == 0:
                parts_list[i].primer_forward = enzms[enzm][0] + reverse_complement(bb_linkers[0]) + DR
                parts_list[i].sequence = enzms[enzm][0] + reverse_complement(bb_linkers[0]) + parts_list[i].sequence
                
                parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[parts_list[i].sequence.find(DR):parts_list[i].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1])
                parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i].sequence[parts_list[i].sequence.find(gg_opt[i]) + len(gg_opt[i])-4:] + DR
                parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i].sequence[parts_list[i].sequence.find(gg_opt[i]) + len(gg_opt[i])-4:] + parts_list[i+1].sequence
                parts_list[i].sequence = parts_list[i].sequence[:parts_list[i].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1]

            # If current part is last part, define reverse primer with right backbone overlap and forward primer ordinarily
            elif i == len(parts_list)-1:
                parts_list[i].primer_reverse = reverse_complement(DR + bb_linkers[-1] + enzms[enzm][1])
                parts_list[i].sequence = parts_list[i].sequence + bb_linkers[-1] + enzms[enzm][1]

            # If current part is not first or last part, do ordinary computation
            else:
                parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[parts_list[i].sequence.find(DR):parts_list[i].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1])
                parts_list[i+1].primer_forward = enzms[enzm][0] + parts_list[i].sequence[parts_list[i].sequence.find(gg_opt[i]) + len(gg_opt[i])-4:] + DR
                parts_list[i+1].sequence = enzms[enzm][0] + parts_list[i].sequence[parts_list[i].sequence.find(gg_opt[i]) + len(gg_opt[i])-4:] + parts_list[i+1].sequence
                parts_list[i].sequence = parts_list[i].sequence[:parts_list[i].sequence.find(gg_opt[i])+len(gg_opt[i])] + enzms[enzm][1]

    elif poltype=='tigRNA':
        for i in range(len(parts_list)):
            if i == 0:
                parts_list[i].primer_forward = enzms[enzm][0] + reverse_complement(bb_linkers[0]) + parts_list[i].sequence[:max_ann_len]
                parts_list[i].sequence = enzms[enzm][0] + reverse_complement(bb_linkers[0]) + parts_list[i].sequence

                if len(parts_list) > 1:
                    current_start = parts_list[i].sequence.find(parts_list[i].template_sequence)
                    current_cut = parts_list[i].sequence.find(gg_opt[i]) + len(gg_opt[i])
                    current_tail = parts_list[i].sequence[current_cut-4:]
                    if getattr(parts_list[i+1], 'full_tigRNA_unit', None) and getattr(parts_list[i+1], 'shared_edge_5', None) and parts_list[i+1].shared_edge_5.lower() not in current_tail:
                        raise InvalidUsage("Shared TasH edge repeat was not retained in the oligo overlap", status_code=400, payload={'pge': 'sequence.html', 'box': 'sequence_spacers'})
                    if getattr(parts_list[i+1], 'shared_edge_junction', None) and parts_list[i+1].shared_edge_junction.lower() not in current_tail + parts_list[i+1].sequence[:len(parts_list[i+1].shared_edge_5)]:
                        raise InvalidUsage("Shared TasA edge junction was not retained in the oligo overlap", status_code=400, payload={'pge': 'sequence.html', 'box': 'sequence_spacers'})
                    parts_list[i+1].shared_edge_overlap = current_tail
                    parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[current_start:current_cut] + enzms[enzm][1])
                    parts_list[i+1].primer_forward = enzms[enzm][0] + current_tail + parts_list[i+1].sequence[:max_ann_len]
                    parts_list[i+1].sequence = enzms[enzm][0] + current_tail + parts_list[i+1].sequence
                    parts_list[i].sequence = parts_list[i].sequence[:current_cut] + enzms[enzm][1]
                else:
                    parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[-max_ann_len:] + bb_linkers[-1] + enzms[enzm][1])
                    parts_list[i].sequence = parts_list[i].sequence + bb_linkers[-1] + enzms[enzm][1]

            elif i == len(parts_list)-1:
                parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[-max_ann_len:] + bb_linkers[-1] + enzms[enzm][1])
                parts_list[i].sequence = parts_list[i].sequence + bb_linkers[-1] + enzms[enzm][1]

            else:
                current_start = parts_list[i].sequence.find(parts_list[i].template_sequence)
                current_cut = parts_list[i].sequence.find(gg_opt[i]) + len(gg_opt[i])
                current_tail = parts_list[i].sequence[current_cut-4:]
                if getattr(parts_list[i+1], 'full_tigRNA_unit', None) and getattr(parts_list[i+1], 'shared_edge_5', None) and parts_list[i+1].shared_edge_5.lower() not in current_tail:
                    raise InvalidUsage("Shared TasH edge repeat was not retained in the oligo overlap", status_code=400, payload={'pge': 'sequence.html', 'box': 'sequence_spacers'})
                if getattr(parts_list[i+1], 'shared_edge_junction', None) and parts_list[i+1].shared_edge_junction.lower() not in current_tail + parts_list[i+1].sequence[:len(parts_list[i+1].shared_edge_5)]:
                    raise InvalidUsage("Shared TasA edge junction was not retained in the oligo overlap", status_code=400, payload={'pge': 'sequence.html', 'box': 'sequence_spacers'})
                parts_list[i+1].shared_edge_overlap = current_tail
                parts_list[i].primer_reverse = reverse_complement(parts_list[i].sequence[current_start:current_cut] + enzms[enzm][1])
                parts_list[i+1].primer_forward = enzms[enzm][0] + current_tail + parts_list[i+1].sequence[:max_ann_len]
                parts_list[i+1].sequence = enzms[enzm][0] + current_tail + parts_list[i+1].sequence
                parts_list[i].sequence = parts_list[i].sequence[:current_cut] + enzms[enzm][1]
                            
    ## Format enzyme cutting site in CAPITAL letters
    for o in parts_list:
        o.primer_forward = o.primer_forward[:len(enzms[enzm][0])].lower() + o.primer_forward[len(enzms[enzm][0]):len(enzms[enzm][0])+4].upper() + o.primer_forward[len(enzms[enzm][0])+4:].lower()
        o.primer_reverse = o.primer_reverse[:len(enzms[enzm][1])].lower() + o.primer_reverse[len(enzms[enzm][1]):len(enzms[enzm][1])+4].upper() + o.primer_reverse[len(enzms[enzm][1])+4:].lower()
    
    ## Assign localisation of parts in whole sequence
    mmry = 0
    for c,part in enumerate(parts_list):
        if c == 0:
            part.localisation = [mmry, mmry + len(part.sequence) - len(enzms[enzm][1])]
            mmry += len(part.sequence) - len(enzms[enzm][1]) - 4
        elif c == len(parts_list)-1:
            part.localisation = [mmry, len(polycistron.sequence)]
            mmry += len(polycistron.sequence)
        else:
            part.localisation = [mmry, mmry + len(part.sequence) - len(enzms[enzm][0]) - len(enzms[enzm][1])]
            mmry += len(part.sequence) - len(enzms[enzm][0]) - len(enzms[enzm][1]) - 4
    
    polycistron.parts = parts_list
    
    
    if poltype == 'tigRNA':
        # tigRNA oligos contain required scarless extension sequence; do not trim them.
        for part in polycistron.parts:
            overlap_tm = primer_pair_overlap_tm(part.primer_forward, part.primer_reverse)
            part.primer_forward_tm = overlap_tm
            part.primer_reverse_tm = overlap_tm
        return polycistron

    
    ## Optimise primer Tm    
    for c,prmr in enumerate(flattn([[part.primer_forward, part.primer_reverse] for part in polycistron.parts])):
        
        prmrRest = prmr[:-max_ann_len]
        prmrRestLen = len(prmrRest)
    
        scr = [1 if i in ['C','G','g','c'] else 0 for i in prmr] # Extract score from sequence for finding GC-clamps
        optns = [prmr[:i] for i in range(len(prmr),len(prmr)-(max_ann_len-18+1), -1) if sum(scr[i-5:i]) in [1,2]] # Find all GC-clamps
        checkStatic = [optn for optn in optns if tm_range[0] <= mt.Tm_NN(optn[prmrRestLen:], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) <= tm_range[1]] # Find all GC-clamp options with static Tm inside provided range
        checkWhole = [optn for optn in checkStatic if tm_range[0] <= mt.Tm_NN(optn, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) <= tm_range[1]] # Find all GC-clamp options with whole Tm inside provided range
        
        if checkWhole: # If there are GC-clamp options with whole Tm, use the shortest one to cut costs
            fin = checkWhole[-1]
        elif checkStatic: # If there are no GC-clamps with whole Tm but such with static Tm, use the shortest, to (1) cut costs and (2) get the Tm of the whole primer as close as possible to range.
            fin = checkStatic[-1]
        elif optns and mt.Tm_NN(optns[-1][prmrRestLen:], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) > tm_range[1]: # If all GC-clamps are above the range, use lowest one
            fin = optns[-1]
        elif optns and mt.Tm_NN(optns[0][prmrRestLen:], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) < tm_range[0]: # If all GC-clamps are below the range, use highest one
            fin = optns[0]
        else: # If there are no GC-clamps, find the primer with static Tm closest to bottom of range for cutting costs
            countUp = prmrRestLen+18
            lastResort = prmr[prmrRestLen:countUp]
            while mt.Tm_NN(lastResort, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) < tm_range[0] and prmrRestLen + len(lastResort) < len(prmr):
                countUp += 1
                lastResort = prmr[prmrRestLen:countUp]
            fin = prmr[:countUp]
        
        if c%2 == 0: # Assign the found primer and corresponding static Tm to the respective part
            polycistron.parts[int(np.floor(c/2))].primer_forward = fin
            polycistron.parts[int(np.floor(c/2))].primer_forward_tm = mt.Tm_NN(fin[prmrRestLen:], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
        else:
            polycistron.parts[int(np.floor(c/2))].primer_reverse = fin
            polycistron.parts[int(np.floor(c/2))].primer_reverse_tm = mt.Tm_NN(fin[prmrRestLen:], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
            
    return polycistron


def pegbldr(sequence, edits, mode='PE2'):
    '''
    Builds the necesary pegRNA and gRNA to accurately introduce a point mutation using the prime editor system.

    PBS length is selected from 8-17 nt by choosing the PBS whose melting
    temperature is closest to 30 C.
    
    :param sequence: The original 5' -> 3' nucleic acid sequence consisting of [ACGTacgt] that is to be edited
    :type sequence: str
    :param edits: An array containing all desired edits in the form [[index of edit, new bases, type of edit],[...]]
    :type edits: list
    :param mode: The type of prime editing mode to design for. One of 'PE2' or 'PE3', defaults to 'PE3'
    :type mode: str, optional
    
    :return: A list of lists containing for each computed guide RNA the name, type, sequence and strand specifications
    :rtype: list
    '''
    
    ## Catching errors
    if sequence == '':
        raise InvalidUsage("No sequence input", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'sequence'})
    elif re.search(r'^[ACGTacgt]*$', sequence) is None:
        raise InvalidUsage("Invalid sequence input ", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'sequence'})
    
    if edits == [['']]:
        raise InvalidUsage("No edits input", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'edits'})
    for e in edits:
        if len(e) != 3:
            raise InvalidUsage("Invalid edits input syntax", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'edits'})
        elif len(e[0].split(',')) != len(e[1].split(',')):
            raise InvalidUsage("Inconsistent number of edits", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'edits'})
        elif e[2] not in ['mut', 'del', 'ins']:
            raise InvalidUsage("Invalid edit type", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'edits'})
        elif e[2] == 'mut' and (any([re.search(r'^[0-9]*$', n) is None for n in e[0].split(',')]) or any([re.search(r'^[ACGTacgt]*$', b) is None for b in e[1].split(',')])):
            raise InvalidUsage("Invalid specifications for edit type \'mut\'", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'edits'})
        elif e[2] == 'del' and any([re.search(r'^[0-9]*$', n) is None for n in e[0].split(',')+e[1].split(',')]):
            raise InvalidUsage("Invalid specifications for edit type \'del\'", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'edits'})
        elif e[2] == 'ins' and (any([re.search(r'^[0-9]*$', n) is None for n in e[0].split(',')]) or any([re.search(r'^[ACGTacgt]*$', b) is None for b in e[1].split(',')])):
            raise InvalidUsage("Invalid specifications for edit type \'ins\'", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'edits'})
    
    
    ## Setting up variables
    msg = None
    sequence = sequence.upper()
    plc_seq = sequence[:]
    plc_rev_comp = reverse_complement(plc_seq)
    
    
    ## Starting the main computation
    out = []
    for c,edt in enumerate(edits):
        
        seq = plc_seq[:]
        rev_comp = plc_rev_comp[:]
        
        inds = [int(i) for i in str(edt[0]).split(',')]
        changes = str(edt[1]).upper().split(',')
        if edt[2] == 'del':
            changes = [str(int(chng)+1) for chng in changes] # deletion stretches should include the ending index for intuitive usability

        
        ## Define pegRNA for editing
        allInd = inds + [int(chng) for chng in changes if edt[2] == 'del']
        maxInd = max(allInd)
        minInd = min(allInd)
        
        
        ## Find all PAM motifs in forward and reverse strand in respective possible regions
        pegPAMrgn_forw = seq[:minInd+6]                       # Define the region where PAM could possibly lie in forw strand
        pegPAMs_forw = re.finditer(r'(?=(.GG))', pegPAMrgn_forw) # Find all PAMs in region
        pegPAMrgn_rev = rev_comp[:len(seq)-1-maxInd+6]
        pegPAMs_rev = re.finditer(r'(?=(.GG))', pegPAMrgn_rev)
        pegPAMs_forw = list(pegPAMs_forw)
        pegPAMs_rev = list(pegPAMs_rev)
        pegPAMs_forw = [i.start() for i in pegPAMs_forw]
        pegPAMs_rev = [i.start() for i in pegPAMs_rev]
        
        
        ## Check if usable PAMs are present
        if pegPAMs_forw == pegPAMs_rev == []:
            raise InvalidUsage("There are no usable PAM motifs around the edit", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'sequence'})
        
        
        ## Check if there are PAMs in only one strand
        elif pegPAMs_forw == []:
            pegPAM_rev = pegPAMs_rev[np.argmin([abs(i-(len(seq)-2-max(inds))) for i in pegPAMs_rev])]
            pegPAM = pegPAM_rev
            pegPAM_strand = 'r'
            seq = rev_comp[:]
            rev_comp = sequence[:]
            inds = [len(seq)-i-1 for i in inds] # Recalculate mutation indices for rev strand. -1 because len(x) == ind(x[-1])+1
            changes = [complement(i) for i in changes]
        elif pegPAMs_rev == []:
            pegPAM_forw = pegPAMs_forw[np.argmin([abs(i-1-min(inds)) for i in pegPAMs_forw])]
            pegPAM = pegPAM_forw
            pegPAM_strand = 'f'
        
        
        ## If there are PAMs in both strands, choose the one closest to edit. inf is used to make sure the list is not empty but will never be chosen.
        else:
            forw_min = min([abs(i-1-inds[0]) for i in pegPAMs_forw] + [np.inf])
            rev_min = min([abs(i-(len(seq)-2-inds[-1])) for i in pegPAMs_rev] + [np.inf])
            pegPAM_forw = pegPAMs_forw[np.argmin([abs(i-1-inds[0]) for i in pegPAMs_forw])]
            pegPAM_rev = pegPAMs_rev[np.argmin([abs(i-(len(seq)-2-inds[-1])) for i in pegPAMs_rev])]
            if forw_min <= rev_min: # Of the closest PAM of each strand, which is closest to mutation? -1 because len(x) == ind(x[-1])+1
                pegPAM = pegPAM_forw
                pegPAM_strand = 'f'
            else:
                pegPAM = pegPAM_rev
                pegPAM_strand = 'r'
                seq = rev_comp[:]
                rev_comp = sequence[:]
                if edt[2] == 'mut':
                    inds = [len(seq)-i-1 for i in inds] # Recalculate mutation indices for rev strand. -1 because len(x) == ind(x[-1])+1
                    changes = [complement(i) for i in changes]
                elif edt[2] == 'ins':
                    inds = [len(seq)-i for i in inds]
                    changes = [complement(i) for i in changes]
                else:
                    changes = [int(chng) for chng in changes] # Change type of changes to int
                    changes_plc = changes[:] # create placeholder
                    changes = [len(seq)-i-1 for i in inds] # switch changes and inds
                    inds = [len(seq)-i-1 for i in changes_plc]
                                  
        if max(inds)-pegPAM > 30:
                msg = "There was no PAM motif in +/- 30 nt proximity of edit " + str(c) + ". Used the nearest one, which was " + str(max(inds)-pegPAM) + " bp away."

        if pegPAM-20 < 0:
            raise InvalidUsage("The provided sequence does not cover enough area around the edit", status_code=400, payload={'pge': 'peg_generation.html', 'box': 'sequence'})
        else:
            pegspacer = seq[pegPAM-20:pegPAM]                    # Spacer should be 20 nt in length and end at PAM
        PBS,pbs_len,pbs_temp = choose_pbs(seq, pegPAM) # PBS must be in opposite direction

        
        ## Calculate RT templates depending on type of edit
        if edt[2] == 'mut':                                      # Check if edit is point mutation
            
            pre_RT_len = max([13, max(inds)-(pegPAM-3)]) # Set default length of RT-template to 13 (recommended by Anzalone et al. 2019) or until edit if further
            post_RT_len = pre_RT_len + re.search(r'[AGT]', seq[pegPAM-3+pre_RT_len:]).start() # From default length find next D
            RT_templ = seq[pegPAM-3:pegPAM-3+post_RT_len+1]      # Retrieve RT-template
            RT_templ = [i for i in RT_templ]
            
            for c_pm,pm in enumerate(inds):                      # If several point mutations, go through each separately
                RT_templ[int(pm)-(pegPAM-3)] = changes[c_pm]     # Include edit in RT-template
                
        elif edt[2] == 'ins':
            
            pre_RT_len = max([13, max(inds)-(pegPAM-3)+6]) # template should have additional length 5' of insert to ensure binding
            post_RT_len = pre_RT_len + re.search(r'[AGT]', seq[pegPAM-3+pre_RT_len:]).start()
            RT_templ = seq[pegPAM-3:pegPAM-3+post_RT_len+1]
            RT_templ = [i for i in RT_templ]
            for c_pm,pm in enumerate(inds):
                RT_templ = RT_templ[:int(pm)-(pegPAM-3)] + [l for l in changes[c_pm]] + RT_templ[int(pm)-(pegPAM-3):]
                
        elif edt[2] == 'del':
            
            changes = [int(chng) for chng in changes]
            len_deltns = 0
            for c_pm,pm in enumerate(inds):
                len_deltns += changes[c_pm]-pm
            pre_RT_len = max([13, max(changes)-(pegPAM-3)+6+len_deltns])
            post_RT_len = pre_RT_len + re.search(r'[AGT]', seq[pegPAM-3+pre_RT_len:]).start()
            RT_templ = seq[pegPAM-3:pegPAM-3+post_RT_len+1]
            RT_templ = [i for i in RT_templ]
            
            mmry = 0
            for c_pm,pm in enumerate(inds):
                del RT_templ[pm-(pegPAM-3)-mmry:changes[c_pm]-(pegPAM-3)-mmry]
                mmry += changes[c_pm]-pm
                
        RT_templ = ''.join(RT_templ)
        
        RT_templ = reverse_complement(RT_templ)        # RT-template must be in opposite direction
        
        pegRNA = pegspacer + scaffld.upper() + RT_templ + PBS
        PEG_PBS_LENGTHS[pegRNA.upper()] = pbs_len
        
        out.append(['pegRNA'+str(c), 'pegRNA', pegRNA, pegPAM_strand])

        
        ## Define gRNA for PE3
        if mode == 'PE3':
            
            nuseq = [i for i in seq]
            mmry = 0

            for c2,ind in enumerate(inds):
                if edt[2] == 'mut':
                    nuseq[ind] = changes[c2]
                elif edt[2] == 'ins':
                    nuseq = nuseq[:ind] + [b for b in changes[c2]] + nuseq[ind:]
                elif edt[2] == 'del':
                    del nuseq[ind:changes[c2]]
            nuseq = ''.join(nuseq)
            
            gPAMrgn = reverse_complement(nuseq[pegPAM-3:])     # gRNA must bind to other strand somewhere downstream of pegPAM
            gPAMs = re.finditer(r'(?=(.GG))', gPAMrgn)     # find all PAMs in that region
            gPAMs = np.array([i.start() for i in gPAMs])
            gPAM = gPAMs[abs(np.array(gPAMs)-(len(gPAMrgn)-44)).argmin()] # Find PAM closest to 47 nt downstream of pegPAM nick (the gRNA nick should be 50 nt downstream of pegPAM nick)

            gspacer = gPAMrgn[gPAM-20:gPAM]                # Define spacer

            gRNA = gspacer

            out.append(['gRNA'+str(c), 'gRNA', gRNA, pegPAM_strand])
    
    return out,msg


def PTGbldr(name, inserts, poltype='ptg'):
    '''
    Takes all desired parts of PTG GG assembly and gives out the respective inserts for PTG. During the 
    process, each part is appended with the same tRNA. The unit of part and tRNA are then treated as one insert.

    Three polycistron architectures are supported:

    * ``ptg`` accepts ``gRNA``, ``pegRNA`` and ``smRNA`` inserts and adds the
      tRNA/scaffold context required by the PTG architecture.
    * ``ca`` accepts ``crRNA`` inserts and adds the Cpf1 direct repeat.
    * ``tigRNA`` accepts only ``tigRNA`` inserts. It treats every sequence as a
      mature TIGR-Tas tigRNA and does not add tRNAs or a Cpf1 direct repeat.
    
    :param inserts: List of parts of the form [['type', 'sequence'],[...]]
    :type inserts: list
    :param poltype: Type of polycistronic architecture to use. Must be one of
        'ptg', 'ca' or 'tigRNA', defaults to 'ptg'
    :type poltype: str, optional
    '''
    
    ## Catching errors
    for e in inserts:
        if len(e) != 2:
            raise InvalidUsage("Invalid input syntax", status_code=400, payload={'pge': 'sequence.html', 'box': 'sequence_spacers'})
        elif e[0] not in ['gRNA', 'pegRNA', 'smRNA', 'crRNA', 'tigRNA']:
            raise InvalidUsage("Invalid RNA type", status_code=400, payload={'pge': 'sequence.html', 'box': 'sequence_spacers'})
        elif re.search(r'^[ACGTacgt]*$', e[1]) is None:
            raise InvalidUsage("Invalid sequence input", status_code=400, payload={'pge': 'sequence.html', 'box': 'sequence_spacers'})
    
    
    if poltype=='ptg':
        PTG_parts = []
        PTG_parts.append(Part(name+'_0', 'tRNA', tRNA))
        c = 1
    
        ## Take each coding sequence of the output and append it with a tRNA to the list
        for prt in inserts:
            if prt[0] not in ['pegRNA', 'gRNA', 'smRNA']:
                raise InvalidUsage("crRNAs and tigRNAs should be processed in their matching polycistron type", status_code=400, payload={'pge': 'sequence.html', 'box': 'poltype_input'})
            if prt[0] == 'pegRNA':
                peg_part = Part(name+'_'+str(c), 'pegRNA', str(prt[1]))
                if str(prt[1]).upper() in PEG_PBS_LENGTHS:
                    peg_part.pbs_len = PEG_PBS_LENGTHS[str(prt[1]).upper()]
                PTG_parts.append(peg_part)
                PTG_parts.append(Part(name+'_'+str(c+1), 'tRNA', tRNA))
                c += 2
            elif prt[0] == 'gRNA':
                PTG_parts.append(Part(name+'_'+str(c), 'gRNA', str(prt[1]) + scaffld + tRNA))
                c += 1
            elif prt[0] == 'smRNA':
                PTG_parts.append(Part(name+'_'+str(c), 'smRNA', str(prt[1]) + tRNA))
                c += 1
    
        return PTG_parts
        
    elif poltype=='ca':
        CA_parts = []
        c = 0
        
        for c,prt in enumerate(inserts):
            if prt[0] != 'crRNA':
                raise InvalidUsage("CA can only process crRNAs", status_code=400, payload={'pge': 'sequence.html', 'box': 'poltype_input'})
            CA_parts.append(Part(name+'_'+str(c), prt[0], DR + str(prt[1])))
            c += 1
        
        CA_parts.append(Part(name+'_'+str(c), 'DR', DR))
	            
        return CA_parts

    elif poltype=='tigRNA':
        tigRNA_parts = []
        previous_edge_3 = None

        for c,prt in enumerate(inserts):
            if prt[0] != 'tigRNA':
                raise InvalidUsage("tigRNA mode can only process tigRNAs", status_code=400, payload={'pge': 'sequence.html', 'box': 'poltype_input'})
            sequence = str(prt[1])
            part = Part(name+'_'+str(c), prt[0], sequence)
            for system,config in TAS_SYSTEMS.items():
                edge_5 = config['edge_5'].replace('U', 'T')
                edge_3 = config['edge_3'].replace('U', 'T')
                if previous_edge_3 == edge_5 and sequence.upper().startswith(edge_5.upper()) and edge_5.upper() == edge_3.upper():
                    part.full_tigRNA_unit = sequence
                    part.sequence = sequence[len(edge_5):]
                    part.shared_edge_5 = edge_5
                    part.shared_edge_5_system = system
                    part.shared_edge_note = 'Leading ' + edge_5 + ' edge repeat is shared with the previous TasH unit and is not duplicated in the assembled array.'
                    break
                elif system == 'TasA' and previous_edge_3 == edge_3 and sequence.upper().startswith(edge_5.upper()) and edge_5.upper() != edge_3.upper():
                    part.shared_edge_3 = edge_3
                    part.shared_edge_5 = edge_5
                    part.shared_edge_5_system = system
                    part.shared_edge_junction = edge_3 + edge_5
                    part.shared_edge_note = 'Upstream ' + edge_3 + ' right edge and downstream ' + edge_5 + ' left edge form the shared TasA junction ' + part.shared_edge_junction + '.'
                    break
            previous_edge_3 = None
            for system,config in TAS_SYSTEMS.items():
                edge_3 = config['edge_3'].replace('U', 'T')
                if sequence.upper().endswith(edge_3.upper()):
                    previous_edge_3 = edge_3
                    break
            tigRNA_parts.append(part)

        return tigRNA_parts
        

def annotatePrimers(polycistron, oligo_prefix='o', oligo_index='0', staticBorderPrimers=False, noBorderPrimers=False, poltype='ptg', enzm='bsai', bb_linkers=['tgcc','gttt'], ad_linkers=[]):
    '''
    Annotates the provided primers in the polycistron

    Reusable/default border oligos are applied only in ``ptg`` mode. ``ca`` and
    ``tigRNA`` assemblies do not use the tRNA-processed PTG architecture, so
    their border oligos are target-specific and always included in the final
    oligo list.
    
    :param polycistron: Polycistron object, to which the primers should be annotated
    :type polycistron: Polycistron object
    :param oligo_prefix: Prefix to use for the oligos
    :type oligo_prefix: str
    :param oligo_index: Index from where to start the oligo numbering
    :type oligo_index: str
    :param staticBorderPrimers: Whether the border primers should be optimised or the static default ones for each architecture and restriction enzyme should be assigned
    :type staticBorderPrimers: Boolean
    :param noBorderPrimers: Whether the border primers should be omitted from the numbering and instead assigned a default name
    :type noBorderPrimers: Boolean
    :param poltype_opt: The type of polycistronic architecture to use. Must be one of 'ptg', 'ca' or 'tigRNA', defaults to 'ptg'
    :type poltype_opt: str, optional
    :param enzm: Type II restriction enzyme to use for the Golden Gate assembly. Defaults to 'bsai'
    :type enzm: str, optional
    :param bb_linkers: Linkers of the destination plasmid flanking the final PTG, defaults to ['tgcc','gttt']
    :type bb_linkers: list, optional
    :param ad_linkers: Additional linkers in the destination plasmid, defaults to []
    :type ad_linkers: list, optional
    '''
    
    ## Setting up variables
    enzms={'bsai': ['gaggtctcg', 'cgagacctc'], 'bsmbi': ['tgcgtctca', 'tgagacgca'], 'btgzi': ['ctgcgatggagtatgtta', 'taacatactccatcgcag'], 'bbsi': ['ttgaagactt', 'aagtcttcaa']} #templates found in pUU080 (bsai), pUPD2 (bsmbi), Ortega-Escalante et al. 2018 (btgzi), Kun (bbsi)
    positions = len(oligo_index)
    collapsed_index = int(oligo_index)
    if poltype != 'ptg':
        staticBorderPrimers = False
        noBorderPrimers = False
    
    
    ## Catching errors
    if oligo_index != '' and re.search(r'^[0-9]*$', oligo_index) is None:
        raise InvalidUsage("Starting index must be a number in string format", status_code=400, payload={'pge': 'sequence.html', 'box': 'oligo_index'})
    for lnk in bb_linkers+ad_linkers:
        if len(lnk) != 4 or re.search(r'^[ACGTacgt]*$', lnk) is None:
            raise InvalidUsage("Invalid linker input", status_code=400, payload={'pge': 'sequence.html', 'box': 'link'})
    
    
    ## If desired, assign the default border primer sequences
    if staticBorderPrimers or noBorderPrimers:
        if poltype == 'ptg':
            polycistron.parts[0].primer_forward = enzms[enzm][0] + reverse_complement(bb_linkers[0].upper()) + 'aacaaagcaccagtggtctagtggtag'
            polycistron.parts[-1].primer_reverse = reverse_complement(enzms[enzm][1]) + reverse_complement(bb_linkers[1].upper()) + 'tgcaccagccgggaatcgaac'
        if poltype == 'ca':
            polycistron.parts[0].primer_forward = enzms[enzm][0] + reverse_complement(bb_linkers[0].upper()) + 'aatttctactgttgtagat'
            polycistron.parts[-1].primer_reverse = reverse_complement(enzms[enzm][1]) + reverse_complement(bb_linkers[1].upper()) + 'atctacaacagtagaaatt'
    
    
    ## Write a complete list of oligo IDs and sequences to the oligos property of the polycistron object
    primerList = flattn([[i.primer_forward, i.primer_reverse] for i in polycistron.parts])
    default_fw_name = 'default_' + poltype + '_' + enzm + '_' + bb_linkers[0].lower() + '_fw'
    default_rv_name = 'default_' + poltype + '_' + enzm + '_' + bb_linkers[-1].lower() + '_rv'
    for c,primer in enumerate(primerList):
        if noBorderPrimers and c == 0:
            polycistron.oligos.append([default_fw_name, primer])
        elif noBorderPrimers and c == len(primerList) - 1:
            polycistron.oligos.append([default_rv_name, primer])
        else:
            polycistron.oligos.append([oligo_prefix+format(collapsed_index,'0' + str(positions)), primer])
            collapsed_index += 1
    
    
    ## Annotate only the oligo sequence that matches the final GenBank construct.
    ## Type IIS recognition sites are left unannotated; 4 bp overhangs are retained.
    for c,part in enumerate(polycistron.parts):
        for primer,oligo,strand,enzm_seq in [(part.primer_forward, polycistron.oligos[2*c][0], 1, enzms[enzm][0]),
                                             (part.primer_reverse, polycistron.oligos[2*c+1][0], -1, enzms[enzm][1])]:
            match = primer_construct_match(primer, polycistron.sequence, enzm_seq, strand)
            if match is not None:
                strt,end = match
                direction = 'forward' if strand == 1 else 'reverse'
                add_feature(polycistron, strt, end, 'primer_bind',
                            label=oligo,
                            note=oligo + ' ' + direction + ' oligo region matching the final construct; Type IIS recognition site is excluded and the 4 bp Golden Gate overhang is included',
                            strand=strand)
    
    return polycistron
