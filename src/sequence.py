from collections import namedtuple
from utils import bio as bioutils
import genome
import gene
import random

class Sequence(namedtuple('Sequence', 'genome, start, end, strand, seq')):
    def __repr__(self):
        return '%s[%d, %d] %s' % ('+' if self.strand == 1 else '-',
                                  self.start, self.end, self.seq)

def new_seq(genome, start, end, strand, seq_check=None):
    s = genome.sequence[start:end]
    if strand == -1:
        s = bioutils.reverse_complement(s)
    assert seq_check is None or s == seq_check, s + '|' + seq_check
    return Sequence(genome, start, end, strand, s)

def length(seq):
    return seq.end - seq.start

def reverse_complement(seq):
    return new_seq(seq.genome, seq.start, seq.end, -seq.strand)

def trim(seq, left, right):
    """Trim sequence from left and right"""
    if seq.strand == 1:
        return new_seq(seq.genome, seq.start+left, seq.end-right, seq.strand)
    else:
        return new_seq(seq.genome, seq.start+right, seq.end-left, seq.strand)

def overlap_test(seq, other):
    """Check if two sequences overlap enough."""
    def get_overlap(loca, locb):
        """Given two locations, return the overlap ratio."""
        overlap_len = max(0, min(loca[1], locb[1]) - max(loca[0], locb[0]))
        return float(overlap_len) / (loca[1]-loca[0]+1)
    loca = (seq.start, seq.end)
    locb = (other.start, other.end)
    overlap_a = get_overlap(loca, locb)
    overlap_b = get_overlap(locb, loca)
    return (genome.accession(seq.genome) == genome.accession(other.genome)
            and (overlap_a + overlap_b) / 2.0 >= 0.75)

def expand(seq):
    """Given a site (assuming in promoter region), expand it to entire promoter
    region"""
    proms = [prom for prom in genome.promoters(seq.genome)
             if prom.start <= seq.start and prom.end >= seq.end]
    assert len(proms) == 1
    return proms[0]

def nearby_gene(seq):
    """Return nearby gene. When the sequence is a binding site, this
    function returns the gene that is likely regulated by the TF binding to
    this site."""
    genes = [opr[0] for opr in seq.genome.operons
             if ((opr[0].strand == 1 and opr[0].end > seq.start) or
                 (opr[0].strand == -1 and opr[0].start < seq.end))]
    return min(*genes, key=lambda g: distance(g, seq))

def nearby_operon(seq):
    """Return the nearby operon (list of genes) of the given sequence"""
    return gene.operon(nearby_gene(seq), seq.genome)

def distance(seq, other):
    """Return the distance between two sequence objects"""
    assert seq.start < seq.end and other.start < other.end
    return max(seq.start, other.start) - min(seq.end, other.end)


def search(seq, other_str):
    """Search the string in the given sequence object"""
    if seq.strand == -1:
        return search(reverse_complement(seq), other_str)
    index = seq.seq.find(other_str)
    if index >= 0:
        return new_seq(seq.genome, seq.start+index,
                       seq.start+index+len(other_str), 1)
    else:
        rindex = seq.seq.find(bioutils.reverse_complement(other_str))
        if rindex >= 0:
            return new_seq(seq.genome, seq.start+rindex,
                           seq.start+rindex+len(other_str), -1)
    return None

def merge_overlapping_seqs(seqs):
    """Given a collection of sequences, merge ones that overlap and return the
    merged list of sequences.
    http://www.geeksforgeeks.org/merging-intervals/"""
    # Make sure they are all from the same genome
    seqs.sort(key=lambda s: s.start)
    stack = [seqs[0]]
    for i in xrange(1, len(seqs)):
        top = stack[-1]
        # If current interval is not overlapping with stack top, push it to stack
        if top.end < seqs[i].start:
            stack.append(seqs[i])
        # Otherwise update the ending of the top
        elif top.end < seqs[i].end:
            top = genome.subseq(top.genome, top.start, seqs[i].end, top.strand)
            stack.pop()
            stack.append(top)
    return stack
            
