from collections import namedtuple
from utils import list as listutils
from utils import bio as bioutils
import sequence
import genome
from Bio import motifs
import yassi
import math
import random
from tqdm import tqdm

class Motif(namedtuple('Motif', 'sites motif_')):
    """Motif class derived from named-tuple. Inherited from named-tuple to
    override '__repr__' function"""
    def __repr__(self):
        return '\n'.join([str(site) for site in self.sites])

def new_motif(sites):
    """Given sites, return motif object"""
    sites = listutils.nub_by(sequence.overlap_test, sites)
    seqs = [site.seq for site in sites]
    motif_ = motifs.Motif(instances=motifs.Instances(seqs))
    motif_.pseudocounts = dict(A=0.25, C=0.25, G=0.25, T=0.25)
    return Motif(sites, motif_)

def reverse_complement(motif):
    """Reverse complement of a motif"""
    sites = [sequence.reverse_complement(site) for site in motif.sites]
    return new_motif(sites)

def permute(motif):
    """Permute the given motif by shuffling its columns"""
    cols = range(length(motif))
    random.shuffle(cols)
    shuffled = [''.join(site[i] for i in cols) for site in seqs(motif)]
    _motif = motifs.Motif(instances=motifs.Instances(shuffled))
    _motif.pseudocounts = pseudocounts(motif)
    return Motif(None, _motif)

def length(motif):
    """Motif length"""
    return motif.motif_.length

def pseudocounts(motif):
    """Return used pseudocounts"""
    return motif.motif_.pseudocounts

def size(motif):
    """Number of sites in the motif"""
    return len(motif.motif_.instances)

def seqs(motif):
    """Return sites as string"""
    return [str(inst) for inst in motif.motif_.instances]

def regulon(motif):
    """Return the list of regulated genes"""
    return listutils.nub([sequence.nearby_operon(site) for site in motif.sites])

def get_genome(motif):
    """Return the genome of the motif"""
    assert len(listutils.nub([site.genome for site in motif.sites])) == 1
    return motif.sites[0].genome

def pssm(motif):
    """Position-specific scoring matrix of the motif"""
    return motif.motif_.pssm

def IC(motif):
    """Information content of the motif"""
    return pssm(motif).mean()

def pssm_col(motif, col):
    """A column of the PSSM"""
    return dict((let, pssm(motif)[let][col]) for let in "ACTG")

def pwm(motif):
    """Position weight matrix of the motif"""
    return motif.motif_.pwm

def pwm_col(motif, col):
    """A column of the PWM"""
    return dict((let, pwm(motif)[let][col]) for let in "ACTG")

def background(motif):
    """Background dist. of the motif"""
    return [motif.motif_.background[let] for let in "ACTG"]

def consensus(motif):
    """Consensus sequence"""
    return motif.motif_.consensus

def weblogo(motif):
    return bioutils.weblogo(motif.motif_)

def pssm_search(motif, seq):
    results = yassi.search(seqs(motif), seq, background(motif))
    mlen = length(motif)
    scores = dict((seq[pos:pos+mlen], scr) for pos, scr in results)
    return scores

def pssm_search_on_regions(motif, regions, nsites):
    """Find high-scoring nsites"""
    seq = ''.join(region.seq for region in regions)
    scores = pssm_search(motif, seq + bioutils.reverse_complement(seq))
    mlen = length(motif)
    sites = []
    threshold = -51
    for region in tqdm(regions):
        rseq = region.seq if region.strand == 1 else bioutils.reverse_complement(region.seq)
        rseq_comp = bioutils.complement(rseq)
        for offset in xrange(len(rseq)-mlen+1):
            site = rseq[offset:offset+mlen]
            rsite = rseq_comp[offset:offset+mlen][::-1]
            if max(scores[site], scores[rsite]) > threshold:
                newsite = genome.subseq(region.genome,
                                        region.start+offset,
                                        region.start+offset+mlen,
                                        1 if scores[site] > scores[rsite] else -1)
                if not any(sequence.overlap_test(newsite, site) for site in sites):
                    sites.append(newsite)
            if len(sites) > nsites:
                sites.sort(key=lambda site: scores[site.seq], reverse=True)
                sites = sites[:nsites]
                threshold = scores[sites[-1].seq]
    assert len(sites) == nsites
    return sites

def alignment(motif, other):
    """Align two motifs that give the maximum IC of the aligned region"""
    max_ic = float('-inf')
    for offset in range(-length(motif) + 1, length(other)):
        if offset < 0:
            ic = ic_at(motif, other, -offset)
        else:
            ic = ic_at(other, motif, offset)

        
        if ic > max_ic:
            max_ic = ic
            max_offset = offset
    return max_offset

def ic_at(motif, other, offset):
    """Return the total IC of two aligned motifs"""
    alignment_len = min(length(motif)-offset, length(other))
    motif_seqs = [site[offset:alignment_len+offset] for site in seqs(motif)]
    other_seqs = [site[:alignment_len] for site in seqs(other)]
    # Create the motif and compute the IC
    amotif = motifs.Motif(instances=motifs.Instances(motif_seqs+other_seqs))
    amotif.pseudocounts = dict(A=0.25, C=0.25, G=0.25, T=0.25)
    return amotif.pssm.mean()


def euclidean(cola, colb):
    """Euclidean distance between two frequency vector"""
    return math.sqrt(sum((cola[let]-colb[let])**2 for let in "ACTG"))

def kl(cola, colb):
    """KL divergence between two frequency vector"""
    safe_log2 = lambda x: math.log(x, 2) if x != 0 else 0.0
    return (sum(cola[l] * safe_log2(cola[l] / colb[l]) for l in "ACTG") +
            sum(colb[l] * safe_log2(colb[l] / cola[l]) for l in "ACTG"))

def distance(motif, other, col_dist_func=euclidean, offset=None):
    """Compute the distance between two motifs, using the given column distance
    function."""
    if offset is None:
        offset = alignment(motif, other)
    if offset < 0:
        return distance(other, motif, col_dist_func, -offset)
    dists = []
    for pos in range(min(length(motif), length(other)-offset)):
        cola = pwm_col(motif, pos)
        colb = pwm_col(other, pos+offset)
        dists.append(col_dist_func(cola, colb))
    return sum(dists) / len(dists)
