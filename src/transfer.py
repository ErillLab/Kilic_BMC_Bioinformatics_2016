
from collections import namedtuple
import genome
import sequence
import motif
import meme
import gene
from utils import list as listutils
from utils import bio as bioutils
import pandas as pd
import lasagna
import random

Config = namedtuple('Config', ['reference',
                               'target',
                               'TF',
                               'meme_path',
                               'meme_mod',
                               'meme_nmotifs',
                               'data_dir',
                               'genome_dir',
                               'operon_dir',
                               'ortholog_dir',
                               'motif_file'])

def sample_config():
    return Config(reference='NC_011333.1',
                  target='NC_002505.1',
                  TF='Fur',
                  meme_path='/home/sefa/meme/bin/meme',
                  meme_mod='zoops', # zero or one per sequence
                  meme_nmotifs=1,   # number of motifs
                  data_dir='/home/sefa/Dropbox/cg/data',
                  genome_dir='/home/sefa/Dropbox/cg/data/genomes',
                  operon_dir='/home/sefa/Dropbox/cg/data/operons',
                  ortholog_dir='/home/sefa/Dropbox/cg/data/orthologs',
                  motif_file='/home/sefa/Dropbox/cg/data/merged_data.tsv')

def new_config(reference, target, TF):
    return Config(reference=reference,
                  target=target,
                  TF=TF,
                  meme_path='/home/sefa/meme/bin/meme',
                  meme_mod='zoops', # zero or one per sequence
                  meme_nmotifs=1,   # number of motifs
                  data_dir='/home/sefa/Dropbox/cg/data',
                  genome_dir='/home/sefa/Dropbox/cg/data/genomes',
                  operon_dir='/home/sefa/Dropbox/cg/data/operons',
                  ortholog_dir='/home/sefa/Dropbox/cg/data/orthologs',
                  motif_file='/home/sefa/Dropbox/cg/data/merged_data.tsv')

def get_reference(config):
    return genome.new_genome(config.reference, config)

def get_target(config):
    return genome.new_genome(config.target, config)

def get_true_motif(genome, config):
    """Return the motif for the given genome and TF specified in config"""
    return get_motif_from_tsv(genome, config)

def get_motif_from_tsv(genom, config):
    """Create a motif by reading csv file"""
    df = pd.read_csv(config.motif_file, sep='\t')
    df = df[(df.TF == config.TF) &
            (df.genome_accession == genome.accession_full(genom))]
    sites = run_lasagna(df, genom)
    return motif.new_motif(sites)

def run_lasagna(df, genom):
    """Given a data frame containing binding sites for a given TF and species,
    if the binding sites are unaligned variable-length sequences, run LASAGNA
    (doi:10.1186/1471-2105-14-108) to align them. Return a data frame with
    aligned sites and their locations on the genome"""
    def fill_gaps(row, aligned_seq, aligned_strand):
        """Given a row of binding site data frame, its aligned sequence and
        strand, return the binding site with gaps filled"""
        is_gap = lambda x: x == '-'
        left_gaps = len(listutils.take_while(is_gap, aligned_seq))
        right_gaps = len(listutils.take_while(is_gap, aligned_seq[::-1]))
        rc = bioutils.reverse_complement
        # Handle each case separately
        if aligned_strand == '+' and row.site_strand == 1:
            sseq = ((row.left_flanking[-left_gaps:] if left_gaps > 0 else '') +
                    aligned_seq.strip('-') + row.right_flanking[:right_gaps])
            sstart = row.site_start - left_gaps
            send = row.site_end + right_gaps
            sstrand = 1
        elif aligned_strand == '-' and row.site_strand == 1:
            sseq = (rc(row.right_flanking[:left_gaps]) +
                    aligned_seq.strip('-') +
                    (rc(row.left_flanking[-right_gaps:])
                     if right_gaps > 0 else ''))
            sstart = row.site_start - right_gaps
            send = row.site_end + left_gaps
            sstrand = -1
        elif aligned_strand == '+' and row.site_strand == -1:
            sseq = ((row.left_flanking[-left_gaps:] if left_gaps > 0 else '') +
                    aligned_seq.strip('-') + row.right_flanking[:right_gaps])
            sstart = row.site_start - right_gaps
            send = row.site_end + left_gaps
            sstrand = -1
        elif aligned_strand == '-' and row.site_strand == -1:
            sseq = (rc(row.right_flanking[:left_gaps]) +
                    aligned_seq.strip('-') +
                    (rc(row.left_flanking[-right_gaps:])
                     if right_gaps > 0 else ''))
            sstart = row.site_start - left_gaps
            send = row.site_end + right_gaps
            sstrand = 1
        else:
            assert False, "impossible case"
        return sequence.new_seq(genom, sstart, send, sstrand, sseq)

    def suggest_trim(aligned, ic_thr=0.0, cov_thr=0.4):
        """Given a list of aligned sequences with gaps, return left and right
        trim lengths, supported by LASAGNA"""
        counts = lasagna.ComputeCounts(aligned, 0, alphabet='acgt')[1]
        IC = lasagna.ComputeIC(aligned, 0, counts=counts, counts_2mer=None)[0]
        coverages = lasagna.ComputeCoverage(counts, len(aligned), alphabet='acgt')
        sugl, cntl, cntr = lasagna.SuggestSize(IC, coverages, ICThres=ic_thr, covThres=cov_thr)
        return cntl, cntr

    sites = [seq for seq in df.site_sequence]
    aligned, idx_aligned, strands = lasagna.LASAGNA([s.lower() for s in sites], 0)
    aligned = [seq.upper() for seq in aligned]
    seqs_filled = [fill_gaps(df.iloc[i], aligned[i], strands[i])
                   for i in xrange(len(sites))]
    lt, rt = suggest_trim(aligned)
    return [sequence.trim(seq, lt, rt) for seq in seqs_filled]


def init():
    config = sample_config()
    ref = get_reference(config)
    target = get_target(config)
    return ref, target, config

def direct_transfer(ref, target, config):
    true_motif = get_true_motif(ref, config)
    nsites = int(1.15 * motif.size(true_motif))
    sites = motif.pssm_search_on_regions(true_motif, genome.promoters(target), nsites)
    return motif.new_motif(sites)

def get_meme_settings(config, ref, seeded):
    m = get_true_motif(ref, config)
    meme_settings = {
        'path': config.meme_path,
        'mod': config.meme_mod,
        'nmotifs': config.meme_nmotifs,
        'minw': int(motif.length(m) / 2.),
        'maxw': int(motif.length(m) * 1.5),
    }
    if seeded:
        meme_settings['cons'] = motif.consensus(m)
    return meme_settings

def meme_on_pssm_searched(ref, target, config, seeded):
    """Scan the target genome with the reference motif. Pass regions with
    putative site to the motif discovery as input sequences"""
    meme_settings = get_meme_settings(config, ref, seeded)
    true_motif = get_true_motif(ref, config)
    nsites = int(2.15 * motif.size(true_motif))
    psites = motif.pssm_search_on_regions(true_motif, genome.promoters(target), nsites)
    regions = listutils.nub([sequence.expand(site) for site in psites])
    for reg in regions:
        print reg.start, reg.end, reg.strand
    motifs = [motif.new_motif(sites)
              for sites in meme.motif_discovery(regions, meme_settings)]
    return motifs

def network_transfer(ref, target, config, seeded):
    """Given reference regulon, identify the regulon that is orthologous to the
    reference and therefore the motif. Extension to the network transfer is to
    perform motif discovery seeded."""
    meme_settings = get_meme_settings(config, ref, seeded)
    true_motif = get_true_motif(ref, config)
    true_regulon = motif.regulon(true_motif)
    # Get orthologs
    orthologs = genome.orthologs(ref, target, config.ortholog_dir)
    # Get target regulon (inferred)
    inf_regulon = listutils.nub([gene.operon(orthologs[g], target)
                                 for opr in true_regulon
                                 for g in opr
                                 if orthologs.get(g, None)])
    promoters = listutils.nub([gene.upstream_region(target, opr[0])
                               for opr in inf_regulon])
    # If less than 3 promoters add random promoters to be able to run MEME.
    while len(promoters) < 3:
        promoters.append(random.choice(genome.promoters(target)))
    motifs = [motif.new_motif(sites)
              for sites in meme.motif_discovery(promoters, meme_settings)]
    return motifs

def random_motif(ref, target, config):
    """Given reference motif M having k sites, generate k random subsequences
    from target genome and build a motif from randomly generated sites"""
    true_motif = get_true_motif(ref, config)
    mlen = motif.length(true_motif)
    msize = motif.size(true_motif)
    sites = [genome.random_seq(target, mlen) for _ in range(msize)]
    return motif.new_motif(sites)

