from collections import namedtuple, defaultdict
from utils import bio as bioutils
import os
import pandas as pd
import gene
import sequence
import random

class Genome(namedtuple('Genome', 'record sequence operons')):
    def __repr__(self):
        return ' '.join(self.record.description.split()[:2])

def new_genome(gacc, config):
    """Create a new genome"""
    record = bioutils.get_genome(gacc, config.genome_dir)
    genes = [gene.new_gene(f) for f in record.features if f.type == 'gene']
    operons = group_genes(gacc, genes, config.operon_dir)
    return Genome(record, str(record.seq), operons)

def subseq(genome, start, end, strand):
    """Return a subsequence from genome"""
    return sequence.new_seq(genome, start, end, strand)

def random_seq(genome, seq_len):
    """Return a randomly generated sequence from the genome."""
    start_pos = random.randint(0, length(genome)-seq_len-1)
    return subseq(genome, start_pos, start_pos+seq_len, random.choice([1, -1]))

def accession(genome):
    """Return genome accession number without the version"""
    return genome.record.id.split('.')[0]

def accession_full(genome):
    """Return genome accession number"""
    return genome.record.id

def description(genome):
    """Return genome description"""
    return genome.record.description

def length(genome):
    """Genome length"""
    return len(genome.sequence)

def read_operons(operon_file):
    """Read operons from given file. Return a dictionary of operons {locus_tag:
    operon_id}
    """
    df = pd.read_csv(operon_file, sep='\t')
    operons = dict(zip(df['Synonym'], df['OperonID']))
    return operons

def group_genes(gacc, genes, operon_dir):
    """Group list of genes by operon"""
    # read operon file
    opr_file = os.path.join(operon_dir, gacc.split('.')[0] + '.opr')
    opr_dict = read_operons(opr_file)
    # group genes into operons
    operons = defaultdict(list)
    not_found = -1              # use as operon id for genes not in opr_dict
    for g in genes:
        opr_id = opr_dict.get(g.locus_tag, None)
        if opr_id:
            operons[opr_id].append(g)
        else:
            operons[not_found] = [g]
            not_found -= 1      # new id for the next not-found
    # sort all operons
    for gs in operons.values():
        gs.sort(key=lambda g: g.start, reverse=gs[0].strand == -1)
    return operons.values()

def promoters(genome):
    """Return the list of upstream regions of all genes."""
    proms = sequence.merge_overlapping_seqs(
        [get_operon_upstream(genome, opr) for opr in genome.operons])
    return proms

def get_operon_upstream(genome, operon):
    """Return the promoter region of the operon"""
    operon.sort(key=lambda gene: gene.start)
    if operon[0].strand == 1:
        return gene.upstream_region(genome, operon[0])
    else:
        return gene.upstream_region(genome, operon[-1])

def orthologs(genome, other, ortholog_dir):
    def ortholog_file():
        """Determine the file name containing orthologs between two orgs"""
        template = os.path.join(ortholog_dir, '%s_%s.tsv')
        alt1 = template % (accession(genome), accession(other))
        alt2 = template % (accession(other), accession(genome))
        return alt1 if os.path.isfile(alt1) else alt2
    print ortholog_file()
    df = pd.read_csv(ortholog_file(), sep='\t')
    lt_orthologs = dict(zip(df['org1_locus_tag'], df['org2_locus_tag']) +
                        zip(df['org2_locus_tag'], df['org1_locus_tag']))
    orthologs = {}
    for opr in genome.operons:
        for g in opr:
            og = lt_orthologs.get(g.locus_tag, None)
            if og:
                orthologs[g] = gene_by_locus_tag(other, og)
    return orthologs

def gene_by_locus_tag(genome, locus_tag):
    """Get gene by locus tag"""
    try:
        return [g for opr in genome.operons for g in opr
                if g.locus_tag == locus_tag][0]
    except:
        print "Invalid locus tag", locus_tag
        return None
