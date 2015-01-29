from collections import namedtuple
import genome

Gene = namedtuple('Gene', 'locus_tag start end strand')

def new_gene(feature):
    """Return a gene object"""
    return Gene(locus_tag=feature.qualifiers['locus_tag'][0],
                start=feature.location.start.position,
                end=feature.location.end.position,
                strand=feature.location.strand)


def upstream_region(genom, gene, up=300, down=50):
    """Upstream region of the gene"""
    if gene.strand == 1:
        start = max(0, gene.start-up)
        end = min(genome.length(genom), gene.start+down)
        return genome.subseq(genom, start, end, 1)
    else:
        start = max(0, gene.end-down)
        end = min(genome.length(genom), gene.end+up)
        return genome.subseq(genom, start, end, -1)

def operon(gene, genom):
    """Return the operon  that the given gene belongs to"""
    opr, = [opr for opr in genom.operons if gene in opr]
    return opr
