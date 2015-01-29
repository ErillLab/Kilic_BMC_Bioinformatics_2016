import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import re
import glob
import transfer
import motif
import pickle
from utils import list as listutils
from utils import bio as bioutils
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

RUNDIR = '/home/sefa/Desktop/runs'

def read_motif_csv():
    """Read the data into a Pandas data frame"""
    df = pd.read_csv('../data/merged_data.tsv', sep='\t')
    return df

def group_sites(df):
    """Group binding sites by TF, TF accession and genome accession"""
    return df.groupby(['TF', 'TF_accession', 'genome_accession'])

def data_by_database(df):
    """Plot the number of sites by database."""
    gs = df.groupby('database').size()
    gs.plot(kind='bar')
    plt.xlabel("Databases")
    plt.ylabel("Number of binding sites with experimental evidence")
    plt.xticks(rotation=45)
    plt.show()

def motifs_by_size(df):
    """The number of motifs by size"""
    groups = group_sites(df)
    sizes = pd.DataFrame(groups.size(), columns=['size'])
    sizes = sizes[sizes['size']<=80]
    sizes.plot(kind='hist', bins=20)
    plt.show()

def find_pairs(df, sz_thr=10):
    """Given data frame and the motif size threshold, return all pairs of
    TF/species pairs for transfer. Each pair has the same TF and different
    species where the transfer is performed from one to the other."""
    groups = group_sites(group_sites(df).filter(lambda grp: len(grp) >= sz_thr))
    pairs = [(ga, gb) for ga, _ in groups for gb, _ in groups
             if (ga[0] == gb[0] and # TF
                 ga[1] != gb[1] and # TF accession
                 ga[2] != gb[2] and # genome accession
                 not ga[0].lower().startswith('sig'))]
    return pairs

def tfs_in_pairs(pairs):
    """Find all TFs in the list of pairs"""
    return set(x[0] for x, _ in pairs)

def num_tests(df):
    """How many tests are possible?"""
    size_range = range(2, 50)
    num_transfers = pd.DataFrame(
        {'threshold': size_range,
         'num_pairs': [len(find_pairs(df, th)) for th in size_range]})
    num_transfers.plot(x='threshold', y='num_pairs', legend=False)
    plt.show()

def strip_version(accession):
    """Remove the version and return the stripped accession number (without
    version)"""
    return accession.split('.')[0]
    
def check_ortholog_files():
    """Check if orthologs are available."""
    df = read_motif_csv()
    ortholog_dir = '../data/orthologs'
    pairs = [(strip_version(ga), strip_version(gb))
             for (_, _, ga), (_, _, gb) in find_pairs(df)]
    for ga, gb in tqdm(pairs):
        alt1 = os.path.join(ortholog_dir, '%s_%s.tsv' % (ga, gb))
        alt2 = os.path.join(ortholog_dir, '%s_%s.tsv' % (gb, ga))
        assert os.path.isfile(alt1) or os.path.isfile(alt2)
        if not os.path.isfile(alt1) and not os.path.isfile(alt2):
            reciprocal_blast.find_orthologs([ga, gb])

def check_operon_files():
    """Check if operons for all genomes are available"""
    df = read_motif_csv()
    operon_dir = '../data/operons'
    orgs = [strip_version(ga) for _, (_, _, ga) in find_pairs(df)]
    for org in orgs:
        assert os.path.isfile(os.path.join(operon_dir, org+'.opr')), org

def get_protein_accessions():
    """Return dictionary {(genome_accession, TF): protein_accession}."""
    all_pairs = find_pairs(read_motif_csv(), 10)
    xs = listutils.nub([p[0] for p in all_pairs] + [p[1] for p in all_pairs])
    return dict(((x[2], x[0]), x[1]) for x in xs)

def protein_distance(acca, accb):
    """Given two protein accessions, return their sequence-based evolutionary
    distance"""
    matrix = matlist.blosum62
    gap_open = -10
    gap_extend = -1
    reca = bioutils.get_protein(acca)
    recb = bioutils.get_protein(accb)
    aln, = pairwise2.align.globalds(str(reca.seq), str(recb.seq), matrix,
                                    gap_open, gap_extend,
                                    one_alignment_only=True)
    ala, alb, score, begin, end = aln
    num_gaps = sum([1.0 for (ba, bb) in zip(ala, alb) if ba == '-' or bb == '-'])
    num_id = sum([1.0 for (ba, bb) in zip(ala, alb) if ba == bb])
    return 1 - (num_id / (len(ala) - num_gaps))

def run_transfer(method):
    df = read_motif_csv()
    pairs = find_pairs(df)
    rundir = os.path.join(RUNDIR, method)
    for i, ((tfa, pa, ga), (tfb, pb, gb)) in enumerate(pairs):
        print "Running %d of %d." % (i, len(pairs))
        print method, tfa, ga, gb
        pf = os.path.join(rundir, '_'.join([tfa, ga, gb]) + '.pkl')
        if os.path.isfile(pf):
            continue
        config = transfer.new_config(ga, gb, tfa)
        ref = transfer.get_reference(config)
        target = transfer.get_target(config)
        if method == 'direct_transfer':
            motif = transfer.direct_transfer(ref, target, config)
        elif method == 'meme_on_pssm_searched':
            motif = transfer.meme_on_pssm_searched(ref, target, config, False)
        elif method == 'meme_on_pssm_searched_seeded':
            motif = transfer.meme_on_pssm_searched(ref, target, config, True)
        elif method == 'network_transfer':
            motif = transfer.network_transfer(ref, target, config, False)
        elif method == 'network_transfer_seeded':
            motif = transfer.network_transfer(ref, target, config, True)
        elif method == 'random':
            motif = transfer.random_motif(ref, target, config)
        elif method == 'permuted':
            from motif import permute
            target_motif = transfer.get_true_motif(target, config)
            motif = permute(target_motif)
        pickle.dump(motif, open(pf, 'w'))

def post_process(method):
    """Post-processing of the inferred motifs"""
    rundir = os.path.join(RUNDIR, method)
    pkls = glob.glob(os.path.join(rundir, '*pkl'))
    protein_accessions = get_protein_accessions()
    df_cols = ['tf',
               'reference_tf_accession',
               'target_tf_accession',
               'reference_genome_accession',
               'target_genome_accession',
               'reference_motif_ic',
               'target_motif_ic',
               'inferred_motif_ic',
               'reference_motif_size',
               'target_motif_size',
               'inferred_motif_size',
               'protein_distance',
               'transfer_method',
               'reference_target_dist',
               'target_inferred_dist']
    df_dict = dict((col, [None]*len(pkls)) for col in df_cols)
    for i, pkl in enumerate(pkls):        
        match = re.match(r'\S+/(\w+)_(NC_\d+\.\d)_(NC_\d+\.\d).pkl', pkl)
        tf, ga, gb = match.groups()
        print method, "%d of %d: %s %s %s" % (i, len(pkls), tf, ga, gb)
        config = transfer.new_config(ga, gb, tf)
        ref = transfer.get_reference(config)
        target = transfer.get_target(config)
        reference_motif = transfer.get_true_motif(ref, config)
        target_motif = transfer.get_true_motif(target, config)
        inferred_motif = pickle.load(open(pkl))
        df_dict['tf'][i] = tf
        df_dict['reference_tf_accession'][i] = protein_accessions[(ga, tf)]
        df_dict['target_tf_accession'][i] = protein_accessions[(gb, tf)]
        df_dict['reference_genome_accession'][i] = ga
        df_dict['target_genome_accession'][i] = gb
        df_dict['reference_motif_ic'][i] = motif.IC(reference_motif)
        df_dict['target_motif_ic'][i] = motif.IC(target_motif)
        df_dict['inferred_motif_ic'][i] = motif.IC(inferred_motif)
        df_dict['reference_motif_size'][i] = motif.size(reference_motif)
        df_dict['target_motif_size'][i] = motif.size(target_motif)
        df_dict['inferred_motif_size'][i] = motif.size(inferred_motif)
        df_dict['protein_distance'][i] = protein_distance(
            protein_accessions[(ga, tf)], protein_accessions[(gb, tf)])
        df_dict['transfer_method'][i] = method
        df_dict['reference_target_dist'][i] \
            = motif.distance(reference_motif, target_motif)
        df_dict['target_inferred_dist'][i] \
            = motif.distance(target_motif, inferred_motif)
    df = pd.DataFrame(df_dict)
    df.to_csv('../results_%s.csv' % method, index=False, cols=df_cols)
