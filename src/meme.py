import os
import tempfile
from subprocess import call
from Bio.Motif.Parsers import MEME
import sequence

def seqs_to_fasta(seqs):
    """Given a list of sequences, return a string of sites in FASTA format."""
    return '\n'.join(">seq_%d\n%s" % (i, seq) for i, seq in enumerate(seqs))

def run_meme(seqs, outdir, meme_settings):
    """Motif discovery on given sequences via MEME. Return MEME record"""
    infile = tempfile.NamedTemporaryFile()
    with open(infile.name, 'w') as f:
        f.write(seqs_to_fasta(seqs))

    print seqs_to_fasta(seqs)

    # Call shell command
    cmd = [meme_settings['path'], infile.name, '-oc', outdir,
           '-dna',              # sequences use DNA alphabet
           '-revcomp',          # allow sites on both strands
           '-maxsize', 10000000] # max dataset size on chars

    # extend meme command
    for opt in ['cons', 'nmotifs', 'mod', 'minw', 'maxw']:
        if opt in meme_settings:
            cmd.extend(['-'+opt, meme_settings[opt]])
    cmd = map(str, cmd)

    #print "Running MEME."
    print ' '.join(cmd)
    call(cmd)
    # Parse text output of the MEME
    with open(os.path.join(outdir, 'meme.txt')) as f:
        rec = MEME.read(f)
    return rec

def motif_discovery(regions, meme_settings):
    """Run MEME on given regions."""
    outdir = tempfile.mkdtemp()
    rec = run_meme([reg.seq for reg in regions], outdir, meme_settings)
    # Create a Seq object for each site in the motif.
    motifs = [filter(lambda x: x is not None,
                     [sequence.search(reg, str(inst))
                      for inst in m.instances for reg in regions])
              for m in rec.motifs]
    return motifs
