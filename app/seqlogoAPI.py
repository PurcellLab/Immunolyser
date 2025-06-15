import sys
import os
from subprocess import call

inp = sys.argv[1]
out = sys.argv[2]
name = sys.argv[3]
nine_mers = sys.argv[4]
total_peptides = sys.argv[5]

# Read motif length from the same folder as input file
motif_length_file = os.path.join(os.path.dirname(inp), "motif_length.txt")
with open(motif_length_file, 'r') as f:
    motif_length = f.read().strip()

call([
    'python2', 'app/tools/seq2logo-2.1/Seq2Logo.py',
    '-f', inp,
    '-o', out,
    '--format', '[JPEG]',
    '-t', '{} based on {} {}-mers'.format(name, nine_mers, motif_length),
    '-S', '2',
    '-I', '2'
])
