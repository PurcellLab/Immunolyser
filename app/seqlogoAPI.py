import sys
import os
from subprocess import call

# === Log all arguments ===
print("Arguments received:", file=sys.stderr)
for i, arg in enumerate(sys.argv):
    print(f"argv[{i}] = {arg}", file=sys.stderr)

inp = sys.argv[1]
out = sys.argv[2]
name = sys.argv[3]
nine_mers = sys.argv[4]
total_peptides = sys.argv[5]

print(f"Input file: {inp}", file=sys.stderr)
print(f"Output folder: {out}", file=sys.stderr)
print(f"Name: {name}", file=sys.stderr)
print(f"Nine mers count: {nine_mers}", file=sys.stderr)
print(f"Total peptides: {total_peptides}", file=sys.stderr)

# Read motif length from the same folder as input file
motif_length_file = os.path.join(os.path.dirname(inp), "motif_length.txt")
print(f"Reading motif length from: {motif_length_file}", file=sys.stderr)

try:
    with open(motif_length_file, 'r') as f:
        motif_length = f.read().strip()
        print(f"Motif length read: {motif_length}", file=sys.stderr)
except Exception as e:
    print(f"Error reading motif length: {e}", file=sys.stderr)
    motif_length = "UNKNOWN"

# Log final call command
cmd = [
    'python2', 'app/tools/seq2logo-2.1/Seq2Logo.py',
    '-f', inp,
    '-o', out,
    '--format', '[JPEG]',
    '-t', '{} based on {} {}-mers'.format(name, nine_mers, motif_length),
    '-S', '2',
    '-I', '2'
]

print("Calling Seq2Logo with:", file=sys.stderr)
print(" ".join(cmd), file=sys.stderr)

# Execute
call(cmd)
print("Seq2Logo call finished", file=sys.stderr)
