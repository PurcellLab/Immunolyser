import sys
import os

# Get command-line arguments
task_id = sys.argv[1]
data_mount = sys.argv[2]
sample = sys.argv[3]
replicate = sys.argv[4] + '_9mer.txt'
cluster = sys.argv[5]
mhc_class = sys.argv[6]
motif_length = int(sys.argv[7])

# Common arguments
arg_l = f'-l {motif_length}'  # Motif length
arg_g = f'-g {cluster}'  # Number of clusters
arg_k = f'-k {os.cpu_count()}'  # Number of processes
arg_T = '-T'  # Use trash cluster
arg_R = f'-R {os.getcwd()}/app/static/images/{task_id}/{sample}/gibbscluster/{replicate[:-9]}'

# Set parameters based on MHC class
if mhc_class == "I":
    # Settings for MHC class I
    arg_C = '-C'  # Single sequence moves every iteration
    additional_args = "-I 0 -D 0 -S 5 -b 0.8 -q 5 -c 0 -z 1 -j 2"
else:
    # Settings for MHC class II
    arg_C = ''  # No -C, meaning shift moves are enabled
    additional_args = "-I 0 -D 0 -S 5 -b 0.8 -q 5 -c 0 -z 1 -j 2 -s 100 -r 20"

# Construct the command
cmd = (
    f'perl ./app/tools/gibbscluster-2.0/GibbsCluster-2.0e_SA_for_seqlogo.pl '
    f'-f {data_mount}/{task_id}/{sample}/{replicate} '
    f'-H R '
    f'-G ./../../../../../../../../app/tools/seq2logo-2.1/Seq2Logo.py '
    f'{arg_g} {arg_k} {arg_T} {arg_C} {arg_l} {additional_args} {arg_R}'
)

# Print and execute the command
print("Gibbs Command:", cmd)
print(os.popen(cmd).read())
