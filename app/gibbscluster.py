import sys
import os
from subprocess import call

# Get command-line arguments
task_id = sys.argv[1]
data_mount = sys.argv[2]
mhc_class = sys.argv[3]
motif_length = int(sys.argv[4])

# Determine the project root
project_root = os.path.dirname(os.path.realpath(os.path.join(__file__, "..")))

# Set seq2logo path
seq2logo_path = os.path.join(project_root, "app/tools/seq2logo-2.1/Seq2Logo.py")

task_path = os.path.join(data_mount, task_id)

# Default settings for both classes
num_clusters = '1-6'
mc_temperature = '1.5'
num_temp_steps = '20'
num_iterations = '10'
num_seeds = '5'
penalty_lambda = '0.8'
small_cluster_weight = '5'
sequence_weighting_type = '0'
background_model = '1'  # Pre-calculated Uniprot
use_trash_cluster = '1'
trash_cluster_threshold = '2'

# Command for MHC class I
if mhc_class == "I":
    max_insertion_length = '1'
    max_deletion_length = '4'
    indel_move_interval = '10'
    shift_move_interval = '20'  # Class I specific
    phase_shift_move_interval = '100'  # Class I specific
    hydrophobic_p1_preference = '0'

    input_file_ends_with = '_8to14mer.txt'

    # Command to run GibbsCluster for MHC Class I
    for sample in os.listdir(task_path):
        sample_path = os.path.join(task_path, sample)
        if not os.path.isdir(sample_path):
            continue
        for replicate in os.listdir(sample_path):
            if replicate.endswith(input_file_ends_with):
                print('Replicate file:', replicate)
                gibbs_command = (
                    f'perl {project_root}/app/tools/gibbscluster-2.0/GibbsCluster-2.0e_SA.pl '
                    f'-f {sample_path}/{replicate} '
                    f'-H R '
                    f'-G {seq2logo_path} '
                    f'-g {num_clusters} '  # 1-5 clusters
                    f'-k {os.cpu_count()} '
                    f'-T -C -l {motif_length} '
                    f'-R {project_root}/app/static/images/{task_id}/{sample}/gibbscluster/{replicate[:-13]} '
                    f'-i {num_iterations} '
                    f'-t {mc_temperature} '
                    f'-n {num_temp_steps} '
                    f'-b {penalty_lambda} '
                    f'-q {small_cluster_weight} '
                    f'-S {num_seeds} '
                    f'-c {sequence_weighting_type} '
                    f'-z {background_model} '
                    f'-T {use_trash_cluster} '
                    f'-j {trash_cluster_threshold} '
                    f'-D {max_deletion_length} '
                    f'-I {max_insertion_length} '
                    f'-u {indel_move_interval} '
                    f'-x {shift_move_interval} '
                    f'-s {phase_shift_move_interval} '
                    f'-p {hydrophobic_p1_preference} '
                )
                
                print("Gibbs Command for MHC Class I:", gibbs_command)
                print(os.popen(gibbs_command).read())

# Command for MHC class II
elif mhc_class == "II":
    max_insertion_length = '0'  # No insertion for Class II
    max_deletion_length = '0'   # No deletion for Class II
    indel_move_interval = '20'  # For Class II, different move interval
    shift_move_interval = '20'  # Shift moves are activated
    phase_shift_move_interval = '100'  # Phase shift interval for Class II
    hydrophobic_p1_preference = '1'  # Class II specific (could be 1 or a different preference)
    arg_C = ''

    input_file_ends_with = '_12to20mer.txt'

    # Command to run GibbsCluster for MHC Class II
    for sample in os.listdir(task_path):
        sample_path = os.path.join(task_path, sample)
        if not os.path.isdir(sample_path):
            continue
        for replicate in os.listdir(sample_path):
            if replicate.endswith(input_file_ends_with):
                print('Replicate file:', replicate)
                gibbs_command = (
                    f'perl {project_root}/app/tools/gibbscluster-2.0/GibbsCluster-2.0e_SA.pl '
                    f'-f {sample_path}/{replicate} '
                    f'-H R '
                    f'-G {seq2logo_path} '
                    f'-g {num_clusters} '  # 1-5 clusters
                    f'-k {os.cpu_count()} '
                    f'-T {arg_C} -l {motif_length} '  # No -C for Class II
                    f'-R {project_root}/app/static/images/{task_id}/{sample}/gibbscluster/{replicate[:-14]} '
                    f'-i {num_iterations} '
                    f'-t {mc_temperature} '
                    f'-n {num_temp_steps} '
                    f'-b {penalty_lambda} '
                    f'-q {small_cluster_weight} '
                    f'-S {num_seeds} '
                    f'-c {sequence_weighting_type} '
                    f'-z {background_model} '
                    f'-T {use_trash_cluster} '
                    f'-j {trash_cluster_threshold} '
                    f'-D {max_deletion_length} '
                    f'-I {max_insertion_length} '
                    f'-u {indel_move_interval} '
                    f'-x {shift_move_interval} '
                    f'-s {phase_shift_move_interval} '
                    f'-p {hydrophobic_p1_preference} '
                )
                
                print("Gibbs Command for MHC Class II:", gibbs_command)
                print(os.popen(gibbs_command).read())
