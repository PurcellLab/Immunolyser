import sys
import os
from subprocess import call

task_id = sys.argv[1]
data_mount = sys.argv[2]

sample = sys.argv[3]
replicate = sys.argv[4]
motif_length = int(sys.argv[5])
mhc_class = sys.argv[6]
       
replicate += '_9mer.txt'

# Argument to specify the motif length
arg_l = f'-l {motif_length}'

# Set additional parameters based on mhc_class
if mhc_class == "I":
    arg_C = '-C'
    arg_r = ''
    arg_s = ''
else:
    arg_C = ''
    arg_r = '-r 20'
    arg_s = '-s 100'

print(os.popen('perl ./app/tools/gibbscluster-2.0/GibbsCluster-2.0e_SA_for_bar.pl \
        -f {}/{}/{}/{} \
        -H R \
        -G ./../../../../../../../../app/tools/seq2logo-2.1/Seq2Logo.py \
        -g 1-5 \
        -k {} \
        -T \
        {} \
        {} \
        {} \
        -R {}/app/static/images/{}/{}/gibbscluster/{} '.format(data_mount,task_id, sample, replicate, os.cpu_count(), arg_C, arg_r, arg_s, arg_l, os.getcwd(), task_id, sample, replicate[:-9])).read())
