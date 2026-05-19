import sys
import os
from subprocess import call

task_id = sys.argv[1]
data_mount = sys.argv[2]
motif_length = int(sys.argv[3])


for sample in os.listdir('{}/{}'.format(data_mount,task_id)):
    if not os.path.isdir(os.path.join(data_mount, task_id, sample)):
        continue
    for replicate in os.listdir('{}/{}/{}'.format(data_mount,task_id,sample)):

        if replicate[-9:] == f'_{motif_length}mer.txt':

            print('Seq2Logo Commands', ['python2', 'app/tools/seq2logo-2.1/Seq2Logo.py', '-f', '{}/{}/{}/{}'.format(data_mount,task_id, sample, replicate), '-o', '{}/{}/{}/seqlogos/{}'.format('app/static/images', task_id, sample, replicate[:-9]), '--format', '[JPEG]', '-t', '{}'.format(replicate[:-9]), '-S', '2', '-I', '2'])

            call(['python2', 'app/tools/seq2logo-2.1/Seq2Logo.py', '-f', '{}/{}/{}/{}'.format(data_mount,task_id, sample, replicate), '-o', '{}/{}/{}/seqlogos/{}'.format('app/static/images', task_id, sample, replicate[:-9]), '--format', '[JPEG]', '-t', '{}'.format(replicate[:-9]), '-S', '2', '-I', '2'])
