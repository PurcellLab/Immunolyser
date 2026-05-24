import plotly, json, re, os, glob, shutil, subprocess
import plotly.graph_objs as go
import pandas as pd
from numpy import std, mean
from statistics import stdev
from subprocess import call, Popen, run, DEVNULL
from distutils.dir_util import copy_tree
from zipfile import ZipFile
from os.path import basename
from app import app
from constants import *
from mhcflurry import Class1PresentationPredictor
from pathlib import Path
from collections import defaultdict
from ProtPeptigram.runner import run_pipeline

project_root = os.path.dirname(os.path.realpath(os.path.join(__file__, "..")))
data_mount = app.config['IMMUNOLYSER_DATA']
regex_find_cureved_brackets = '\(.*?\)'
regex_find_square_brackets = '\[.*?\]'

def plot_lenght_distribution(samples, hist="percent", taskId=None):
    fig = go.Figure()
    export_dir = os.path.join(project_root, "app" ,"static", "images", taskId,  "export", "peptide_length_distribution")
    os.makedirs(export_dir, exist_ok=True)  # Make sure the folder exists

    for sample_name, sample in samples.items():
        peptideProportion = {}

        if hist == 'percent':
            for replicate, data in sample.items():
                peptideProportion[replicate] = data.groupby('Length').count()['Peptide'] / data.shape[0] * 100
            yaxis_label = '% Peptides'
            file_suffix = 'percentage'
        else:
            for replicate, data in sample.items():
                peptideProportion[replicate] = data.groupby('Length').count()['Peptide']
            yaxis_label = 'Number of Peptides'
            file_suffix = 'absolute'

        bardatacombined = pd.concat(peptideProportion, axis=1).apply(lambda x: mean(x), axis=1)
        bardatacombined = bardatacombined.to_frame().reset_index().rename(columns={0: 'Count'})

        # Save CSV
        csv_filename = f"{sample_name}_{file_suffix}.csv"
        csv_path = os.path.join(export_dir, csv_filename)
        bardatacombined.to_csv(csv_path, index=False)

        # Optional: calculate error bars
        if len(peptideProportion) > 1:
            # Uncomment if using error bars later
            errors = pd.concat(peptideProportion, axis=1).std(axis=1)
            fig.add_trace(go.Bar(
                x=bardatacombined['Length'],
                y=bardatacombined['Count'],
                name=sample_name,
                error_y=dict(
                    type='data',
                    array=errors,  # Optional
                    color='green',
                    thickness=1,
                    width=3,
                )
            ))
        else:
            fig.add_trace(go.Bar(
                x=bardatacombined['Length'],
                y=bardatacombined['Count'],
                name=sample_name
            ))

    fig.update_layout(
        xaxis=dict(title='<i>Length</i>'),
        yaxis=dict(title=f'<i>{yaxis_label}</i>')
    )

    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    return graphJSON
    
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
# The following method filters the data to remove contamination.
# This method is specific to a PEAKS output file for peptides.
def filterPeaksFile(samples, minLen=1, maxLen=133):
    
    for file_name,sample in samples.items():

        print('---Filtering {} file---'.format(file_name))
        print('Total number of peptides : {}'.format(sample.shape[0]))

#       Temporary variable to store the changes done because of filtering process
        temp = sample
        
#       Dropping null Peptides
        temp = temp[temp['Peptide'].apply(lambda x: isinstance(x, str) and x.strip() != '')]

#       Removing contamincation founf from accession number
        if temp.columns.__contains__('Accession'):
            temp = temp[temp.apply(lambda x : str(x['Accession']).find('CONTAM') == -1,axis=1)]
            temp = temp[temp.apply(lambda x : str(x['Accession']).find('DECOY') == -1,axis=1)]
            print('Number of peptides after removing peptides with accession marked as #CONTAM or #DECOY : {}'.format(temp.shape[0]))

#       Removing PTMs in Peptide column
        if temp.shape[0] == 0:
            raise Exception(f"After filtering for CONTAM/DECOY, no valid peptides remain in file '{file_name}'.")

        # Now safe to apply omitPTMContent
        temp['Peptide'] = temp['Peptide'].apply(omitPTMContent)

#       Generating Length Column
        temp['Length']= temp.apply(lambda x : len(x['Peptide']), axis=1)

#       Filtering on the basis of the peptide length
        temp = temp[temp.apply(lambda x : x['Length'] in range(minLen,maxLen),axis=1)]

        print('Number of peptides after keeping peptides with lenght from {} to {} : {}'.format(minLen, maxLen, temp.shape[0]))

#       Final check: validate amino acids
        invalid_peptides = temp[~temp['Peptide'].apply(lambda p: set(p).issubset(STANDARD_AA))]
        if not invalid_peptides.empty:
            raise Exception(f"File '{file_name}' contains peptides with invalid amino acids.")
    
        samples[file_name] = temp

    return samples

def omitPTMContent(x):
    if re.search(r'[(].+[)]',x) != None:
        x = re.sub(regex_find_cureved_brackets,'' ,x)

    if re.search(r'[[].+[]]',x) != None:
        x = re.sub(regex_find_square_brackets,'' ,x)

    return x

def saveNmerData(location, samples, peptideLength = 9, unique = True):

    for file_name, data in samples.items():
        for replicate_name, replicate_data in data.items():

            # Keeping only unique peptides
            if unique == True:
                replicate_data = replicate_data.drop_duplicates('Peptide', keep='first')
                file_extension = 'mer.txt'
            else:
                file_extension = 'merwithduplicates.txt'

            if type(peptideLength) == int:
                replicate_data[replicate_data.Length == peptideLength]['Peptide'].to_csv(os.path.join(location, file_name, replicate_name[:-4]+'_'+str(peptideLength)+file_extension), header=False, index=False)
            else:
                replicate_data[replicate_data['Length'].between(peptideLength[0], peptideLength[1], inclusive='both')]['Peptide'].to_csv(os.path.join(location, file_name, replicate_name[:-4]+'_'+str(peptideLength[0])+'to'+str(peptideLength[1])+file_extension), header=False, index=False)


def getSeqLogosImages(samples_data, task_id, motif_length, logger):
    seqlogos = {}

    for sample, replicates in samples_data.items():
        sample_dir = os.path.join(data_mount, task_id, sample)
        pattern = f'*_{motif_length}mer.txt'
        matched_files = glob.glob(os.path.join(sample_dir, pattern))

        if matched_files:
            peptide_file = matched_files[0]
            try:
                with open(peptide_file, 'r') as f:
                    lines = [line for line in f if line.strip()]
                    num_peptides = len(lines)
            except Exception as e:
                logger.exception(f"Failed to read peptide file: {peptide_file}")
                num_peptides = 0
        else:
            logger.warning(f"No peptide file found for sample={sample}, motif_length={motif_length}")
            num_peptides = 0

        # Create the list of logo image filenames + peptide count
        seqlogos[sample] = [
            [replicate[:-4] + '-001.jpg', num_peptides]
            for replicate in sorted(replicates)
        ]

    return seqlogos

def getGibbsImages(logger, taskId, samples_data):
    logger.info(f'getGibbsImages method called with taskId: {taskId} and {len(samples_data)} samples.')
    
    gibbsImages = {}

    os.chdir(project_root)

    # This approach has to be modified as the cluster is picked from the files(JPG) present in results.
    # It should be linked with gibbscluster directly to get the results.
    for sample, replicates in dict(sorted(samples_data.items())).items():
        gibbsImages[sample] = dict()

        for replicate in sorted(replicates.keys()):

            logger.info(f'Path for Barplots: app/static/images/{taskId}/{sample}/gibbscluster/{replicate[:-4]}/*/images/*.png')

            bar_plot = [x[len('app/static/'):] for x in glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate[:-4]}/*/images/*.barplot.png')]

            # Processing only if Bar Plot was generated for the input
            if len(bar_plot) == 0:
                logger.warning(f'No bar plot found at the directory for sample {sample}, replicate {replicate}')
                continue

            # Finding the best cluster
            tab_files = glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate[:-4]}/*/images/gibbs.KLDvsClusters.tab')
            bestCluster = pd.read_table(tab_files[0])
            bestCluster = bestCluster[bestCluster.columns].sum(axis=1).idxmax()

            clusters = [[x[len('app/static/'):], "Number of peptides in core could not be calculated", [], None, None] for x in sorted(glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate[:-4]}/*/logos/gibbs_logos_*of{bestCluster}*-001.png'))]

            # Finding the number of records used for the cluster
            findNumberOfPeptidesInCore(clusters, taskId, sample, replicate)

            # Append predicted allele information to the object
            appendPredictedAllelesInfo(clusters, taskId, sample, replicate)

            # Updating gibbsImages
            gibbsImages[sample][replicate[:-4]] = dict()
            gibbsImages[sample][replicate[:-4]][bar_plot[0]] = clusters

    return gibbsImages


def getGibbsImagesAll(logger, taskId, samples_data):
    """Like getGibbsImages but returns every available cluster count for offline export."""
    logger.info(f'getGibbsImagesAll called with taskId: {taskId}')
    os.chdir(project_root)

    gibbs_all = {}

    for sample, replicates in dict(sorted(samples_data.items())).items():
        gibbs_all[sample] = {}

        for replicate in sorted(replicates.keys()):
            bar_plots = [x[len('app/static/'):] for x in glob.glob(
                f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate[:-4]}/*/images/*.barplot.png')]
            if not bar_plots:
                continue

            bar_plot = bar_plots[0]

            tab_files = glob.glob(
                f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate[:-4]}/*/images/gibbs.KLDvsClusters.tab')
            if not tab_files:
                continue
            best_df = pd.read_table(tab_files[0])
            best_cluster = best_df[best_df.columns].sum(axis=1).idxmax()

            # Determine max cluster count from logo filenames
            all_logos = glob.glob(
                f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate[:-4]}/*/logos/gibbs_logos_*-001.png')
            max_n = 0
            for lf in all_logos:
                try:
                    max_n = max(max_n, int(os.path.basename(lf).split('of')[1].split('-')[0]))
                except Exception:
                    pass

            gibbs_all[sample][replicate[:-4]] = {}

            def _make_clusters(n):
                entries = [
                    [x[len('app/static/'):], "Number of peptides in core could not be calculated", [], None, None]
                    for x in sorted(glob.glob(
                        f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate[:-4]}/*/logos/gibbs_logos_*of{n}*-001.png'))
                ]
                findNumberOfPeptidesInCore(entries, taskId, sample, replicate)
                appendPredictedAllelesInfo(entries, taskId, sample, replicate)
                return entries

            gibbs_all[sample][replicate[:-4]][''] = {'bar_plot': bar_plot, 'clusters': _make_clusters(best_cluster)}

            for n in range(1, min(max_n, 6) + 1):
                logos_n = glob.glob(
                    f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate[:-4]}/*/logos/gibbs_logos_*of{n}*-001.png')
                if logos_n:
                    gibbs_all[sample][replicate[:-4]][str(n)] = {'bar_plot': bar_plot, 'clusters': _make_clusters(n)}

    return gibbs_all


# Method to calculate the peptides present in cluster
def findNumberOfPeptidesInCore(clusters, taskId, sample, replicate):
    print(f'findNumberOfPeptidesInCore : Clusters passed={clusters}')

    for cluster in clusters:
        cluster_attempt = os.path.basename(cluster[0]).split("_")[2].split("-")[0]

        try:
            path_for_core = f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate[:-4]}/*/cores/*{cluster_attempt}*'
            print(f'findNumberOfPeptidesInCore : Searching for Core at={path_for_core}')

            core_files = glob.glob(path_for_core)
            print(f'Available cores={core_files}')

            num_peptides = pd.read_table(core_files[0], header=None).shape[0]

            cluster[1] = num_peptides  # Update peptide count
            cluster.append(cluster_attempt)  # Append cluster_attempt so next method can use it

        except Exception as e:
            print(f"Error in findNumberOfPeptidesInCore for cluster {cluster_attempt}: {e}")
            continue

    print(f'Updated clusters with peptide count and cluster_attempt: {clusters}')

# Method to append predicted allele information to the clusters (MHC-TP)
# cluster[2] is set to a list of top-3 prediction dicts: [{hla, score, motif_url}, ...]
def appendPredictedAllelesInfo(clusters, taskId, sample, replicate):
    print(f'appendPredictedAllelesInfo : Clusters passed={clusters}')

    for cluster in clusters:
        # Use the cluster_attempt passed from findNumberOfPeptidesInCore (index 5)
        cluster_attempt = cluster[5] if len(cluster) > 5 else os.path.basename(cluster[0]).split("_")[2].split("-")[0]

        try:
            path_corr_matrix = f'app/static/images/{taskId}/{sample}/hla_clust_output/{replicate[:-4]}/clust_result/corr-data/corr_matrix.csv'
            df = pd.read_csv(path_corr_matrix)

            matching_rows = df[df['Cluster'] == cluster_attempt]

            if not matching_rows.empty:
                top_rows = matching_rows.sort_values(by='Correlation', ascending=False).head(3)
                predictions = []
                for _, row in top_rows.iterrows():
                    hla = row['HLA']
                    ref_path = f'app/static/images/{taskId}/{sample}/hla_clust_output/{replicate[:-4]}/clust_result/allotypes-img/{hla}.png'
                    motif_url = (
                        f'/static/images/{taskId}/{sample}/hla_clust_output/{replicate[:-4]}/clust_result/allotypes-img/{hla}.png'
                        if os.path.exists(ref_path) else None
                    )
                    predictions.append({
                        'hla': hla,
                        'score': round(float(row['Correlation']), 2),
                        'motif_url': motif_url,
                    })
                cluster[2] = predictions
                print(f"Cluster: {cluster_attempt}, top predictions: {predictions}")
            else:
                print(f"No matching rows found for Cluster: {cluster_attempt}")

        except Exception as e:
            print(f"Error processing cluster {cluster_attempt}: {e}")
            continue

    print(f'appendPredictedAllelesInfo : Updated clusters: {clusters}')

# This method will generate the binding predictions.
# User can select the binding prediction to be used.
# User can enter the names of alleles of interest.
def generateBindingPredictions(taskId, alleles_unformatted, method, ALLELE_DICTIONARY):
    
    print('Generating Binding Predictions for task {} for {} alleles using {}.'.format(taskId,alleles_unformatted,method.short_name))

    # Ensuring program in the right directory
    os.chdir(project_root)

    # Load the allele compatibility matrix (assuming it's a CSV file)
    compatibility_matrix_path = os.path.join('app', 'static', 'images', taskId, 'allele_compatibility_matrix.csv')
    compatibility_matrix = pd.read_csv(compatibility_matrix_path, index_col=0)

    for sample in os.listdir('{}/{}'.format(data_mount,taskId)):
        if not os.path.isdir(os.path.join(data_mount, taskId, sample)):
            continue
        for replicate in os.listdir('{}/{}/{}'.format(data_mount,taskId,sample)):

            # Loading data in dataframe for the use of predictors
            if replicate[-12:] == '8to14mer.txt':
                data = pd.read_csv('{}/{}/{}/{}.csv'.format(data_mount,taskId,sample,replicate[:-13]))


            if sample != 'Control':
                if replicate[-12:] == '8to14mer.txt':

            # Loading data in dataframe for the use of predictors
                    data = pd.read_csv('{}/{}/{}/{}'.format(data_mount,taskId,sample,replicate), header=None)
                    input_peptides = data[0].tolist() 

                    # Check if the method (prediction tool) is compatible with each allele
                    if method.short_name == Class_One_Predictors.MixMHCpred.short_name:
                        for allele in alleles_unformatted.split(","):
                            # Check if the allele is compatible with the current tool
                            if compatibility_matrix.at[method.full_name, allele] == 'Yes':  # or 'No', depending on your matrix values
                                # Run the command for compatible alleles
                                mixmhc_outdir = f'{project_root}/app/static/images/{taskId}/{sample}/MixMHCpred/{replicate[:-13]}/{allele.replace(":", "_")}'
                                os.makedirs(mixmhc_outdir, exist_ok=True)
                                mixmhc_allele = get_allele_name_tool_specific(allele, 'mixMHCpred 3.0', MHC_Class.One, ALLELE_DICTIONARY)
                                print(f"  Running MixMHCpred: allele={mixmhc_allele}, input={data_mount}/{taskId}/{sample}/{replicate}")
                                mixmhc_env = os.environ.copy()
                                mixmhc_env['PATH'] = f'{project_root}/lenv/bin:' + mixmhc_env.get('PATH', '/usr/bin:/bin')
                                result = subprocess.run(
                                    [
                                        f'{project_root}/app/tools/MixMHCpred/MixMHCpred',
                                        '-i', f'{data_mount}/{taskId}/{sample}/{replicate}',
                                        '-o', f'{mixmhc_outdir}/{replicate}',
                                        '-a', mixmhc_allele
                                    ],
                                    capture_output=True, text=True, env=mixmhc_env
                                )
                                if result.returncode != 0:
                                    print(f"  MixMHCpred ERROR (rc={result.returncode}): {result.stderr[:500]}")

                    elif(method.short_name==Class_One_Predictors.NetMHCpan.short_name):

                        for allele in alleles_unformatted.split(","):
                            # Check if the allele is compatible with the current tool
                            if compatibility_matrix.at[method.full_name, allele] == 'Yes':  # or 'No', depending on your matrix values
                                # Run the command for compatible alleles
                                subprocess.run(
                                    ['{}/app/tools/netMHCpan-4.2/netMHCpan'.format(project_root),
                                    '-xls', '-BA', '-p',
                                    '{}/{}/{}/{}'.format(data_mount, taskId, sample, replicate),
                                    '-a', get_allele_name_tool_specific(allele, 'netMHCpan 4.2 b', MHC_Class.One, ALLELE_DICTIONARY),
                                    '-xlsfile', '{}/app/static/images/{}/{}/NetMHCpan/{}/{}/{}'.format(
                                        project_root, taskId, sample, replicate[:-13], allele.replace(':', '_'), replicate)],
                                    stdout=DEVNULL,  # Suppress standard output
                                )

                    # Check if the method (prediction tool) is 'MHCflurry' and process accordingly
                    if method.short_name == Class_One_Predictors.MHCflurry.short_name:
                        predictor = Class1PresentationPredictor.load()

                        for allele in alleles_unformatted.split(','):
                            # Check if the allele is compatible with MHCflurry
                            if compatibility_matrix.at[method.full_name, allele] == 'Yes':  # or 'No', depending on your matrix values
                                # Predict and save the result only for compatible alleles
                                mhc_flurry_prediction_result = predictor.predict(
                                    peptides=input_peptides,
                                    alleles=[get_allele_name_tool_specific(allele, 'MHCflurry 2.0', MHC_Class.One, ALLELE_DICTIONARY)],
                                    verbose=1
                                )
                                
                                # Save the prediction result
                                result_path = f'{project_root}/app/static/images/{taskId}/{sample}/{Class_One_Predictors.MHCflurry.short_name}/{replicate[:-13]}/{allele.replace(":", "_")}/{replicate}'
                                os.makedirs(os.path.dirname(result_path), exist_ok=True)
                                mhc_flurry_prediction_result.to_csv(result_path, index=False)

                elif replicate[-13:] == '12to20mer.txt':
                    # Check if the method (prediction tool) is 'MixMHC2pred' and process accordingly
                    if method.short_name == Class_Two_Predictors.MixMHC2pred.short_name:
                        for allele in alleles_unformatted.split(','):
                            # Check if the allele is compatible with MixMHC2pred
                            if compatibility_matrix.at[Class_Two_Predictors.MixMHC2pred.full_name, allele] == 'Yes':  # or 'No', depending on your matrix values
                                # Run MixMHC2pred-2.0 command
                                mixmhc2_outdir = f'{project_root}/app/static/images/{taskId}/{sample}/MixMHC2pred/{replicate[:-14]}/{allele.replace(":", "_")}'
                                os.makedirs(mixmhc2_outdir, exist_ok=True)
                                command = [
                                    f'{project_root}/app/tools/MixMHC2pred-2.0/MixMHC2pred_unix',
                                    '-i', f'{data_mount}/{taskId}/{sample}/{replicate}',
                                    '-o', f'{mixmhc2_outdir}/{replicate}',
                                    '-a', get_allele_name_tool_specific(allele, 'MixMHC2pred-2.0', MHC_Class.Two, ALLELE_DICTIONARY),
                                    '--no_context'
                                ]
                                call(command)

                    # Check if the method (prediction tool) is 'NetMHCpanII' and process accordingly
                    if method.short_name == Class_Two_Predictors.NetMHCpanII.short_name:
                        for allele in alleles_unformatted.split(','):
                            # Check if the allele is compatible with NetMHCpanII
                            if compatibility_matrix.at[Class_Two_Predictors.NetMHCpanII.full_name, allele] == 'Yes':  # or 'No', depending on your matrix values
                                # Prepare the command to run NetMHCpanII for compatible alleles
                                command = [
                                    f'{project_root}/app/tools/netMHCIIpan-4.3/netMHCIIpan', '-xls', '-inptype', '1',
                                    '-f', '{}/{}/{}/{}'.format(data_mount, taskId, sample, replicate),
                                    '-a', get_allele_name_tool_specific(allele, 'netMHCIIpan 4.3 e', MHC_Class.Two, ALLELE_DICTIONARY),
                                    '-xlsfile', f'{project_root}/app/static/images/{taskId}/{sample}/{Class_Two_Predictors.NetMHCpanII}/{replicate[:-14]}/{allele.replace(":", "_")}/{replicate}'
                                ]

                                # Run the command for the compatible allele
                                run(command, stdout=DEVNULL)  # Suppress standard output
        
            os.chdir(project_root)

def saveBindersData(taskId, alleles, method, mhcclass):

    # Taking out the peptides from the uploaded control data to tag binders present in control group.
    control_peptides = set()

    # Load the allele compatibility matrix (assuming it's a CSV file)
    compatibility_matrix_path = os.path.join(project_root,'app', 'static', 'images', taskId, 'allele_compatibility_matrix.csv')
    compatibility_matrix = pd.read_csv(compatibility_matrix_path, index_col=0)

    if mhcclass == MHC_Class.One:
        print("Entering IF: MHC Class I detected")
        control_replicates = glob.glob(f'{data_mount}/{taskId}/Control/*8to14mer.txt')
    elif mhcclass == MHC_Class.Two:
        print("Entering IF: MHC Class II detected")
        control_replicates = glob.glob(f'{data_mount}/{taskId}/Control/*12to20mer.txt')

    if len(control_replicates) != 0:
        print(f"Entering IF: Found {len(control_replicates)} control replicate files")
        for control_replicate in control_replicates:
            f = open(control_replicate,'r')
            for peptide in f.readlines():
                control_peptides.add(peptide.replace("\n",""))
            f.close()

    print('Number of pre-processed peptides from control group:', len(control_peptides))

    for sample in os.listdir('{}/{}'.format(data_mount,taskId)):
        if not os.path.isdir(os.path.join(data_mount, taskId, sample)):
            continue
        for replicate in os.listdir('{}/{}/{}'.format(data_mount,taskId,sample)):
            if sample != 'Control' and (replicate[-12:] == '8to14mer.txt' or replicate[-13:]=='12to20mer.txt'):
                print(f"Entering IF: Processing sample={sample}, replicate={replicate}")

                # Original upload file used to derive all other columns present in the input file
                if replicate[-12:] == '8to14mer.txt':
                    print("Entering IF: Detected 8to14mer replicate")
                    input_file = pd.read_csv('{}/{}/{}/{}.csv'.format(data_mount,taskId,sample,replicate[:-13]))
                elif replicate[-13:]=='12to20mer.txt':
                    print("Entering IF: Detected 12to20mer replicate")
                    input_file = pd.read_csv('{}/{}/{}/{}.csv'.format(data_mount,taskId,sample,replicate[:-14]))

                # Dropping null Peptides
                input_file = input_file[input_file['Peptide'].apply(lambda x: isinstance(x, str) and x.strip() != '')]

                # Adding Colunm to represen the peptides without the PTM changes
                input_file['StrippedPeptide'] = input_file.apply(lambda x : omitPTMContent(x['Peptide']),axis=1)

                # Adding PTM detected method
                input_file['PTM detected'] = input_file.apply(lambda x: 'N' if x['Peptide'] == x['StrippedPeptide'] else 'Y', axis=1)

                # Initialsing the allele and binders collection
                alleles_dict = {}

                # MHCflurry case
                if method.short_name == Class_One_Predictors.MHCflurry.short_name:
                    print("Entering IF: Running MHCflurry case")

                    for allele in alleles.split(','):
                        if compatibility_matrix.at[Class_One_Predictors.MHCflurry.full_name, allele] == 'Yes':
                            print(f"  Compatible allele found for MHCflurry: {allele}")
                            f = pd.read_csv(f'{project_root}/app/static/images/{taskId}/{sample}/{Class_One_Predictors.MHCflurry}/{replicate[:-13]}/{allele.replace(":", "_")}/{replicate}')

                            f['Binding Level'] = ""
                            f['Control'] = ""

                            f['Binding Level'] = f['presentation_percentile'].apply(
                                lambda x: 'SB' if float(x) <= 0.2 else ('WB' if float(x) <= 2 else '')
                            )
                            f['Control'] = f['peptide'].apply(lambda x : 'Y' if x in control_peptides else '')
                            f.rename(columns={'peptide': 'StrippedPeptide'}, inplace=True)

                            outpath = 'app/static/images/{}/{}/{}/{}/binders/{}/{}_{}_{}_binders.csv'.format(
                                taskId, sample, method.short_name, replicate[:-13], allele.replace(':', '_'),
                                replicate[:-13], allele.replace(':', '_'), method.short_name)
                            f.sort_values(by=['presentation_percentile'])[['StrippedPeptide', 'presentation_percentile', 'Binding Level', 'affinity', 'Control']]\
                                .merge(input_file, on='StrippedPeptide', how='left')\
                                .to_csv(outpath, index=False)
                            print(f"Saved file: {outpath}")

                # MixMHCpred case
                if method.short_name == Class_One_Predictors.MixMHCpred.short_name:
                    print("Entering IF: Running MixMHCpred case")

                    for allele in alleles.split(','):
                        if compatibility_matrix.at[Class_One_Predictors.MixMHCpred.full_name, allele] == 'Yes':
                            print(f"  Compatible allele found for MixMHCpred: {allele}")
                            f = pd.read_csv(f'{project_root}/app/static/images/{taskId}/{sample}/MixMHCpred/{replicate[:-13]}/{allele.replace(":", "_")}/{replicate}', skiprows=11, sep='\t')

                            f['Binding Level'] = ""
                            f['Control'] = ""
                            f['Binding Level'] = f['%Rank_bestAllele'].apply(
                                lambda x: 'SB' if float(x) <= 2 else ('WB' if float(x) <= 10 else '')
                            )
                            f['Control'] = f['Peptide'].apply(lambda x : 'Y' if x in control_peptides else '')
                            f.rename(columns={'Peptide': 'StrippedPeptide'}, inplace=True)

                            outpath = f'{project_root}/app/static/images/{taskId}/{sample}/{method.short_name}/{replicate[:-13]}/binders/{allele.replace(":", "_")}/{replicate[:-13]}_{allele.replace(":", "_")}_{method.short_name}_binders.csv'
                            f.sort_values(by=['%Rank_bestAllele'])[['StrippedPeptide', '%Rank_bestAllele', 'Binding Level', 'Control']]\
                                .merge(input_file, on='StrippedPeptide', how='left')\
                                .to_csv(outpath, index=False)
                            print(f"Saved file: {outpath}")

                # MixMHC2pred case
                if method.short_name == Class_Two_Predictors.MixMHC2pred.short_name:
                    print("Entering IF: Running MixMHC2pred case")

                    for allele in alleles.split(','):
                        if compatibility_matrix.at[Class_Two_Predictors.MixMHC2pred.full_name, allele] == 'Yes':
                            print(f"  Compatible allele found for MixMHC2pred: {allele}")
                            f = pd.read_csv(f'{project_root}/app/static/images/{taskId}/{sample}/MixMHC2pred/{replicate[:-14]}/{allele.replace(":", "_")}/{replicate}', skiprows=19, sep='\t')

                            f['Binding Level'] = ""
                            f['Control'] = ""
                            f['Binding Level'] = f['%Rank_best'].apply(
                                lambda x: 'SB' if float(x) <= 2 else ('WB' if float(x) <= 10 else '')
                            )
                            f['Control'] = f['Peptide'].apply(lambda x : 'Y' if x in control_peptides else '')
                            f.rename(columns={'Peptide': 'StrippedPeptide'}, inplace=True)

                            s = f\
                                .sort_values(by=['%Rank_best'])[['StrippedPeptide','Core_best','%Rank_best','Binding Level','Control']]\
                                .merge(input_file, on='StrippedPeptide',how='left')

                            # Adding special column to hold both StrippedPeptide and Core_best
                            s['Peptides : StrippedPeptide : Core_best'] = s['Peptide'] + ' : ' + s['StrippedPeptide'] + ' : ' + s['Core_best']

                            outpath = 'app/static/images/{}/{}/{}/{}/binders/{}/{}_{}_{}_binders.csv'.format(taskId,sample,method.short_name,replicate[:-14],allele.replace(':', '_'),replicate[:-14],allele.replace(':', '_'),method.short_name)
                            s.to_csv(outpath, index=False)
                            print(f"Saved file: {outpath}")

                # NetMHCpanII case
                if method.short_name == Class_Two_Predictors.NetMHCpanII.short_name:
                    print("Entering IF: Running NetMHCpanII case")

                    for allele in alleles.split(','):
                        if compatibility_matrix.at[Class_Two_Predictors.NetMHCpanII.full_name, allele] == 'Yes':
                            print(f"  Compatible allele found for NetMHCpanII: {allele}")
                            f = pd.read_table(f'{project_root}/app/static/images/{taskId}/{sample}/{Class_Two_Predictors.NetMHCpanII}/{replicate[:-14]}/{allele.replace(":", "_")}/{replicate}', skiprows=2)

                            f['Binding Level'] = ""
                            f['Control'] = ""
                            f['Binding Level'] = f['EL_rank'].apply(
                                lambda x: 'SB' if float(x) <= 1 else ('WB' if float(x) <= 5 else '')
                            )
                            f['Control'] = f['Peptide'].apply(lambda x : 'Y' if x in control_peptides else '')
                            f.rename(columns={'Peptide': 'StrippedPeptide'}, inplace=True)

                            s = f.sort_values(by=['EL_rank'])[['StrippedPeptide','Core','EL_rank','Binding Level','Control']]\
                                .merge(input_file, on='StrippedPeptide',how='left')
                            s['Peptides : StrippedPeptide : Core'] = s['Peptide'] + ' : ' + s['StrippedPeptide'] + ' : ' + s['Core']

                            outpath = f'{project_root}/app/static/images/{taskId}/{sample}/{method.short_name}/{replicate[:-14]}/binders/{allele.replace(":", "_")}/{replicate[:-14]}_{allele.replace(":", "_")}_{method.short_name}_binders.csv'
                            s.to_csv(outpath, index=False)
                            print(f"Saved file: {outpath}")

                            nine_mer_path = os.path.join(data_mount, taskId, sample, replicate[:-14]+'_9mer.txt')
                            s[['Core']].drop_duplicates(subset='Core').to_csv(nine_mer_path, header=False, index=False)
                            print(f"Saved file: {nine_mer_path}")

                # NetMHCpan case
                if method.short_name == Class_One_Predictors.NetMHCpan.short_name:
                    print("Entering IF: Running NetMHCpan case")

                    for allele in alleles.split(','):
                        if compatibility_matrix.at[Class_One_Predictors.NetMHCpan.full_name, allele] == 'Yes':
                            print(f"  Compatible allele found for NetMHCpan: {allele}")
                            f = pd.read_table(f'{project_root}/app/static/images/{taskId}/{sample}/NetMHCpan/{replicate[:-13]}/{allele.replace(":", "_")}/{replicate}', skiprows=2)

                            f['Binding Level'] = ""
                            f['Control'] = ""
                            f['Binding Level'] = f['EL_rank'].apply(
                                lambda x: 'SB' if float(x) <= 0.5 else ('WB' if float(x) <= 2 else '')
                            )
                            f['Control'] = f['Peptide'].apply(lambda x : 'Y' if x in control_peptides else '')
                            f.rename(columns={'Peptide': 'StrippedPeptide'}, inplace=True)

                            outpath = f'{project_root}/app/static/images/{taskId}/{sample}/{method.short_name}/{replicate[:-13]}/binders/{allele.replace(":", "_")}/{replicate[:-13]}_{allele.replace(":", "_")}_{method.short_name}_binders.csv'
                            f.sort_values(by=['EL_rank'])[['StrippedPeptide', 'EL_rank', 'Binding Level', 'BA_score', 'core', 'Control']]\
                                .merge(input_file, on='StrippedPeptide', how='left')\
                                .to_csv(outpath, index=False)
                            print(f"Saved file: {outpath}")

def getPredictionResuslts(taskId,alleles,methods_passed,samples):

    alleles = alleles.split(',')

    # Using a method local copy as object is not udpated. Hence not conflicting with other funcatinlitites.
    methods = methods_passed.copy()
    # Adding Majority Voted Binder as a method. Easy to fetch from UI
    methods.append('Majority_Voted')

    predicted_binders = {}

    for sample in samples:
        if sample !='Control':

            os.chdir(project_root)
            predicted_binders[sample] = {}
            
            for allele in alleles:
                
                allele = allele.replace(':', '_')

                predicted_binders[sample][allele] = {}

                for method in methods:
                    predicted_binders[sample][allele][method] = {}

    # for sample in os.listdir(f'{data_mount}/{taskId}/'):
    for sample in samples:

        for method in methods:
            os.chdir(project_root)
            for replicate in os.listdir(f'{data_mount}/{taskId}/{sample}/'):
                if replicate[-12:] == '8to14mer.txt':
                    for allele in alleles:

                        allele = allele.replace(':', '_')
                        os.chdir('app/')
                        temp = glob.glob(f'static/images/{taskId}/{sample}/{method}/{replicate[:-13]}/binders/{allele}/*.csv')
                        if len(temp) != 0:
                            predicted_binders[sample][allele][method][replicate[:-13]]= temp[0]
                        os.chdir(project_root)

                elif replicate[-13:] == '12to20mer.txt':
                    for allele in alleles:

                        allele = allele.replace(':', '_')
                        os.chdir('app/')
                        temp = glob.glob(f'static/images/{taskId}/{sample}/{method}/{replicate[:-14]}/binders/{allele}/*.csv')
                        if len(temp) != 0:
                            predicted_binders[sample][allele][method][replicate[:-14]]= temp[0]
                        os.chdir(project_root)

    # print(predicted_binders)      
    return predicted_binders

def getPredictionResusltsForUpset(taskId,alleles,methods,samples):

    alleles = alleles.split(',')

    predicted_binders = {}


    for allele in alleles:

        allele = allele.replace(':', '_')
        
        os.chdir(project_root)

        predicted_binders[allele] = {}

        for sample in samples:
            if sample !='Control':
                predicted_binders[allele][sample] = []
                

                for replicate in os.listdir(f'{data_mount}/{taskId}/{sample}/'):
                    if replicate[-12:] == '8to14mer.txt':
                        # os.chdir('app/')
                        # temp = glob.glob(f'static/images/{taskId}/{sample}/{method}/{replicate[:-13]}/binders/{allele}/*.csv')
                        # if len(temp) != 0:
                        #     predicted_binders[allele][sample][method][replicate[:-13]]= temp[0]
                        # os.chdir(project_root)
                        predicted_binders[allele][sample].append(replicate[:-13])
                    elif replicate[-13:] == '12to20mer.txt':
                        predicted_binders[allele][sample].append(replicate[:-14])

    return predicted_binders

def getOverLapData(samples_data):

    overlap = {}

    for sample,replicates in samples_data.items():
        overlap[sample] = [replicate[:-4] for replicate, data in replicates.items()]

    return overlap         
    
def get_allele_name_tool_specific(allele, predictor, class_, df):
    # Print the parameters (input)
    print(f"Fetching allele name tool specific for allele: {allele}, predictor: {predictor}, class: {class_}")

    # Inline logic for filtering
    result = df.loc[
        (df['Allele name standardised'] == allele) &
        (df['Class'] == class_) &
        (df['Predictor'] == predictor),
        'Allele name tool specific'
    ]

    # Check the result and print
    if not result.empty:
        print(f"Found match: {result.iloc[0]}")
        return result.iloc[0]
    else:
        print(f"No match found for allele: {allele}, predictor: {predictor}, class: {class_}")
        return None

def saveMajorityVotedBinders(taskId, data, predictionTools, alleles_unformatted, ALLELE_DICTIONARY):
    # Create directories to store majority binding prediction results
    for sample, replicates in data.items():
        for predictionTool in predictionTools:
            for replicate in replicates:
                if alleles_unformatted != "":
                    for allele in alleles_unformatted.split(','):
                        try:
                            if sample != 'Control':
                                path = os.path.join(
                                    'app', 'static', 'images', taskId, sample,
                                    'Majority_Voted', replicate[:-4], 'binders',
                                    allele.replace(':', '_')
                                )
                                if not os.path.exists(path):
                                    Path(path).mkdir(parents=True, exist_ok=True)
                                    print(f"Directory Created : {path}")
                        except FileExistsError:
                            print(f"Directory already exists {path}")

    # Load allele compatibility matrix
    compatibility_matrix_path = os.path.join('app', 'static', 'images', taskId, 'allele_compatibility_matrix.csv')
    compatibility_matrix = pd.read_csv(compatibility_matrix_path, index_col=0)

    for sample, replicates in data.items():
        if sample != 'Control':
            for replicate in replicates:
                for allele in compatibility_matrix.columns:
                    # Get tools compatible with this allele
                    compatible_tools = compatibility_matrix.index[compatibility_matrix[allele] == "Yes"].tolist()

                    binder_files = []
                    for predictionTool in compatible_tools:
                        predictor_short_name = get_short_name(predictionTool)
                        allele_sanitized = allele.replace(':', '_')
                        search_path = f'app/static/images/{taskId}/{sample}/{predictor_short_name}/{replicate[:-4]}/binders/{allele_sanitized}/*.csv'
                        print(f"Searching binders is: {search_path}")

                        files = glob.glob(search_path)
                        binder_files.extend([(f, predictor_short_name) for f in files])
                        print(f"Found files: {binder_files}")

                    # Majority voting logic
                    peptide_counts = defaultdict(int)
                    all_data = []
                    extra_cols = []

                    for binder_file, tool_name in binder_files:
                        df = pd.read_csv(binder_file)

                        # Remove intermediate columns (but don't drop yet)
                        cols_to_remove = df.columns[
                            df.columns.get_loc('StrippedPeptide') + 1 : df.columns.get_loc('Control')
                        ]
                        renamed_cols = {col: f"{tool_name}_{col}" for col in cols_to_remove}
                        extra_df = df[['StrippedPeptide'] + list(cols_to_remove)].rename(columns=renamed_cols)
                        extra_cols.append(extra_df)

                        # --- Voting logic goes here (Binding Level still exists) ---
                        peptides_in_file = set(
                            df.loc[df['Binding Level'].notna() & (df['Binding Level'] != ''), 'StrippedPeptide']
                            .dropna()
                            .astype(str)
                        )
                        for peptide in peptides_in_file:
                            peptide_counts[peptide] += 1

                        # Now it's safe to drop the intermediate columns
                        df = df.drop(columns=cols_to_remove)
                        all_data.append(df)

                    majority_threshold = len(binder_files) // 2
                    majority_peptides = [
                        pep for pep, count in peptide_counts.items()
                        if count > majority_threshold
                    ]

                    # Combine all main data
                    combined_df = pd.concat(all_data, ignore_index=True)

                    combined_df['Is Majority Voted Binder'] = combined_df['StrippedPeptide'].apply(
                        lambda x: 'Y' if x in majority_peptides else 'N'
                    )

                    # Merge back the extra columns (outer join by StrippedPeptide)
                    for extra_df in extra_cols:
                        combined_df = combined_df.merge(extra_df, on='StrippedPeptide', how='left')

                    # Final deduplication — remove exact duplicate rows
                    filtered_df = combined_df.drop_duplicates()

                    output_path = os.path.join(
                        project_root, 'app', 'static', 'images', taskId, sample,
                        'Majority_Voted', replicate[:-4], 'binders', allele.replace(':', '_'),
                        f"{replicate[:-4]}_{allele.replace(':', '_')}_majority_voted_binders.csv"
                    )
                    filtered_df.to_csv(output_path, index=False)

def runHLAClust(taskId, data, species=None, use_mhc_tp_full_DB=None, mhcclass=None, logger=None):

    logger.info(f'Running HLA Clust for task {taskId}.')

    # Determine effective species for MHC-TP (Class II uses separate reference db per species)
    if mhcclass == MHC_Class.Two:
        db_species = f"{species.lower()}_classii" if species else None
        allele_file = None  # No restricted list for Class II — search full db
    else:
        db_species = species
        allele_file = os.path.join(project_root, 'app', 'static', 'mhc-tp-default-search-alleles.csv')

    # Creating directories to store majority binding prediction results
    for sample, replicates in data.items():
        for replicate in replicates:
                    try:
                        if sample != 'Control':

                            # Path to store user friendly binders data
                            path = os.path.join(project_root, 'app', 'static', 'images', taskId, sample, 'hla_clust_output', replicate[:-4])

                            Path(path).mkdir(parents=True, exist_ok=True)
                            logger.info(f'Directory Created : {path}')

                            # Running the tool for every replicate
                            gibbs_base = os.path.join(project_root, 'app', 'static', 'images', taskId, sample, 'gibbscluster', replicate[:-4])
                            gibbs_subdirs = sorted([d for d in os.listdir(gibbs_base) if os.path.isdir(os.path.join(gibbs_base, d))])
                            input_file = os.path.join(gibbs_base, gibbs_subdirs[0]) if gibbs_subdirs else gibbs_base
                            ref_file = os.path.join(project_root, 'app', 'tools', 'HLA-PepClust', 'data', 'ref_data')

                            run_clust_search(
                                input_file=input_file,
                                ref_file=ref_file,
                                output_dir=path,
                                species=species,
                                db_species=db_species,
                                use_mhc_tp_full_DB=use_mhc_tp_full_DB,
                                allele_file=allele_file,
                                logger=logger
                            )

                    except FileExistsError:
                        logger.info(f'Directory already exists {path}')

def run_clust_search(input_file, ref_file, output_dir, species, db_species=None, use_mhc_tp_full_DB=None, allele_file=None, logger=None):
    try:
        effective_species = db_species or species

        # Graceful skip if reference db not installed for this species
        db_file = os.path.join(ref_file, f'{effective_species}.db')
        if not os.path.exists(db_file):
            if logger:
                logger.warning(f"MHC-TP: reference db not found for '{effective_species}' at {db_file}. Skipping.")
            return {"skipped": f"Reference database not installed for species '{effective_species}'"}

        # Construct base command
        command = [
            f"{project_root}/app/tools/HLA-PepClust/hlapepclust-env/bin/clust-search",
            input_file,
            ref_file,
            "-im",
            "--output", output_dir,
            "--processes", str(os.cpu_count()),
            "--NumbaDB", ref_file,
            "-s", effective_species, "-t", "0.1",
        ]

        # If human Class I and restricted allele list requested, add --hla flag
        if (
            effective_species.lower() == "human"
            and use_mhc_tp_full_DB
            and use_mhc_tp_full_DB.lower() == "no"
            and allele_file
        ):
            import csv
            with open(allele_file, newline='') as csvfile:
                reader = csv.DictReader(csvfile)
                allele_list = [row['Allele name standardised'] for row in reader if row['Allele name standardised']]
                hla_types_arg = ",".join(allele_list)
                command.extend(["-hla", hla_types_arg])
                if logger:
                    logger.info(f"Using restricted allele list for Human: {hla_types_arg}")

        # Log the command
        if logger:
            logger.info(f"Running HLA Clust with command: {' '.join(command)}")

        # Run the command
        result = subprocess.run(command, capture_output=True, text=True)

        # Check for errors
        if result.returncode != 0:
            raise RuntimeError(f"Command failed: {result.stderr}")

        return {"success": "Clustering completed", "output": result.stdout}

    except Exception as e:
        return {"error": str(e)}
    
def getHLAClustResults(taskId, data):
    print(f'getMajorityBindingImages method called with taskId: {taskId}')

    bindingImages = {}

    for sample, replicates in data.items():
        if sample == 'Control':
            continue

        sample_dict = {}

        for replicate in replicates:
            path = os.path.join(
                'app', 'static', 'images', taskId, sample, 'hla_clust_output', replicate[:-4]
            )
            
            html_files = sorted(glob.glob(os.path.join(path, '**', '*result.html'), recursive=True))
            html_files = [os.path.relpath(f, 'app/static') for f in html_files]

            if not html_files:
                print(f'No result.html files found for sample {sample}, replicate {replicate}')
                continue

            sample_dict[replicate[:-4]] = html_files

        # 🔹 Only add the sample if it actually has replicates with files
        if sample_dict:
            bindingImages[sample] = sample_dict

    # 🔹 Return False if everything was empty
    if not bindingImages:
        print("bindingImages is empty")
        return False

    print("bindingImages:", bindingImages)
    return bindingImages

def generate_peptigram(csv_path, fasta_path, protein_ids, output_dir):
    """
    Generate and save a peptigram visualization for given proteins.

    Parameters:
        csv_path (str): Path to the CSV file with peptide peaks.
        fasta_path (str): Path to the FASTA file with protein sequences.
        protein_ids (list): List of protein accession IDs to visualize.
        output_dir (str): Directory to save the output PNG file.
        output_filename (str): Output file name (default: 'protein_visualization.png').
    """

    # Ensure output directory exists
    if not os.path.isdir(output_dir):
        raise FileNotFoundError(f"Output directory does not exist: {output_dir}")

    run_pipeline(
        csv_path =csv_path,
        fasta_path = fasta_path,
        output_dir= output_dir,
        protein_list= protein_ids,
        intensity_threshold = 0.0,
        min_samples = 1,
    )