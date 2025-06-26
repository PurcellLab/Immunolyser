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
from ProtPeptigram.DataProcessor import PeptideDataProcessor
from ProtPeptigram.viz import ImmunoViz

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
            title = 'The relative frequency distribution of the peptide lengths'
            yaxis_label = '% Peptides'
            file_suffix = 'percentage'
        else:
            for replicate, data in sample.items():
                peptideProportion[replicate] = data.groupby('Length').count()['Peptide']
            title = 'The frequency distribution of the peptide lengths'
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
            # errors = pd.concat(peptideProportion, axis=1).apply(lambda x: stdev(x), axis=1)
            fig.add_trace(go.Bar(
                x=bardatacombined['Length'],
                y=bardatacombined['Count'],
                name=sample_name,
                error_y=dict(
                    type='data',
                    # array=errors,  # Optional
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
        title=title,
        xaxis=dict(title='<i>Length</i>'),
        yaxis=dict(title=f'<i>{yaxis_label}</i>')
    )

    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    return graphJSON
    
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
        temp['Peptide'] = temp.apply(lambda x : omitPTMContent(x['Peptide']),axis=1)

#       Generating Length Column
        temp['Length']= temp.apply(lambda x : len(x['Peptide']), axis=1)

#       Filtering on the basis of the peptide length
        temp = temp[temp.apply(lambda x : x['Length'] in range(minLen,maxLen),axis=1)]

        print('Number of peptides after keeping peptides with lenght from {} to {} : {}'.format(minLen, maxLen, temp.shape[0]))

#       Filtering the control peptides out
        # temp = temp[temp.apply(lambda x : x['Peptide'] not in control_peptides,axis=1)]        
    
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

            logger.info(f'Path for Barplots: app/static/images/{taskId}/{sample}/gibbscluster/{replicate[:-4]}/images/*.JPG')

            bar_plot = [os.path.basename(x) for x in glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate[:-4]}/images/*.JPG')]
            
            # Processing only if Bar Plot was generated for the input
            if len(bar_plot) == 0:
                logger.warning(f'No bar plot found at the directory for sample {sample}, replicate {replicate}')
                continue
            
            # Finding the best cluster
            bestCluster = pd.read_table(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate[:-4]}/images/gibbs.KLDvsClusters.tab')
            bestCluster = bestCluster[bestCluster.columns].sum(axis=1).idxmax()

            clusters = [[os.path.basename(x), "Number of peptides in core could not be calculated", "Allele not predicted", "PCC score associated with allele", "Url of reference motif"] for x in sorted(glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate[:-4]}/logos/gibbs_logos_*of{bestCluster}*.jpg'))]

            # Finding the number of records used for the cluster
            findNumberOfPeptidesInCore(clusters, taskId, sample, replicate)

            # Append predicted allele information to the object
            appendPredictedAllelesInfo(clusters, taskId, sample, replicate)

            # Updating gibbsImages
            gibbsImages[sample][replicate[:-4]] = dict()
            gibbsImages[sample][replicate[:-4]][bar_plot[0]] = clusters

    return gibbsImages

# Method to calculate the peptides present in cluster
def findNumberOfPeptidesInCore(clusters, taskId, sample, replicate):

    print(f'findNumberOfPeptidesInCore : Clusters passed={clusters}')

    # Finding the number of records used for the cluster
    for cluster in clusters:
        cluster_attempt = os.path.basename(cluster[0]).split("_")[2].split("-")[0]

        try :
            path_for_core = f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate[:-4]}/cores/*{cluster_attempt}*'
            print(f'findNumberOfPeptidesInCore : Cluster={clusters}, Searching for Core at={path_for_core}')

            available_cores = [os.path.basename(x) for x in glob.glob(path_for_core)]
            print(f'findNumberOfPeptidesInCore : Cluster={clusters}, Available cores={available_cores}')
            core = available_cores[0]
            cluster[1] = pd.read_table(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate[:-4]}/cores/{core}', header=None).shape[0]
        except Exception as e:
            print(e)
            continue

    print(f'findNumberOfPeptidesInCore : Updaeted gibbs clusters information: {clusters}')

# Method to append predicted allele information to the clusters (MHC-TP)
def appendPredictedAllelesInfo(clusters, taskId, sample, replicate):
    print(f'appendPredictedAllelesInfo : Clusters passed={clusters}')

    # Finding the number of records used for the cluster
    for cluster in clusters:
        cluster_attempt = os.path.basename(cluster[0]).split("_")[2].split("-")[0]  # It will store, say '1of2'

        try:
            # Read csv file generated by MHC-TP
            path_corr_matrix = f'app/static/images/{taskId}/{sample}/hla_clust_output/{replicate[:-4]}/clust_result/corr-data/corr_matrix.csv'

            # Read the CSV into a pandas DataFrame
            df = pd.read_csv(path_corr_matrix)

            # Filter the DataFrame for the rows where Cluster matches cluster_attempt
            matching_rows = df[df['Cluster'] == cluster_attempt]

            if not matching_rows.empty:
                # Sort the matching rows by Correlation in descending order and pick the first one
                highest_corr_row = matching_rows.sort_values(by='Correlation', ascending=False).iloc[0]

                # Extract the HLA and Correlation values
                hla = highest_corr_row['HLA']
                correlation = highest_corr_row['Correlation']
                cluster[2] = hla
                cluster[3] = correlation
                print(f"Cluster: {cluster_attempt}, HLA: {hla}, Correlation: {correlation}")

                # Extract the URL of the reference motif
                reference_motif_url = [f"/static/images/{taskId}/{sample}/hla_clust_output/{replicate[:-4]}/clust_result/allotypes-img/{os.path.basename(x)}" for x in glob.glob(f'app/static/images/{taskId}/{sample}/hla_clust_output/{replicate[:-4]}/clust_result/allotypes-img/{hla}.png')]
                cluster[4] = reference_motif_url[0] if reference_motif_url else "No URL found"
                print(f"Reference motif URL: {cluster[4]}")


            else:
                print(f"No matching rows found for Cluster: {cluster_attempt}")

        except Exception as e:
            print(f"Error processing cluster {cluster_attempt}: {e}")
            continue

    print(f'appendPredictedAllelesInfo : Updated gibbs clusters information: {clusters}')


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
                                call(
                                    [
                                        f'{project_root}/app/tools/MixMHCpred/MixMHCpred',
                                        '-i', f'{data_mount}/{taskId}/{sample}/{replicate}',
                                        '-o', f'{project_root}/app/static/images/{taskId}/{sample}/MixMHCpred/{replicate[:-13]}/{allele.replace(":", "_")}/{replicate}',
                                        '-a', get_allele_name_tool_specific(allele, 'mixMHCpred 3.0', MHC_Class.One, ALLELE_DICTIONARY)
                                    ]
                                )

                    elif(method.short_name==Class_One_Predictors.NetMHCpan.short_name):

                        for allele in alleles_unformatted.split(","):
                            # Check if the allele is compatible with the current tool
                            if compatibility_matrix.at[method.full_name, allele] == 'Yes':  # or 'No', depending on your matrix values
                                # Run the command for compatible alleles
                                run(
                                    ['{}/app/tools/netMHCpan-4.1/netMHCpan'.format(project_root), '-xls', '-p', 
                                    '{}/{}/{}/{}'.format(data_mount, taskId, sample, replicate),
                                    '-a', get_allele_name_tool_specific(allele, 'netMHCpan 4.1 b', MHC_Class.One, ALLELE_DICTIONARY),
                                    '-xlsfile', '{}/app/static/images/{}/{}/NetMHCpan/{}/{}/{}'.format(project_root,taskId, sample, replicate[:-13], allele.replace(':', '_'), replicate)],
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
                                # Define result path
                                result_path = f'{project_root}/app/static/images/{taskId}/{sample}/{Class_One_Predictors.MHCflurry.short_name}/{replicate[:-13]}/{allele.replace(":", "_")}/{replicate}'
                                mhc_flurry_prediction_result.to_csv(result_path, index=False)

                elif replicate[-13:] == '12to20mer.txt':
                    # Check if the method (prediction tool) is 'MixMHC2pred' and process accordingly
                    if method.short_name == Class_Two_Predictors.MixMHC2pred.short_name:
                        for allele in alleles_unformatted.split(','):
                            # Check if the allele is compatible with MixMHC2pred
                            if compatibility_matrix.at[Class_Two_Predictors.MixMHC2pred.full_name, allele] == 'Yes':  # or 'No', depending on your matrix values
                                # Prepare the command to run MixMHC2pred for compatible alleles
                                # Run MixMHC2pred-2.0 command
                                command = [
                                    f'{project_root}/app/tools/MixMHC2pred-2.0/MixMHC2pred_unix',
                                    '-i', f'{data_mount}/{taskId}/{sample}/{replicate}',
                                    '-o', f'{project_root}/app/static/images/{taskId}/{sample}/MixMHC2pred/{replicate[:-14]}/{allele.replace(":", "_")}/{replicate}',
                                    '-a', get_allele_name_tool_specific(allele, 'MixMHC2pred-2.0', MHC_Class.Two, ALLELE_DICTIONARY),
                                    '--no_context'
                                ]
                                # Run the command for the compatible allele
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
        control_replicates = glob.glob(f'{data_mount}/{taskId}/Control/*8to14mer.txt')
    elif mhcclass == MHC_Class.Two:
        control_replicates = glob.glob(f'{data_mount}/{taskId}/Control/*12to20mer.txt')

    if len(control_replicates) != 0:
        for control_replicate in control_replicates:
            f = open(control_replicate,'r')
            
            for peptide in f.readlines():
                control_peptides.add(peptide.replace("\n",""))
            
            f.close()

    print('Number of pre-processed peptides from control group:', len(control_peptides))

    for sample in os.listdir('{}/{}'.format(data_mount,taskId)):
        for replicate in os.listdir('{}/{}/{}'.format(data_mount,taskId,sample)):
            if sample != 'Control' and (replicate[-12:] == '8to14mer.txt' or replicate[-13:]=='12to20mer.txt'):

                # Original upload file used to derive all other columns present in the input file
                if replicate[-12:] == '8to14mer.txt':
                    input_file = pd.read_csv('{}/{}/{}/{}.csv'.format(data_mount,taskId,sample,replicate[:-13]))
                elif replicate[-13:]=='12to20mer.txt':
                    input_file = pd.read_csv('{}/{}/{}/{}.csv'.format(data_mount,taskId,sample,replicate[:-14]))

                # Dropping null Peptides
                input_file = input_file[input_file['Peptide'].apply(lambda x: isinstance(x, str) and x.strip() != '')]

                # Adding Colunm to represen the peptides without the PTM changes
                input_file['PlainPeptide'] = input_file.apply(lambda x : omitPTMContent(x['Peptide']),axis=1)

                # Adding PTM detected method
                input_file['PTM detected'] = input_file.apply(lambda x: 'N' if x['Peptide'] == x['PlainPeptide'] else 'Y', axis=1)


                # Initialsing the allele and binders collection
                alleles_dict = {}
                # MHCflurry case
                if method.short_name == Class_One_Predictors.MHCflurry.short_name:

                    for allele in alleles.split(','):
                        
                        # Check if the allele is compatible with MHCflurry
                        if compatibility_matrix.at[Class_One_Predictors.MHCflurry.full_name, allele] == 'Yes':  # or 'No', depending on your matrix values
                            
                            f = pd.read_csv(f'{project_root}/app/static/images/{taskId}/{sample}/{Class_One_Predictors.MHCflurry}/{replicate[:-13]}/{allele.replace(":", "_")}/{replicate}')

                            f['Binding Level'] = ""
                            f['Control'] = ""

                            # Filtering peptides with presentation_percentile (%) less than or equal to 2
                            f = f[f.apply(lambda x : float(x['presentation_percentile']) <= 2, axis=1)]

                            # Tagging each binder as SB (Strong binder) or WB (Weak binder)
                            f['Binding Level'] = f['presentation_percentile'].apply(lambda x : 'SB' if float(x) <= 0.2 else 'WB')

                            # Tagging binders present in control group
                            f['Control'] = f['peptide'].apply(lambda x : 'Y' if x in control_peptides else '')

                            # Renaming the peptide column to PlainPeptide
                            f.rename(columns={'peptide': 'PlainPeptide'}, inplace=True)

                            # Saving the filtered and tagged binders to a CSV
                            f\
                                .sort_values(by=['presentation_percentile'])[['PlainPeptide', 'presentation_percentile', 'Binding Level', 'Control']]\
                                .merge(input_file, on='PlainPeptide', how='left')\
                                .to_csv('app/static/images/{}/{}/{}/{}/binders/{}/{}_{}_{}_binders.csv'.format(taskId, sample, method.short_name,
                                                                                                        replicate[:-13], allele.replace(':', '_'),
                                                                                                        replicate[:-13], allele.replace(':', '_'), method.short_name), index=False)

                # MixMHCpred case
                if method.short_name == Class_One_Predictors.MixMHCpred.short_name:

                    for allele in alleles.split(','):
                        
                        # Check if the allele is compatible with MixMHCpred
                        if compatibility_matrix.at[Class_One_Predictors.MixMHCpred.full_name, allele] == 'Yes':  # or 'No', depending on your matrix values

                            f = pd.read_csv(f'{project_root}/app/static/images/{taskId}/{sample}/MixMHCpred/{replicate[:-13]}/{allele.replace(":", "_")}/{replicate}', skiprows=11, sep='\t')

                            f['Binding Level'] = ""
                            f['Control'] = ""
                                            
                            # Keep only strong and weak binders
                            f = f[f.apply(lambda x : float(x['%Rank_bestAllele']) <= 10, axis=1)]

                            # Tagging each binder as SB (Strong binder) or WB (Weak binder)
                            f['Binding Level'] = f['%Rank_bestAllele'].apply(lambda x : 'SB' if float(x) <= 2 else 'WB')
                                            
                            # Tagging binders present in control group
                            f['Control'] = f['Peptide'].apply(lambda x : 'Y' if x in control_peptides else '')

                            # Renaming the peptide column to PlainPeptide
                            f.rename(columns={'Peptide': 'PlainPeptide'}, inplace=True)

                            f.sort_values(by=['%Rank_bestAllele'])[['PlainPeptide', '%Rank_bestAllele', 'Binding Level', 'Control']] \
                                .merge(input_file, on='PlainPeptide', how='left') \
                                .to_csv(f'{project_root}/app/static/images/{taskId}/{sample}/{method.short_name}/{replicate[:-13]}/binders/{allele.replace(":", "_")}/{replicate[:-13]}_{allele.replace(":", "_")}_{method.short_name}_binders.csv', index=False)

                # MixMHC2pred case
                if method.short_name == Class_Two_Predictors.MixMHC2pred.short_name:

                    for allele in alleles.split(','):
                            
                        # Check if the allele is compatible with MixMHCpred
                        if compatibility_matrix.at[Class_Two_Predictors.MixMHC2pred.full_name, allele] == 'Yes':  # or 'No', depending on your matrix values
                    
                            f = pd.read_csv(f'{project_root}/app/static/images/{taskId}/{sample}/MixMHC2pred/{replicate[:-14]}/{allele.replace(":", "_")}/{replicate}', skiprows=19, sep='\t')

                            f['Binding Level'] = ""
                            f['Control'] = ""
                                
                            # Keeping strong and weak binders only
                            f = f[f.apply(lambda x : float(x['%Rank_best'])<=10,axis=1)]

                            # Tagging each binder as SB(Strong binder) or WB(Weak binder)
                            f['Binding Level'] = f['%Rank_best'].apply(lambda x : 'SB' if float(x)<=2 else 'WB')
                                
                            # Tagging binders present in control group
                            f['Control'] = f['Peptide'].apply(lambda x : 'Y' if x in control_peptides else '')

                            # Updating the name of binding results column Peptide to PlainPeptide
                            f.rename(columns={'Peptide': 'PlainPeptide'}, inplace=True)

                # NetMHCpanII case
                if method.short_name == Class_Two_Predictors.NetMHCpanII:

                    for allele in alleles.split(','):
                        
                        # Check if the allele is compatible with MixMHCpred
                        if compatibility_matrix.at[Class_Two_Predictors.NetMHCpanII.full_name, allele] == 'Yes':  # or 'No', depending on your matrix values

                            f = pd.read_table(f'{project_root}/app/static/images/{taskId}/{sample}/{Class_Two_Predictors.NetMHCpanII}/{replicate[:-14]}/{allele.replace(":", "_")}/{replicate}', skiprows=1)

                            f['Binding Level'] = ""
                            f['Control'] = ""
                        
                            # Keeping strong and weak binders only
                            f = f[f.apply(lambda x : float(x['Rank'])<=5,axis=1)]
                            
                            f['Binding Level'] = f['Rank'].apply(lambda x : 'SB' if float(x)<=1 else 'WB')

                            # Tagging binders present in control group
                            f['Control'] = f['Peptide'].apply(lambda x : 'Y' if x in control_peptides else '')

                            # Updating the name of binding results column Peptide to PlainPeptide
                            f.rename(columns={'Peptide': 'PlainPeptide'}, inplace=True)

                            f.sort_values(by=['Rank'])[['PlainPeptide', 'Rank', 'Binding Level', 'Control']] \
                                .merge(input_file, on='PlainPeptide', how='left') \
                                .to_csv(f'{project_root}/app/static/images/{taskId}/{sample}/{method.short_name}/{replicate[:-14]}/binders/{allele.replace(":", "_")}/{replicate[:-13]}_{allele.replace(":", "_")}_{method.short_name}_binders.csv', index=False)

                            s = f\
                                .sort_values(by=['Rank'])[['PlainPeptide','Core','Rank','Binding Level','Control']]\
                                .merge(input_file, on='PlainPeptide',how='left')

                            # Adding special column to hold both PlainPeptide and Core_best
                            s['Peptides : PlainPeptide : Core'] = s['Peptide'] + ' : ' + s['PlainPeptide'] + ' : ' + s['Core']

                            s.to_csv(f'{project_root}/app/static/images/{taskId}/{sample}/{method.short_name}/{replicate[:-14]}/binders/{allele.replace(":", "_")}/{replicate[:-14]}_{allele.replace(":", "_")}_{method.short_name}_binders.csv', index=False)

                            # Saving the predicted core and saving it in 9mer file
                            s[['Core']]\
                                .drop_duplicates(subset='Core')\
                                .to_csv(os.path.join(data_mount, taskId, sample, replicate[:-14]+'_9mer.txt'), header=False, index=False)

                # netMHCpan case
                if method.short_name == Class_One_Predictors.NetMHCpan.short_name:

                    for allele in alleles.split(','):

                        # Check if the allele is compatible with MixMHCpred
                        if compatibility_matrix.at[Class_One_Predictors.NetMHCpan.full_name, allele] == 'Yes':  # or 'No', depending on your matrix values

                            f = pd.read_table(f'{project_root}/app/static/images/{taskId}/{sample}/NetMHCpan/{replicate[:-13]}/{allele.replace(":", "_")}/{replicate}', skiprows=1)

                            f['Binding Level'] = ""
                            f['Control'] = ""
                        
                            # Keeping strong and weak binders only
                            f = f[f.apply(lambda x : float(x['EL_Rank'])<=2,axis=1)]
                            
                            f['Binding Level'] = f['EL_Rank'].apply(lambda x : 'SB' if float(x)<=0.5 else 'WB')

                            # Tagging binders present in control group
                            f['Control'] = f['Peptide'].apply(lambda x : 'Y' if x in control_peptides else '')

                            # Updating the name of binding results column Peptide to PlainPeptide
                            f.rename(columns={'Peptide': 'PlainPeptide'}, inplace=True)

                            f.sort_values(by=['EL_Rank'])[['PlainPeptide', 'EL_Rank', 'Binding Level', 'Control']] \
                                .merge(input_file, on='PlainPeptide', how='left') \
                                .to_csv(f'{project_root}/app/static/images/{taskId}/{sample}/{method.short_name}/{replicate[:-13]}/binders/{allele.replace(":", "_")}/{replicate[:-13]}_{allele.replace(":", "_")}_{method.short_name}_binders.csv', index=False)

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

    # Creating directories to store majority binding prediction results
    for sample, replicates in data.items():
        for predictionTool in predictionTools:
            for replicate in replicates:
                if alleles_unformatted != "":
                    for allele in alleles_unformatted.split(','):
                        try:
                            if sample != 'Control':

                                # Path to store user friendly binders data
                                path = os.path.join('app', 'static', 'images', taskId, sample, 'Majority_Voted', replicate[:-4], 'binders', allele.replace(':', '_'))    

                            if not os.path.exists(path):
                                # os.makedirs(directory)
                                Path(path).mkdir(parents=True, exist_ok=True)
                                print("Directory Created : {}".format(path))
                                
                        except FileExistsError:
                            print("Directory already exists {}".format(path))

    # Creating CSV files to store majority binding prediction results

    # Load the allele compatibility matrix (assuming it's a CSV file)
    compatibility_matrix_path = os.path.join('app', 'static', 'images', taskId, 'allele_compatibility_matrix.csv')
    compatibility_matrix = pd.read_csv(compatibility_matrix_path, index_col=0)

    for sample, replicates in data.items():

        for replicate in replicates:

            # Iterating selected alleles
            for allele in compatibility_matrix.columns:

                    # Filter tools that say "Yes" for this allele
                    compatible_tools = compatibility_matrix.index[compatibility_matrix[allele] == "Yes"].tolist()

                    binder_files = []
                    for predictionTool in compatible_tools:

                        # Store the short name in a variable
                        predictor_short_name = get_short_name(predictionTool)
                        # Remove : from allele name
                        allele = allele.replace(':', '_')

                        search_path = f'app/static/images/{taskId}/{sample}/{predictor_short_name}/{replicate[:-4]}/binders/{allele}/*.csv'
                        print(f"Searching binders is: {search_path}")  # Print where it's searching
                        
                        files = glob.glob(search_path)

                        binder_files.extend(files)  

                    print(f"  Found files: {binder_files}")  # Print found files

                    # Dictionary to hold the peptide occurrences across files
                    peptide_counts = defaultdict(int)
                    all_data = []
                    # Perform majority voting based and save in files
                    # Iterate over all files in the list
                    for binder_file in binder_files:
                        # Read the CSV file into a DataFrame
                        df = pd.read_csv(binder_file)
                        
                        # Remove unwanted columns between 'PlainPeptide' and 'Control'
                        columns_to_remove = df.columns[df.columns.get_loc('PlainPeptide')+1 : df.columns.get_loc('Control')]
                        df.drop(columns=columns_to_remove, inplace=True)
                        
                        # Add data from the file to the all_data list
                        all_data.append(df)
                        
                        # Track peptide occurrences
                        for peptide in df['PlainPeptide']:
                            peptide_counts[peptide] += 1
                    
                    # Determine majority threshold (e.g., for 3 files, threshold is 2)
                    majority_threshold = len(binder_files) // 2

                    # Filter peptides that meet the majority voting rule
                    majority_peptides = [peptide for peptide, count in peptide_counts.items() if count > majority_threshold]
                    
                    # Concatenate the data from all files
                    combined_df = pd.concat(all_data, ignore_index=True)
                    
                    # Filter the rows to keep only the majority peptides
                    filtered_df = combined_df[combined_df['PlainPeptide'].isin(majority_peptides)].drop_duplicates(subset=['PlainPeptide'])

                    # Write the filtered data to a new CSV file
                    filtered_df.to_csv(f'{project_root}/app/static/images/{taskId}/{sample}/Majority_Voted/{replicate[:-4]}/binders/{allele.replace(":", "_")}/{replicate[:-4]}_{allele.replace(":", "_")}_majority_voted_binders.csv', index=False)

def runHLAClust(taskId, data, species=None, logger=None):

    logger.info(f'Running HLA Clust for task {taskId}.')

    # Creating directories to store majority binding prediction results
    for sample, replicates in data.items():
        for replicate in replicates:
                    try:
                        if sample != 'Control':

                            # Path to store user friendly binders data
                            path = os.path.join(project_root, 'app', 'static', 'images', taskId, sample, 'hla_clust_output', replicate[:-4])    

                            # Running the tool for every replicate
                            input_file = os.path.join(project_root, 'app', 'static', 'images', taskId, sample, 'gibbscluster', replicate[:-4])
                            ref_file = os.path.join(project_root, 'app', 'tools', 'HLA-PepClust', 'data', 'ref_data')
                            output_dir = path

                            run_clust_search(input_file=input_file, ref_file=ref_file, output_dir=output_dir, species=species,logger=logger)

                        if not os.path.exists(path):
                            # os.makedirs(directory)
                            Path(path).mkdir(parents=True, exist_ok=True)
                            logger.info(f'Directory Created : {path}')
                            
                    except FileExistsError:
                        logger.info(f'Directory already exists {path}')

def run_clust_search(input_file, ref_file, output_dir, species, logger=None):
    try:
        # Validate species and set the corresponding flag
        species_flag = []
        if species.lower() == "mouse":
            species_flag = ["-s", "mouse"]
        elif species.lower() == "human":
            species_flag = ["-s", "human"]
        else:
            raise ValueError("Invalid species. Choose either 'mouse' or 'human'.")

        # Construct base command
        command = [
            f"{project_root}/app/tools/HLA-PepClust/hlapepclust-env/bin/clust-search",
            input_file,
            ref_file,
            "-im",
            "--output", output_dir,
            "--processes", str(os.cpu_count())
        ] + species_flag  # Append species flag if applicable

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

        bindingImages[sample] = {}

        for replicate in replicates:
            path = os.path.join('app', 'static', 'images', taskId, sample, 'hla_clust_output', replicate[:-4])
            
            # Collect PNG images
            png_files = sorted(glob.glob(os.path.join(path, '**', '*result.html'), recursive=True))
            png_files = [os.path.relpath(f, 'app/static') for f in png_files]  # Relative paths

            if not png_files:
                print(f'No PNG files found for sample {sample}, replicate {replicate}')
                continue

            bindingImages[sample][replicate[:-4]] = png_files

    return bindingImages

def generate_peptigram(csv_path, fasta_path, protein_ids, output_dir, output_filename="protein_visualization.png"):
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

    # Initialize and load data
    processor = PeptideDataProcessor()
    processor.load_peaks_data(csv_path)
    processor.load_protein_sequences(fasta_path)

    # Process data
    formatted_data = processor.filter_and_format_data(
        filter_contaminants=True,
        intensity_threshold=0,
        min_samples=1
    )

    # Create visualizations
    viz = ImmunoViz(formatted_data)
    fig, _ = viz.plot_peptigram(
        protein_ids=protein_ids,
        group_by="Sample",
        color_by="protein",
        protein_cmap="Set1",
        title="HLA peptide source protein visualization with highlighted highest density regions",
    )

    # Save output
    output_path = os.path.join(output_dir, output_filename)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Peptigram saved to: {output_path}")