from pandas.tseries.offsets import Second
from app import app
from flask import render_template, request, jsonify
from app import sample
from app.forms import InitialiserForm, ParentForm
import pandas as pd
import os
import subprocess
from werkzeug.utils import secure_filename
from app.utils import plot_lenght_distribution, filterPeaksFile, saveNmerData, getSeqLogosImages, getGibbsImages, generateBindingPredictions, saveBindersData, getPredictionResuslts, getPredictionResusltsForUpset, findNumberOfPeptidesInCore, getOverLapData
from app.sample import Sample
import time
from pathlib import Path
import glob
import shutil
from app.Pepscan import PepScan
import re
from collections import Counter,OrderedDict
import uuid

project_root = os.path.dirname(os.path.realpath(os.path.join(__file__, "..")))

# DEMO Task ID
DEMO_TASK_ID = app.config['DEMO_TASK_ID']

data_mount = app.config['IMMUNOLYSER_DATA']

@app.route("/initialiser", methods=["POST", "GET"])
def initialiser():    
    samples = []

    # Have to take this input from user
    maxLen = 30
    minLen = 5

    taskId = getTaskId()
    dirName = os.path.join(data_mount, taskId)
    try:
        # Create target Directory
        os.makedirs(dirName)
        print("Directory " , dirName ,  " Created ") 
    except FileExistsError:
        print("Directory " , dirName ,  " already exists")


    data = {}
    control = list()

    # Creating folders to store images
    for key, value  in request.files.items():
        sample_name = request.form[key]

        is_valid, message = validate_sample_name(sample_name)

        if not is_valid:
            print("Sample name is not valid. " + message)

        # Skipping Control Data
        # if sample_name == "Control":
        #     continue

        # Creating sub directories to store sample data
        try:
            # for seqlogos
            path_for_logos = os.path.join('app', 'static', 'images', taskId, sample_name, 'seqlogos')
            if not os.path.exists(path_for_logos):
                # os.makedirs(directory)
                Path(path_for_logos).mkdir(parents=True, exist_ok=True)
                print("Directory Created: ", path_for_logos) 
        except FileExistsError:
            print("Directory already exists")    

    # Saving the data and loading into the dictionary
    for key, value  in request.files.items():
        sample_name = request.form[key]

        # Creating sub directories to store sample data
        try:
            os.mkdir(os.path.join(dirName, sample_name))
            print("Directory Created: ", os.path.join(dirName, sample_name)) 
        except FileExistsError:
            print("Directory already exists")    

        # Not including the control group in data dict 
        # if sample_name != "Control":
        data[sample_name] = list()

        replicates = request.files.getlist(key)
        print("replics name : {}".format(replicates))
        for replicate in replicates:
            file_filename = secure_filename(replicate.filename)

            # If there is no control file uploaded then there is no point to save it.
            if file_filename != "":
                replicate.save(os.path.join(dirName, sample_name, file_filename))

            # Not including the control group in sample data dict
            # if sample_name != "Control":
                data[sample_name].append(file_filename)
            # elif file_filename !="":
                # control.append(file_filename)


        # If control data is not uploaded, then deleting the sample from the data dictionary
        temp = data.copy()
        for sample_name, replicates in temp.items():
            if len(replicates) == 0:
                data.pop(sample_name)

    # Samples and file uploaded
    print("Samples and files uploaded", data)

    if len(data) == 0:
        return render_template("initialiser.html", initialiser=True)

#     # Storing in classes
#     # sample1 = Sample(sample_one_name)
#     # sample2 = Sample(sample_two_name)

#     # Or in variables
#     sample1 = {}
#     sample2 = {}

    mhcclass = request.form.get('MHCClass')
    print('MHC Class of Interest', mhcclass)

    alleles_unformatted = request.form.get('alleles')

    if (mhcclass == 'mhc2'):
        valid_alleles_present, message = croos_check_the_allele(alleles_unformatted, os.path.join('app', 'references data', "MHC_Class2_allele_names_unformatted.txt"))
    else :
        valid_alleles_present, message = croos_check_the_allele(alleles_unformatted, os.path.join('app', 'references data', "MHC_allele_names_unformatted.txt"))

    if alleles_unformatted != "" and not valid_alleles_present:
        return "Valid alleles not passed for the job."

    # saving mhc class selected in a file
    mhcclass_selected_file = open(os.path.join('app', 'static', 'images', taskId, "mhcclass.txt"), "w")
    mhcclass_selected_file.write(mhcclass)
    mhcclass_selected_file.close()

    # saving alleles selected in a file
    allele_file = open(os.path.join('app', 'static', 'images', taskId, "selectedalleles.txt"), "w")
    allele_file.write(alleles_unformatted)
    allele_file.close()

    # Prediction tools selected by the user
    if (mhcclass == 'mhc2'):
        predictionTools = ['MixMHC2pred']
    else :
        predictionTools = request.form.getlist('predictionTools')
    print("Prediction tools selected: {}".format(predictionTools))

    # saving used prediction tools
    predictiontools = open(os.path.join('app', 'static', 'images', taskId, "predictiontools.txt"), "w")
    predictiontools.write(','.join(predictionTools))
    predictiontools.close()

    # Creating directories to store binding prediction results
    for sample, replicates in data.items():
        for predictionTool in predictionTools:
            for replicate in replicates:
                if alleles_unformatted != "":
                    for allele in alleles_unformatted.split(','):
                        try:
                            if sample != 'Control':
                                path = os.path.join('app', 'static', 'images', taskId, sample, predictionTool, replicate[:-4], 'binders',allele)
                            else:
                                path = os.path.join('app', 'static', 'images', taskId, sample)
                            if not os.path.exists(path):
                                # os.makedirs(directory)
                                Path(path).mkdir(parents=True, exist_ok=True)
                                print("Directory Created : {}".format(path))
                        except FileExistsError:
                            print("Directory already exists {}".format(path))

    # Converting alleles from A0203 format to HLA-A02:03
    alleles = ''
    if alleles_unformatted != "":
        temp = list()
        for allele in alleles_unformatted.split(','):

            # Appedning in genral format, if mhc class 1 is of interest
            if mhcclass == 'mhc1':
                temp.append('HLA-{}:{}'.format(allele[:3],allele[3:]))
            else:
                temp.append(allele) 

        alleles = ",".join(temp)
        del temp
                
    sample_data = {}
    # control_data = {}
    
    # Loading sample data in pandas frames
    for sample_name, file_names in data.items():

        sample_data[sample_name] = dict()
        for replicate in file_names:
            sample_data[sample_name][replicate] = pd.read_csv(os.path.join(dirName, sample_name, replicate))

    # Loading control data in pandas frames
    # for control_replicate in control:
        # control_data[control_replicate] = pd.read_csv(os.path.join(dirName, "Control", control_replicate))

    # Have to later add the user input for length
    for sample_name, sample in sample_data.items():
        sample_data[sample_name] = filterPeaksFile(sample, minLen=minLen, maxLen=maxLen)


    bar_percent = plot_lenght_distribution(sample_data, hist='percent')
    bar_density = plot_lenght_distribution(sample_data, hist='density')

    # Saving 8 to 14 nmers for mhc1 predictions or 12 to 20 for mhc2 predictions
    if mhcclass == 'mhc1':
        minLenForPrediction = 8
        maxLenForPrediction = 14
    elif mhcclass == 'mhc2':
        minLenForPrediction = 12
        maxLenForPrediction = 20

    saveNmerData(dirName, sample_data, peptideLength=[minLenForPrediction,maxLenForPrediction], unique = True)

    for i in range(minLenForPrediction,maxLenForPrediction+1):
        saveNmerData(dirName, sample_data, peptideLength=i, unique = True)

    # Saving 9mers in both mhc 1 and mhc 2 class of analysis to be used for seqlogo and gibbs plot
    saveNmerData(dirName, sample_data, peptideLength=9, unique=True)
   
    # Generating binding predictions
    if alleles!="":    
        for predictionTool in predictionTools:
            generateBindingPredictions(taskId, alleles, predictionTool)

    # Fetching the binders from the results
    if alleles!="":    
        for predictionTool in predictionTools:
            saveBindersData(taskId, alleles_unformatted, predictionTool, mhcclass)


    predicted_binders = None
    # Method to get the prediction results
    if alleles!="":    
        predicted_binders = getPredictionResuslts(taskId,alleles_unformatted,predictionTools,sample_data.keys())
    
    upsetLayout = getPredictionResusltsForUpset(taskId,alleles_unformatted,predictionTools,sample_data.keys())

    # Calling script to generate sequence logos
    # subprocess.call('sudo python3 {} {} {}'.format(os.path.join('app','seqlogo.py'), taskId, data_mount), shell=True)

    # # Method to return names of png files of seqlogos
    # # This value is supposed to be returned from saveNmerDate method but for now writting
    # # temporary script to return names of seqlogos pngs files in a dictionary.
    
    # Generating Seq2Logos and GibbsCluster only if MHC 2 class of interest
    # Calling script to generate sequence logos
    subprocess.check_call(['python3', os.path.join('app','seqlogo.py'), taskId, data_mount], shell=False)

    seqlogos = getSeqLogosImages(sample_data)
    # seqlogos = {}

    # Calling script to generate gibbsclusters
    subprocess.check_call(['python3', os.path.join('app', 'gibbscluster.py'), taskId, data_mount, mhcclass], shell=False)

    # Getting names of the gibbscluster
    gibbsImages = getGibbsImages(taskId, sample_data)
    # gibbsImages = {}

    showSeqLogoandGibbsSection = True
    
    # Do show Majority Voted option when MHC Class 2 analysis
    if mhcclass == 'mhc2':
        hideMajorityVotedOption = False
    else :
        hideMajorityVotedOption = True

    # Data required to plot upset plot to show peptides overlap
    overlapLayout = {}
    overlapLayout = getOverLapData(sample_data)

    return render_template('analytics.html',overlapLayout=overlapLayout, taskId=taskId, peptide_percent=bar_percent, peptide_density=bar_density, seqlogos = seqlogos, gibbsImages = gibbsImages, analytics=True,predicted_binders=predicted_binders, predictionTools = predictionTools,upsetLayout=upsetLayout, showSeqLogoandGibbsSection=showSeqLogoandGibbsSection, hideMajorityVotedOption=hideMajorityVotedOption)

@app.route("/analytics")
def analytics():

    return render_template("error.html",analytics=True, msg = 'initialiser')

@app.route('/<uuid:taskId>')
def getExistingReport(taskId):

    global DEMO_TASK_ID
    demo = False
    # Static ID for the demo
    if str(taskId) == DEMO_TASK_ID:
        demo = True
        pass
    elif is_valid_uuid(taskId) == False:
        return f"The given ID is not a valid task ID."

    taskId = str(taskId)

    # Confirming the project root is correct
    os.chdir(project_root)

    # predictionTools = ['MixMHCpred','ANTHEM','NetMHCpan']
    with open(os.path.join('app', 'static', 'images', taskId, "predictiontools.txt")) as f:
        predictionTools = f.readline().split(',')

    # MHC Class of Interest
    with open(os.path.join('app', 'static', 'images', taskId, "mhcclass.txt")) as f:
        mhcclass = f.readline()

    data = {}
    maxLen = 30
    minLen = 5
    sample_data = {}
    dirName = os.path.join(data_mount, taskId)
    predicted_binders = None
    with open(os.path.join('app', 'static', 'images', taskId, "selectedalleles.txt")) as f:
        alleles_unformatted = f.readline() 

    samples =[ f.name for f in os.scandir(dirName) if f. is_dir()]

    # Saving the data and loading into the dictionary
    for sample_name  in samples:

        # Not including the control group in data dict 
        # if sample_name != "Control":
        data[sample_name] = list()

        filenames = os.listdir(os.path.join(dirName,sample_name))
        replicates = [ filename for filename in filenames if filename.endswith( ".csv" ) ]

        for file_filename in replicates:
            data[sample_name].append(file_filename)
            
        # If control data is not uploaded, then deleting the sample from the data dictionary
        temp = data.copy()
        for sample_name, replicates in temp.items():
            if len(replicates) == 0:
                data.pop(sample_name)


    # Loading sample data in pandas frames
    for sample_name, file_names in data.items():
        sample_data[sample_name] = dict()
        for replicate in file_names:
            sample_data[sample_name][replicate] = pd.read_csv(os.path.join(dirName, sample_name, replicate))

    # Loading control data in pandas frames
    # for control_replicate in control:
        # control_data[control_replicate] = pd.read_csv(os.path.join(dirName, "Control", control_replicate))

    # Have to later add the user input for length
    for sample_name, sample in sample_data.items():
        sample_data[sample_name] = filterPeaksFile(sample, minLen=minLen, maxLen=maxLen)


    bar_percent = plot_lenght_distribution(sample_data, hist='percent')
    bar_density = plot_lenght_distribution(sample_data, hist='density')
    
    seqlogos = getSeqLogosImages(sample_data)
    gibbsImages = getGibbsImages(taskId, sample_data)
    # seqlogos = {}
    # gibbsImages = {}

    showSeqLogoandGibbsSection = True
    
    # Do show Majority Voted option when MHC Class 2 analysis
    if mhcclass == 'mhc2':
        hideMajorityVotedOption = False
    else :
        hideMajorityVotedOption = True

    if alleles_unformatted != '':
        predicted_binders = getPredictionResuslts(taskId,alleles_unformatted,predictionTools,sample_data.keys())

    upsetLayout = getPredictionResusltsForUpset(taskId,alleles_unformatted,predictionTools,sample_data.keys())

    # Data required to plot upset plot to show peptides overlap
    overlapLayout = {}
    overlapLayout = getOverLapData(sample_data)

    return render_template('analytics.html', overlapLayout=overlapLayout, taskId=taskId, analytics=True, demo=demo, peptide_percent=bar_percent, peptide_density=bar_density, seqlogos =seqlogos, gibbsImages=gibbsImages, upsetLayout=upsetLayout, predicted_binders=predicted_binders,predictionTools=predictionTools, showSeqLogoandGibbsSection=showSeqLogoandGibbsSection, hideMajorityVotedOption=hideMajorityVotedOption)

# Method to manage experiment ID
def getTaskId():
    # Generate a random UUID
    unique_id = uuid.uuid4()

    # Convert UUID to string
    unique_id_str = str(unique_id)

    return unique_id_str

def is_valid_uuid(submission_id):
    try:
        uuid.UUID(str(submission_id))  # always cast to string first
        return True
    except ValueError:
        return False

# This method is to create the bar graphs for an input file not created already
@app.route("/api/generateGibbs", methods=["POST"])
def createGibbsBar():
    
    cluster = request.form['cluster']
    taskId = request.form['taskId']
    if is_valid_uuid(taskId) == False:
        return f"The ID '{taskId}' is not a valid task ID."
    replicate = request.form['replicate']
    sample = request.form['sample']

    print(f'generateGibbs : Passed params : Cluster={cluster}, taskId={taskId}, replicate={replicate}, sample={sample}')

    # MHC Class of Interest
    with open(os.path.join('app', 'static', 'images', taskId, "mhcclass.txt")) as f:
        mhcclass = f.readline()

    barLocation = glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}/images/gibbs.KLDvsCluster.barplot.JPG')

    if len(barLocation) ==1:
        barLocation = barLocation[0][4:]

    else:

        # Deleting previous meta files present
        dirpath = f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}'
        if os.path.exists(dirpath) and os.path.isdir(dirpath):
            shutil.rmtree(dirpath)

        subprocess.check_call(['python3', os.path.join('app', 'gibbsclusterBarGraph.py'), taskId, data_mount, sample, replicate], shell=False)

        barLocation = glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}/images/gibbs.KLDvsCluster.barplot.JPG')

        if len(barLocation) ==1:
            barLocation = barLocation[0][4:]
        
        else: 
            barLocation = f'/static/others/gibbsBarNotFound.JPG'
        
    if len(cluster) == 0:
        bestCluster = pd.read_table(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}/images/gibbs.KLDvsClusters.tab')
        bestCluster = bestCluster[bestCluster.columns].sum(axis=1).idxmax()

        print(f"generateGibbs : Best Cluster for {sample}'s {replicate} : {bestCluster}")

        seqClusters = [[x[4:],"Could not be calculated"] for x in sorted(glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}/logos/gibbs_logos_*of{bestCluster}*.jpg'))]

    else:
        seqClusters = [[x[4:], "Could not be calculated"] for x in sorted(glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}/logos/gibbs_logos_*of{cluster}*.jpg'))]

        if len(seqClusters) != int(cluster):
            subprocess.check_call(['python3', os.path.join('app', 'gibbsclusterSeqLogo.py'), taskId, data_mount, sample, replicate, cluster, mhcclass], shell=False)
            seqClusters = [[x[4:], "Could not be calculated"] for x in sorted(glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}/logos/gibbs_logos_*of{cluster}*.jpg'))]

    # Adding information regarding number of peptides in the core
    findNumberOfPeptidesInCore(seqClusters, taskId, sample, replicate+'.txt')

    return {barLocation:seqClusters}

# This method is used to get the found binders in different combinations
@app.route("/api/getBinders", methods=["POST"])
def getBinders():

    tool = request.form['tool']
    taskId = request.form['taskId']
    if is_valid_uuid(taskId) == False:
        return f"The ID '{taskId}' is not a valid task ID."
    allele = request.form['allele']
    listonly = request.form['list']
    replicates = request.form['replicates']

    print(f'getBinder Post request: tool={tool}, taskId={taskId}, allele={allele}, listonly={listonly}, replicates={replicates}')

    predictionTools = ['ANTHEM','NetMHCpan','MixMHCpred']

    samples = {}


    replicates = replicates.split(',')
    
    for i in replicates:
        
        entries = i.split(';')

        if entries[0] == allele:

            if entries[1] not in samples.keys():
                samples[entries[1]] = []
                samples[entries[1]].append(entries[2])
            
            else:
                samples[entries[1]].append(entries[2])

    print('getBinder Post request: sample structure : ', samples)

    if listonly == "":

        return 'bindersFile'
    
    else:
        res = []
        for sample,replicates in samples.items():
            if tool !="":
                binder_files = []
            else:
                binder_files = {}
                for i in predictionTools:
                    binder_files[i] = []

            for replicate in replicates:
                
                if tool == "":

                    binding = getPredictionResuslts(alleles=allele,taskId=taskId,methods=predictionTools,samples=[sample])

                    for method in predictionTools:

                        try:
                            binder_files[method].append(binding[sample][allele][method][replicate])
                        except KeyError:
                            continue

                else:
                    binding = getPredictionResuslts(alleles=allele,taskId=taskId,methods=[tool],samples=[sample])
                    try:
                        binder_files.append(binding[sample][allele][tool][replicate])
                    except KeyError:
                        continue

            binders = []

            print('getBinder Post request: Binder Files : ', binder_files)

            if tool == "":
                binders = {}

                # i is tool and j is file location
                for i,j in binder_files.items():
                    binders[i] = []

                    for k in j:
                        binders[i].extend(pd.read_csv(os.path.join('app',k))['Peptide'].to_list())    

                first = set(binders[predictionTools[0]]).intersection(set(binders[predictionTools[1]]))
                second = set(binders[predictionTools[1]]).intersection(set(binders[predictionTools[2]]))
                third = set(binders[predictionTools[0]]).intersection(set(binders[predictionTools[2]]))

                binders = first.union(second).union(third)

            elif tool=="MixMHC2pred":	
                for i in binder_files:	
                    binders.extend(pd.read_csv(os.path.join('app',i))['Peptides : PlainPeptide : Core_best'].dropna().to_list())	
                binders = set(binders)
                
            else:
                for i in binder_files:
                    binders.extend(pd.read_csv(os.path.join('app',i))['Peptide'].to_list())
                binders = set(binders)

            res.append({'name':sample,'elems':list(binders)})
    
    return jsonify(res)

@app.route("/api/getOverlapPeptides", methods=["POST"])
def getOverLapPeptides():

    taskId = request.form['taskId']
    if is_valid_uuid(taskId) == False:
        return f"The ID '{taskId}' is not a valid task ID."
    replicates = request.form['replicates']

    res = []
    replicates = replicates.split(',')

    for sample in os.listdir('{}/{}'.format(data_mount,taskId)):

        peptides = set()

        for replicate in os.listdir('{}/{}/{}'.format(data_mount,taskId,sample)):
            if replicate[-12:] == '8to14mer.txt' or replicate[-13:]=='12to20mer.txt':

                replicate_name = ""    
                # Original upload file used to derive all other columns present in the input file
                if replicate[-12:] == '8to14mer.txt':
                    replicate_name = replicate[:-13]
                elif replicate[-13:]=='12to20mer.txt':
                    replicate_name = replicate[:-14]

                if replicate_name != "":
                    
                    for i in replicates:
                        if i.split(';')[0] == sample and i.split(';')[1] == replicate_name:
                            peptides.update(pd.read_csv('{}/{}/{}/{}'.format(data_mount,taskId,sample,replicate), header=None)[0].to_list())
                            break

        res.append({'name':sample,'elems':list(peptides)})
        
    return jsonify(res)

@app.route("/api/getSeqLogo", methods=["POST"])

def getSeqLogo():

    name = request.form['name']

    # removing '(' and ')' from name if present (causing problem for subprocess call)
    print(name)
    name = str(name).replace('∩','and').replace('(','').replace(')','').replace(' ','_').strip()


    taskId = request.form['taskId']    
    if is_valid_uuid(taskId) == False:
        return f"The ID '{taskId}' is not a valid task ID."
    elems = request.form['elems']

    # print(type(elems))
    # peptides = json.loads(elems)
    peptides = elems.split(',')
    peptides_location_forseqlogo = os.path.join(project_root,'app','static','images',taskId,'selected-9mer-binders-for-seqlogo.txt')
    binders_location = os.path.join(project_root,'app','static','images',taskId,'selected-binders.txt')

    peptides = pd.DataFrame(peptides)
    peptides.columns = ['peptide']

    # In case of mhc 2 class, saving 2 columns: Peptide and Binding core
    if peptides.shape[0] > 0 and peptides[peptides['peptide'].str.contains(':')].shape[0]>0:
        total_peptides = peptides.shape[0]
        peptideswithcores = peptides['peptide'].str.split(' : ',expand=True)
        peptideswithcores.columns = ['Peptide' ,'PlainPeptide','Core']

        # Saving all binders which can be downloaded in text files
        peptideswithcores[['Peptide','Core']].to_csv(binders_location,index=False)

        nine_mers = peptideswithcores.drop_duplicates(subset='Core').shape[0]
        # Saving 9mer binder cores to be used for seq2logo generation
        peptideswithcores[['Core']].drop_duplicates(subset='Core').to_csv(peptides_location_forseqlogo,index=False,header=False)

    # else using only peptide column
    else:
        total_peptides = peptides.shape[0]
        # Saving all binders which can be downloaded in text files
        peptides.to_csv(binders_location,index=False,header=False)
        peptides =  peptides[peptides.peptide.apply(lambda x: True if len(x)==9 else False)]
        nine_mers = peptides.shape[0]
        # Saving 9mer binders to be used for seq2logo generation
        peptides.to_csv(peptides_location_forseqlogo,index=False,header=False)

    seqLogoLocation = os.path.join(project_root,'app','static','images',taskId,'seqLogoApi')

    print('python3 {} {} {} {} {} {}'.format(os.path.join('app','seqlogoAPI.py'),peptides_location_forseqlogo,seqLogoLocation,name,nine_mers,total_peptides))

    subprocess.check_call(['python3', os.path.join('app','seqlogoAPI.py'),peptides_location_forseqlogo,seqLogoLocation,name,str(nine_mers),str(total_peptides)], shell=False)

    return os.path.join('static','images',taskId,'seqLogoApi-001.jpg')

@app.route("/help")
def help():
    
    # Checking to platform, if it is windows, wsl will be initialised
    #if request.user_agent.platform =='windows':
        #setUpWsl()
    return render_template("help.html", help=True)

@app.route("/pepscanner", methods=["GET"])
def pepscanner():

    return render_template("pepscanner.html", pepscanner=True,pep=False)

@app.route("/api/pepscanner", methods=["POST"])
def generatePepscanner(demo=False):

    taskId = getTaskId()

    dirName = os.path.join('app', 'static', 'images', taskId,'protienandepeptides')
    try:
        # Create target Directory
        os.makedirs(dirName)
        print("Directory " , dirName ,  " Created ") 
    except FileExistsError:
        print("Directory " , dirName ,  " already exists")

    # Deleing any existing heatmap
    myfile=os.path.join(project_root,'app/static/images/pepscanner.png')

    ## If file exists, delete it ##
    if os.path.isfile(myfile):
        os.remove(myfile)
    else:    ## Show an error ##
        print("Error: %s file not found" % myfile)

    demo = request.form.get('demo')

    print(demo == "true")
    if (demo== "true"):
        peptides = 'Type proteins of interest (optional; human proteins only)'  
        fileName = 'elutiondata.csv'

        # Input peptide file
        peptides_file = os.path.join(project_root,'app','static',fileName)

    else: 

        # Extracing passed file and the peptides list
        uploaded_file = request.files['file']
        peptides = request.form['peptides']    
        fileName = uploaded_file.filename.replace('C:\\fakepath\\',"")

        # Input peptide file
        peptides_file = os.path.join(project_root,'app','static','images',taskId,fileName)
        
        if uploaded_file.filename != '':
            uploaded_file.save(peptides_file)

    scanner = PepScan()

    # Attaching reference proteome file
    ref_proteome = os.path.join(project_root,'app','references data','uniprot-proteome_UP000005640.fasta')

    # Ignore this file for now
    accessionsidfile = os.path.join(project_root,'app','references data','accessionids.csv')

    inputFile = pd.read_csv(peptides_file)
    metadata = findMostOccuringAccessionIds(inputFile, taskId, fileName)

    scanner.search_proteome(peptide_file=peptides_file, proteome_file=ref_proteome,accessionsids=accessionsidfile)

    if peptides != 'Type proteins of interest (optional; human proteins only)':
        # Getting the list of peptides entered (and preprocessing, e.g., removing empty strings)
        peptides = peptides.replace(' ','').split(',')
        while("" in peptides) :
            peptides.remove("")
    else:
        peptides = list(metadata['top_protiens'].keys())

        if len(peptides) > 5:
            peptides = peptides[:5]

    print('Peptides passed for pepscanner: {}'.format(peptides))

    scanner.peptide_dist(peptides, taskId)

    metadata['taskId'] = taskId
    metadata['fileName'] = fileName

    return jsonify(metadata)

def findMostOccuringAccessionIds(inputFile, taskId, inputFileName):
    
    accessions = inputFile['Accession'].to_list()
    accessionIds = []
    metadata = {}
    
    for i in accessions:
        for j in str(i).split(':'):
            found = re.search(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}",j)

            if found:
                accessionIds.append(found[0])
    

    metadata['unique_peptides'] = len(set(accessionIds))
    metadata['top_protiens'] = {}

    if metadata['unique_peptides'] >= 10:
        accessionIds = Counter(accessionIds).most_common(10)
    else:
        accessionIds = Counter(accessionIds).most_common(metadata['unique_peptides'])

    # Reading mapping file
    mapping = pd.read_csv(os.path.join(project_root,'app','references data','proteinmapping.csv'))
        
    odict = OrderedDict(accessionIds)

    for key, value in odict.items():
        
        if (len (mapping[mapping['Protein'] == key])>0):
            gn = mapping.loc[mapping['Protein'] == key, "GN"].iloc[0]
            species = mapping.loc[mapping['Protein'] == key, "Species"].iloc[0]
        else:
            gn = ""
            species = ""

        odict[key]= [value, gn, species]

    metadata['top_protiens'] = odict

    for accessiondId in metadata['top_protiens'].keys():
        fileName = inputFileName+ ' ' + accessiondId + '.csv'
        subFile = inputFile[inputFile['Accession'].str.contains(accessiondId, na=False)]
        subFile.to_csv(os.path.join(project_root,'app','static','images',taskId,'protienandepeptides',fileName))

    return metadata

# Following method return the pre-run job using the specified task id.
@app.route("/demo")
def demo():
    global DEMO_TASK_ID
    return getExistingReport(DEMO_TASK_ID)

def validate_sample_name(input_text):
    # Check if the input is not null
    if not input_text:
        return False, "Name cannot be null."

    # Check if the input is not more than 30 characters
    if len(input_text) > 30:
        return False, "Nmae must not exceed 30 characters."

    # Check if the input contains only alphanumeric characters
    if not re.match("^[a-zA-Z0-9_]+$", input_text):
        return False, "Name must contain only alphanumeric characters."

    # If all checks pass, the input is valid
    return True, "Name is valid."

def croos_check_the_allele(items, file_path):
    with open(file_path, 'r') as file:
        # Read lines from the file
        file_content = file.readlines()

        # Remove newline characters from each line
        file_content = [line.strip() for line in file_content]

        # Split the input items into a list
        input_items = items.split(',')

        # Check if each item is present in the file
        for item in input_items:
            if item not in file_content:
                return False, f"Allele '{item}' not found in the file."

    return True, "All alleles are present in the file."