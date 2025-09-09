from app import app, celery
from flask import render_template, request, jsonify, redirect, url_for, send_file, make_response, abort
from werkzeug.utils import secure_filename
from werkzeug.exceptions import NotFound
from app.utils import *
from pathlib import Path
from app.Pepscan import PepScan
from collections import Counter,OrderedDict
import uuid, logging, base64, re, shutil, glob, os, pandas as pd, subprocess, io, requests, zipfile, json, smtplib, datetime, traceback
from Bio import SeqIO
from constants import *
from email.mime.text import MIMEText
from app.email_registry import save_email, get_email
from app.job_registry import insert_job, update_job_status

project_root = os.path.dirname(os.path.realpath(os.path.join(__file__, "..")))

# DEMO Task ID
DEMO_TASK_ID = app.config['DEMO_TASK_ID']

data_mount = app.config['IMMUNOLYSER_DATA']
logger = logging.getLogger(__name__)
# Configure logging format and level as needed
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

# Load Allele dictionary
ALLELE_DICTIONARY = pd.read_csv(os.path.join(project_root,'app','static','Immunolyser2.0_Allele_Dictionary.csv'))

@app.route('/submit_email/<job_id>', methods=['POST'])
def submit_email(job_id):
    data = request.get_json()
    email = data.get('email')

    if not email or '@' not in email:
        logging.warning(f"Invalid email submitted for job {job_id}: {email}")
        return "Invalid email", 400

    save_email(job_id, email)
    logging.info(f"Email {email} registered for job {job_id}")
    return "Email registered", 200

def send_email(to_email, job_id, success=True, error_msg=None):

    from_email = os.getenv("EMAIL_ADDRESS")
    password = os.getenv("EMAIL_APP_PASSWORD")

    if not from_email or not password:
        print("Email credentials are not set in environment variables.")
        return

    if success:
        subject = f"Immunolyser job {job_id} Completed"
        body = f"Your job {job_id} has completed successfully!\nAccess the results here: https://dev.immunolyser.cloud.edu.au/{job_id}"
    else:
        subject = f"Immunolyser job {job_id} Failed"
        body = f"Your job {job_id} failed with the following error:\n{error_msg}"

    msg = MIMEText(body)
    msg['Subject'] = subject
    msg['From'] = from_email
    msg['To'] = to_email

    try:
        print("Connecting to SMTP server...")
        server = smtplib.SMTP_SSL('smtp.gmail.com', 465)
        print("Logging in...")
        server.login(from_email, password)
        print("Sending email...")
        server.sendmail(from_email, [to_email], msg.as_string())
        print("Email sent.")
        server.quit()
    except Exception as e:
        print(f"Failed to send email: {e}")

@app.route("/initialiser", methods=["POST", "GET"])
def initialiser():

    if request.method == 'GET':
        # Handle GET request
        return render_template("initialiser.html", 
                                initialiser=True,
                                sample_name_max_length=app.config['SAMPLE_NAME_MAX_LENGTH'],
                                max_samples=app.config['MAX_SAMPLES'],
                                max_total_peptides=app.config['MAX_TOTAL_PEPTIDES'],
                                max_alleles=app.config['MAX_ALLELES'])
    
    elif request.method == 'POST':
        # Handle POST request

        # Create list of sample names and the files information
        samples = {}
        files_info = {}
        for key, file_list  in request.files.items():

            # Processing content first
            replicates_object = {}
            replicates = request.files.getlist(key)
            for replicate in replicates:
                file_filename = secure_filename(replicate.filename)
                file_content = replicate.read()
                file_content_base64 = base64.b64encode(file_content).decode('utf-8')
                replicates_object[file_filename] = file_content_base64

            samples[request.form[key]] = replicates_object

        # Raise error if no samples uploaded
        if len(samples) == 0:
            return f"No Sample uploaded!"

    #         filename = secure_filename(file.filename)
    # file_content = file.read()  # R

        motif_length = request.form.get('motif_length')
        mhcclass = request.form.get('mhc_class')
        alleles_unformatted = request.form.get('alleles')
        species = request.form.get('species')
        use_mhc_tp_full_DB = request.form.get('useFullDB', 'no')
        
        ip_address = request.headers.get('X-Forwarded-For', request.remote_addr)
        user_agent = request.headers.get('User-Agent')
        referrer = request.referrer

        # Prediction tools selected by the user
        if mhcclass == MHC_Class.Two:
            predictionTools = [
                Class_Two_Predictors.MixMHC2pred.to_dict(),
                Class_Two_Predictors.NetMHCpanII.to_dict(),
            ]
        else:
            predictionTools = [
                Class_One_Predictors.NetMHCpan.to_dict(),
                Class_One_Predictors.MixMHCpred.to_dict(),
                Class_One_Predictors.MHCflurry.to_dict(),
            ]

        task = submit_job.delay(samples, motif_length, mhcclass, alleles_unformatted, predictionTools, species, use_mhc_tp_full_DB)

        insert_job(
            job_id=task.id,
            ip_address=ip_address,
            mhc_class=mhcclass,
            species=species,
            alleles=alleles_unformatted,
            user_agent=user_agent,
            referrer=referrer,
            status="SUBMITTED"
        )

        return redirect(url_for('job_confirmation', task_id=task.id))

@celery.task(name='app.submit_job', bind=True)
def submit_job(self, samples, motif_length, mhcclass, alleles_unformatted, predictionTools, species, use_mhc_tp_full_DB):    

    try:
        # Have to take this input from user
        maxLen = 30
        minLen = 5
        logger.info('Preferred Motif Length: %s', motif_length)
        logger.info('MHC Class of Interest: %s', mhcclass)
        logger.info('alleles_unformatted: %s', alleles_unformatted)
        logger.info('Species: %s', species)
        logger.info('Use Full Database Search: %s', use_mhc_tp_full_DB)

        # Deserialize `predictionTools`
        predictionTools = [Predictor.from_dict(tool) for tool in predictionTools]
        logger.info('Prediction tools selected: : %s', predictionTools)

        total_peptides = 0
        max_rows = app.config['MAX_TOTAL_PEPTIDES']

        taskId = self.request.id

        dirName = os.path.join(data_mount, taskId)
        try:
            # Create target Directory
            os.makedirs(dirName)
            logger.info("Directory %s Created", dirName) 
        except FileExistsError:
            logger.info("Directory %s already exists", dirName)


        data = {}
        control = list()

        # Creating folders to store images
        for sample_name in samples.keys():

            is_valid, message = validate_sample_name(sample_name)

            if not is_valid:
                logger.info("Sample name is not valid. %s", message)

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
                    logger.info("Directory Created: %s", path_for_logos) 
            except FileExistsError:
                logger.info("Directory already exists: %s", path_for_logos)    


        # Saving the data and loading into the dictionary
        for sample_name, replicates in samples.items():

            # Creating sub directories to store sample data
            try:
                os.mkdir(os.path.join(dirName, sample_name))
                logger.info("Directory Created: %s", os.path.join(dirName, sample_name)) 
            except FileExistsError:
                logger.info("Directory already exists: %s", os.path.join(dirName, sample_name))

            # Not including the control group in data dict 
            # if sample_name != "Control":
            data[sample_name] = list()

            files_to_save = {}
            # First pass: Accumulate row counts
            for sample_name, replicates in samples.items():
                files_to_save[sample_name] = {}

                for file_filename, file_content_base64 in replicates.items():
                    replicate = base64.b64decode(file_content_base64)

                    if file_filename != "":
                        # Read the CSV content using pandas
                        df = pd.read_csv(io.BytesIO(replicate))

                        # Count rows excluding the header
                        row_count = len(df)
                        total_peptides += row_count

                        # Store file information for potential saving
                        files_to_save[sample_name][file_filename] = replicate

            # Check if total peptides exceed the limit
            if total_peptides > max_rows:
                raise Exception(f"Total peptides {total_peptides} exceed the maximum allowed {max_rows}.")

            # Second pass: Save files if within the limit
            for sample_name, replicates in files_to_save.items():
                # Create directories if they do not exist
                try:
                    os.mkdir(os.path.join(dirName, sample_name))
                    logger.info("Directory Created: %s", os.path.join(dirName, sample_name))
                except FileExistsError:
                    logger.info("Directory already exists: %s", os.path.join(dirName, sample_name))

                data[sample_name] = list()
                for file_filename, replicate in replicates.items():
                    # Save the file
                    with open(os.path.join(dirName, sample_name, file_filename), 'wb') as f:
                        f.write(replicate)

                    # Storing the filename in data dictionary
                    data[sample_name].append(file_filename)

            # If control data is not uploaded, then deleting the sample from the data dictionary
            temp = data.copy()
            for sample_name, replicates in temp.items():
                if len(replicates) == 0:
                    data.pop(sample_name)

        # Samples and file uploaded
        logger.info("Samples and files uploaded: %s", data)

        valid_alleles_present, message = cross_check_the_allele(alleles_unformatted, ALLELE_DICTIONARY)

        if alleles_unformatted != "" and not valid_alleles_present:
            raise Exception(f"Valid alleles not passed for the job.")

        # saving motif length selected in a file
        motif_length_file = open(os.path.join('app', 'static', 'images', taskId, "motif_length.txt"), "w")
        motif_length_file.write(motif_length)
        motif_length_file.close()

        # saving mhc class selected in a file
        mhcclass_selected_file = open(os.path.join('app', 'static', 'images', taskId, "mhcclass.txt"), "w")
        mhcclass_selected_file.write(mhcclass)
        mhcclass_selected_file.close()

        # Save allele compatibility matrix based on alleles and MHC class of preference selected.
        # if alleles_unformatted != "":
        # Split alleles only if alleles_unformatted is not an empty string
        alleles = alleles_unformatted.split(',') if alleles_unformatted else []

        # Convert predictionTools to a list of full names
        predictionToolNames = [tool.full_name for tool in predictionTools]

        # Create the DataFrame with rows as prediction tools and columns as alleles
        # If alleles is empty, DataFrame will have no columns
        allele_compatibility_matrix = pd.DataFrame(index=predictionToolNames, columns=alleles if alleles else [])

        # Populate the DataFrame
        for tool in predictionTools:
            for allele in alleles:
                match = ALLELE_DICTIONARY[
                    (ALLELE_DICTIONARY["Allele name standardised"] == allele) &
                    (ALLELE_DICTIONARY["Predictor"] == tool.full_name)
                ]
                allele_compatibility_matrix.at[tool.full_name, allele] = "Yes" if not match.empty else "No"

        # Save the DataFrame to a CSV file
        output_path = os.path.join('app', 'static', 'images', taskId, "allele_compatibility_matrix.csv")
        allele_compatibility_matrix.to_csv(output_path, index=True)

        # Creating directories to store binding prediction results
        for sample, replicates in data.items():
            for predictionTool in predictionTools:
                for replicate in replicates:
                    if alleles_unformatted != "":
                        for allele in alleles_unformatted.split(','):
                            try:
                                if sample != 'Control':

                                    # Path to store user friendly binders data
                                    path = os.path.join('app', 'static', 'images', taskId, sample, predictionTool.short_name, replicate[:-4], 'binders',allele.replace(':', '_'))

                                    # Path to store raw binder tool output
                                    path_predictor_output = os.path.join('app', 'static', 'images', taskId, sample, predictionTool.short_name, replicate[:-4],allele.replace(':', '_'))
                                else:
                                    path = os.path.join('app', 'static', 'images', taskId, sample)

                                if not os.path.exists(path):
                                    # os.makedirs(directory)
                                    Path(path).mkdir(parents=True, exist_ok=True)
                                    print("Directory Created : {}".format(path))

                                if not os.path.exists(path_predictor_output):
                                    # os.makedirs(directory)
                                    Path(path_predictor_output).mkdir(parents=True, exist_ok=True)
                                    print("Directory Created : {}".format(path_predictor_output))
                                    
                            except FileExistsError:
                                print("Directory already exists {}".format(path))
                    
        sample_data = {}
        # control_data = {}
        
        # Loading sample data in pandas frames
        for sample_name, file_names in data.items():

            sample_data[sample_name] = dict()
            for replicate in file_names:
                sample_data[sample_name][replicate] = pd.read_csv(os.path.join(dirName, sample_name, replicate))

        # Have to later add the user input for length
        for sample_name, sample in sample_data.items():
            sample_data[sample_name] = filterPeaksFile(sample, minLen=minLen, maxLen=maxLen)

        # Saving 8 to 14 nmers for class one predictions or 12 to 20 for class two predictions
        if mhcclass == MHC_Class.One:
            minLenForPrediction = 8
            maxLenForPrediction = 14
        elif mhcclass == MHC_Class.Two:
            minLenForPrediction = 12
            maxLenForPrediction = 20

        saveNmerData(dirName, sample_data, peptideLength=[minLenForPrediction,maxLenForPrediction], unique = True)

        for i in range(minLenForPrediction,maxLenForPrediction+1):
            saveNmerData(dirName, sample_data, peptideLength=i, unique = True)
    
        # Generating binding predictions
        if alleles_unformatted!="":    
            for predictionTool in predictionTools:
                generateBindingPredictions(taskId, alleles_unformatted, predictionTool, ALLELE_DICTIONARY)

        # Fetching the binders from the results
        if alleles_unformatted!="":    
            for predictionTool in predictionTools:
                saveBindersData(taskId, alleles_unformatted, predictionTool, mhcclass)

            # Store majority voting results
            # Calling method to generate csv file with Majority Voted binders
            saveMajorityVotedBinders(taskId, data, predictionTools, alleles_unformatted, ALLELE_DICTIONARY)
        

        # Do not generate Seq2Logo for Class II, if not Allele is selected
        if mhcclass == MHC_Class.One or (mhcclass == MHC_Class.Two and alleles_unformatted != ''):
            # Calling script to generate sequence logos
            subprocess.check_call(['python3', os.path.join('app','seqlogo.py'), taskId, data_mount, motif_length], shell=False)

        # Calling script to generate gibbsclusters
        subprocess.check_call(['python3', os.path.join('app', 'gibbscluster.py'), taskId, data_mount, mhcclass, motif_length], shell=False)

        # Run HLA-Clust to generate heatmap
        if mhcclass == MHC_Class.One:
            runHLAClust(
                taskId,
                data,
                species=species,
                use_mhc_tp_full_DB=use_mhc_tp_full_DB,
                logger=logger
            )

        # On job success: update status first
        update_job_status(job_id=taskId, status='SUCCESS', error_message=None, logger=logger)
        logging.info(f"Job {taskId} status updated to SUCCESS.")

        # On job success
        email = get_email(taskId)
        if email:
            logging.info(f"Found email '{email}' for completed task '{taskId}', attempting to send notification.")
            try:
                send_email(email, taskId, success=True)
                logging.info(f"Email successfully sent to {email} for job {taskId}.")
            except Exception as e:
                logging.error(f"Failed to send success email to {email} for job {taskId}: {e}")
        else:
            logging.info(f"No email found for task {taskId}; skipping email notification.")

    except Exception as main_exception:
        tb = traceback.format_exc()
        error_msg = f"{main_exception}\n{tb}"
        update_job_status(job_id=taskId, status='FAILURE', error_message=error_msg, logger=logger)
        logging.error(f"Job {taskId} failed with error: {error_msg}")

        email = get_email(taskId)
        if email:
            try:
                send_email(email, taskId, success=False)
                logging.info(f"Failure email successfully sent to {email} for job {taskId}.")
            except Exception as e:
                logging.error(f"Failed to send failure email to {email} for job {taskId}: {e}")
        else:
            logging.info(f"No email found for task {taskId}; skipping failure notification.")

        raise
    
@app.route("/analytics")
def analytics():

    return render_template("error.html",analytics=True, msg = 'initialiser')

@app.route('/job-confirmation/<task_id>')
def job_confirmation(task_id):
    message = f'Request for Immunolyser report has been received. Task ID is {task_id}'
    return render_template('onSubmission.html', message=message)

# GET method to check the status of the job. Job state is managed by Celery
@app.route('/check_status/<job_id>', methods=['GET'])
def check_status(job_id):
    job = submit_job.AsyncResult(job_id)
    if job.state == 'SUCCESS':
        return jsonify({'status': 'success'}), 200
    elif job.state == 'FAILURE':
        return jsonify({'status': 'failure', 'traceback': str(job.traceback)}), 200
    elif job.state == 'PENDING':
        return jsonify({'status': 'pending', 'traceback': str(job.traceback)}), 200
    else:
        return jsonify({'status': job.state}), 200

@app.route('/<taskId>')
def getExistingReport(taskId):

    global DEMO_TASK_ID
    demo = False
    # Static ID for the demo
    if str(taskId) == DEMO_TASK_ID:
        demo = True
        pass
    elif is_valid_uuid(taskId) == False:
        return f"The given ID is not a valid task ID."

    # Confirming the project root is correct
    os.chdir(project_root)

    # Read the allele compatibility matrix
    output_path = os.path.join('app', 'static', 'images', taskId, "allele_compatibility_matrix.csv")

    # Check if the file exists
    if not os.path.exists(output_path):
        # Raise an exception with a custom message
        return f"Due to some recent changes, the existing jobs cannot be accessed. The data is still there. If it is really important or you cannot submit another job, please email the developer with your job ID to access your job."

    allele_compatibility_matrix = pd.read_csv(output_path, index_col=0)

    # Fetch all predictors dynamically
    all_predictors = get_all_predictors()

    # Create predictionTools list by matching full names from the CSV index
    predictionTools = [
        predictor for predictor in all_predictors if predictor.full_name in allele_compatibility_matrix.index
    ]

    # MHC Class of Interest
    with open(os.path.join('app', 'static', 'images', taskId, "mhcclass.txt")) as f:
        mhcclass = f.readline()

    # Selected motif length
    with open(os.path.join('app', 'static', 'images', taskId, "motif_length.txt")) as f:
        motif_length = f.readline()

    data = {}
    maxLen = 30
    minLen = 5
    sample_data = {}
    dirName = os.path.join(data_mount, taskId)
    predicted_binders = None
    
    # Extract alleles as a comma-separated string
    alleles_unformatted = ','.join(set(allele_compatibility_matrix.columns))

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


    bar_percent = plot_lenght_distribution(sample_data, hist='percent', taskId=taskId)
    bar_density = plot_lenght_distribution(sample_data, hist='density', taskId=taskId)
    
    seqlogos = getSeqLogosImages(sample_data, task_id=taskId, motif_length=motif_length, logger=logger)
    gibbsImages = getGibbsImages(logger, taskId, sample_data)
    # seqlogos = {}
    # gibbsImages = {}

    showSeqLogoSection = True
    showGibbsSection = True
    
    # Do not show Majority Voted option when MHC Class 2 analysis
    if mhcclass == MHC_Class.Two:
        hideMajorityVotedOption = False

        if alleles_unformatted == '': # Hiding Motifs results when no alleles was selected to run Class 2 analysis
            showSeqLogoSection = False 

    else : # Need fixing as it is should be visible for Class 2
        hideMajorityVotedOption = True

    if alleles_unformatted != '':
        predicted_binders = getPredictionResuslts(taskId,alleles_unformatted,predictionTools,sample_data.keys())

    upsetLayout = getPredictionResusltsForUpset(taskId,alleles_unformatted,predictionTools,sample_data.keys())

    # Data required to plot upset plot to show peptides overlap
    overlapLayout = {}
    overlapLayout = getOverLapData(sample_data)

    # Assuming 'predictionTools' is a list of Predictor objects
    predictionTools = [tool.short_name for tool in predictionTools]

    bindingImages = getHLAClustResults(taskId, data)
    print(bindingImages)

    zip_path = zip_job_exports(taskId)
    zip_filename = os.path.basename(zip_path)

    return render_template(
        'analytics.html', 
        overlapLayout=overlapLayout, 
        taskId=taskId, 
        analytics=True, 
        demo=demo, 
        peptide_percent=bar_percent, 
        peptide_density=bar_density, 
        seqlogos=seqlogos, 
        gibbsImages=gibbsImages, 
        upsetLayout=upsetLayout, 
        predicted_binders=predicted_binders, 
        predictionTools=predictionTools,  # List of short_names here
        showSeqLogoSection=showSeqLogoSection,
        showGibbsSection=showGibbsSection, 
        hideMajorityVotedOption=hideMajorityVotedOption,
        bindingImages=bindingImages,
        zip_filename=zip_filename
    )

# Method to manage experiment ID
def getTaskId():
    # Generate a random UUID
    unique_id = uuid.uuid4()

    # Convert UUID to string
    unique_id_str = str(unique_id)

    return unique_id_str

def is_valid_uuid(submission_id):
    try:
        # Try to create a UUID object from the given string
        uuid_obj = uuid.UUID(submission_id)
        return True
    except ValueError:
        # ValueError will be raised if the string is not a valid UUID
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

    # Motif length
    with open(os.path.join('app', 'static', 'images', taskId, "motif_length.txt")) as f:
        motif_length = f.readline()

    # MHC Class of Interest
    with open(os.path.join('app', 'static', 'images', taskId, "mhcclass.txt")) as f:
        mhcclass = f.readline()

    barLocation = glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}/images/gibbsKLDvsCluster.barplot.JPG')

    if len(barLocation) ==1:
        barLocation = barLocation[0][4:]

    else:

        # Deleting previous meta files present
        dirpath = f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}'
        # if os.path.exists(dirpath) and os.path.isdir(dirpath):
        #     shutil.rmtree(dirpath)

        # subprocess.check_call(['python3', os.path.join('app', 'gibbsclusterBarGraph.py'), taskId, data_mount, sample, replicate, motif_length, mhcclass], shell=False)

        barLocation = glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}/images/gibbsKLDvsCluster.barplot.JPG')

        if len(barLocation) ==1:
            barLocation = barLocation[0][4:]
        
        else: 
            barLocation = f'/static/others/gibbsBarNotFound.JPG'
        
    if len(cluster) == 0:
        bestCluster = pd.read_table(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}/images/gibbs.KLDvsClusters.tab')
        bestCluster = bestCluster[bestCluster.columns].sum(axis=1).idxmax()

        print(f"generateGibbs : Best Cluster for {sample}'s {replicate} : {bestCluster}")

        seqClusters = [[x[4:],"Number of peptides in core could not be calculated", "Allele not predicted", "Score not calculated", "Url of reference motif"] for x in sorted(glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}/logos/gibbs_logos_*of{bestCluster}*.jpg'))]

    else:
        seqClusters = [[x[4:], "Number of peptides in core could not be calculated", "Allele not predicted", "Score not calculated", "Url of reference motif"] for x in sorted(glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}/logos/gibbs_logos_*of{cluster}*.jpg'))]

        if len(seqClusters) != int(cluster):
            # subprocess.check_call(['python3', os.path.join('app', 'gibbsclusterSeqLogo.py'), taskId, data_mount, sample, replicate, cluster, mhcclass, motif_length], shell=False)
            seqClusters = [[x[4:], "Number of peptides in core could not be calculated", "Allele not predicted", "Score not calculated", "Url of reference motif"] for x in sorted(glob.glob(f'app/static/images/{taskId}/{sample}/gibbscluster/{replicate}/logos/gibbs_logos_*of{cluster}*.jpg'))]

    # Adding information regarding number of peptides in the core
    findNumberOfPeptidesInCore(seqClusters, taskId, sample, replicate+'.txt')

    # Append predicted allele information to the object
    appendPredictedAllelesInfo(seqClusters, taskId, sample, replicate+'.txt')

    return {barLocation: seqClusters}

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

    predictionTools = ['MHCflurry','NetMHCpan','MixMHCpred']

    samples = {}


    replicates = replicates.split(',')
    
    # Build a dictionary 'samples' where keys are sample names (entries[1]) 
    # and values are lists of replicate IDs (entries[2]) that match the specified allele.
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

            # If tool is equal to an empty string: That is request for majority voted binder from the client side.
            if tool !="":
                binder_files = []
            else:
                binder_files = {}
                for i in predictionTools:
                    binder_files[i] = []

            for replicate in replicates:
                
                # If tool string is empty: Fetch resuls for the tools in predictionTools list
                if tool == "":

                    binding = getPredictionResuslts(alleles=allele,taskId=taskId,methods_passed=predictionTools,samples=[sample])

                    for method in predictionTools:

                        try:
                            # Appeding to binder_fikes dictionary
                            binder_files[method].append(binding[sample][allele][method][replicate])
                        except KeyError:
                            continue

                else:
                    binding = getPredictionResuslts(alleles=allele,taskId=taskId,methods_passed=[tool],samples=[sample])
                    try:
                        # Appeding to binder_fikes list
                        binder_files.append(binding[sample][allele][tool][replicate])
                    except KeyError:
                        continue

            binders = []

            print('getBinder Post request: Binder Files : ', binder_files)

            # If tool is not selected: User asked for majority voted results. In following if, intersections are derived across samples.
            if tool == "":
                binders = {}

                # i is tool and j is file location
                for i,j in binder_files.items():
                    binders[i] = []

                    for k in j:
                        binders[i].extend(pd.read_csv(os.path.join('app',k))['Peptide'].to_list())    


                # first: Common binder from first and second tool
                # second: Common binder from second and third tool
                # first: Common binder from first and third tool
                first = set(binders[predictionTools[0]]).intersection(set(binders[predictionTools[1]]))
                second = set(binders[predictionTools[1]]).intersection(set(binders[predictionTools[2]]))
                third = set(binders[predictionTools[0]]).intersection(set(binders[predictionTools[2]]))

                # binders is list of union of first, second and thirdl
                binders = first.union(second).union(third)

            # For class 2. Only MixMHC2pred are considered
            elif tool=="MixMHC2pred":	
                for i in binder_files:	
                    binders.extend(pd.read_csv(os.path.join('app',i))['Peptides : PlainPeptide : Core_best'].dropna().to_list())	
                binders = set(binders)
                
            # Else. For specific prediction tool.
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

        dir_path = '{}/{}/{}'.format(data_mount, taskId, sample)
        if os.listdir(dir_path):  # Check if the directory is not empty
            for replicate in os.listdir(dir_path):
                if replicate[-12:] == '8to14mer.txt' or replicate[-13:] == '12to20mer.txt':
                    replicate_name = ""    
                    
                    # Determine replicate name
                    if replicate[-12:] == '8to14mer.txt':
                        replicate_name = replicate[:-13]
                    elif replicate[-13:] == '12to20mer.txt':
                        replicate_name = replicate[:-14]

                    if replicate_name != "":
                        for i in replicates:
                            if i.split(';')[0] == sample and i.split(';')[1] == replicate_name:
                                peptides.update(
                                    pd.read_csv('{}/{}/{}/{}'.format(data_mount, taskId, sample, replicate), header=None)[0].to_list()
                                )
                                break

            # Append result only if the directory has been processed
            res.append({'name': sample, 'elems': list(peptides)})

        
    return jsonify(res)

@app.route("/api/downloadOverlapPeptides", methods=["POST"])
def download_overlap_peptides():
    task_id = request.form['taskId']
    selected_samples = request.form.getlist('samples[]')

    if not is_valid_uuid(task_id):
        return f"The ID '{task_id}' is not a valid task ID.", 400

    peptide_sets = {}

    for sample in selected_samples:
        dir_path = os.path.join(data_mount, task_id, sample)
        peptides = set()
        if os.path.exists(dir_path):
            for fname in os.listdir(dir_path):
                if fname.endswith('8to14mer.txt') or fname.endswith('12to20mer.txt'):
                    full_path = os.path.join(dir_path, fname)
                    peptides.update(pd.read_csv(full_path, header=None)[0].to_list())
        peptide_sets[sample] = peptides

    # Intersect all selected sets
    if peptide_sets:
        intersected = set.intersection(*peptide_sets.values())
    else:
        intersected = set()

    response = make_response('\n'.join(sorted(intersected)))
    response.headers["Content-Disposition"] = "attachment; filename=overlap_peptides.txt"
    response.headers["Content-Type"] = "text/plain"
    return response

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
        # Read motif length from file
        motif_length_file_path = os.path.join('app', 'static', 'images', taskId, "motif_length.txt")
        with open(motif_length_file_path, 'r') as f:
            motif_length = int(f.read().strip())

        # Filter peptides based on the motif length
        peptides = peptides[peptides.peptide.apply(lambda x: len(x) == motif_length)]
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
    run_prot_peptigram = request.form.get('runProtPeptigram', '').lower() == 'true'

    print('run_prot_peptigram:', run_prot_peptigram)

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

    if (demo== "true"):
        protiens = ''  
        fileName = 'Liver Set1 DDA.csv'
        run_prot_peptigram = True

        # Input peptide file
        peptides_file = os.path.join(project_root,'app','static',fileName)
        # Background human proteome
        ref_proteome = os.path.join(project_root,'app','references data','uniprotkb_proteome_UP000000589_2024_08_14.fasta')

    else: 

        # Extracting passed file, background, and the peptides list
        uploaded_file = request.files['file']
        uploaded_background_file = request.files.get('background')  # Use get() to handle the case where 'background' might not be present
        protiens = request.form['protein_ids']
        fileName = uploaded_file.filename.replace('C:\\fakepath\\', "")

        # Input peptide file
        peptides_file = os.path.join(project_root, 'app', 'static', 'images', taskId, fileName)

        if uploaded_background_file:
            # Background file exists
            background_file_contents = uploaded_background_file.read().decode('utf-8')
            # Validate the uploaded FASTA file
            is_valid, message = validate_fasta(background_file_contents)
            if not is_valid:
                return jsonify({"error": message}), 400
            
            # Save the file if it is valid
            background_filename = secure_filename(uploaded_background_file.filename.replace('C:\\fakepath\\', ""))
            ref_proteome = os.path.join(project_root, 'app', 'static', 'images', taskId, background_filename)
            with open(ref_proteome, 'w') as file:
                file.write(background_file_contents)
        
        else: # Else use existing human background proteome
            ref_proteome = os.path.join(project_root,'app','references data','uniprot-proteome_UP000005640.fasta')

        # Saving the input file
        if uploaded_file.filename != '':
            uploaded_file.save(peptides_file)

    scanner = PepScan()

    inputFile = pd.read_csv(peptides_file)
    metadata = findMostOccuringAccessionIds(inputFile, taskId, fileName)

    scanner.search_proteome(peptide_file=peptides_file, proteome_file=ref_proteome)

    if protiens != '':
        # Getting the list of protiens entered (and preprocessing, e.g., removing empty strings)
        protiens = protiens.replace(' ','').split(',')
        while("" in protiens) :
            protiens.remove("")
    else:
        protiens = list(metadata['top_protiens'].keys())

        if len(protiens) > 5:
            protiens = protiens[:5]

    print('Proteins passed for pepscanner: {}'.format(protiens))

    scanner.peptide_dist(protiens, taskId)

    metadata['taskId'] = taskId
    metadata['fileName'] = fileName

    # Run ProtPeptigram only if opted by the user
    if run_prot_peptigram is True:
        peptigram_images = []

        output_path = os.path.join(project_root, 'app', 'static', 'images', taskId)
        # Generate peptigram images
        generate_peptigram(
            csv_path=peptides_file,
            fasta_path=ref_proteome,
            protein_ids=','.join(protiens),
            output_dir=output_path
        )

        # Collect relative image paths like: static/images/<taskId>/<filename>
        peptigram_images = [
            os.path.join('static', 'images', taskId, os.path.basename(f))
            for f in glob.glob(os.path.join(output_path, 'prot-peptigram_*.png'))
        ]
                # Store in metadata
        metadata['peptigram'] = peptigram_images

    return jsonify(metadata)

def validate_fasta(file_contents):
    """
    Validate a FASTA file using Biopython.

    Parameters:
        file_contents (str): The contents of the uploaded FASTA file.

    Returns:
        bool: True if the file is a valid FASTA file, False otherwise.
        str: Message indicating the validation result.
    """
    try:
        # Create a StringIO object from the file contents
        file_stream = io.StringIO(file_contents)

        # Try to parse the file stream
        records = list(SeqIO.parse(file_stream, "fasta"))

        # Check if any records were found
        if not records:
            return False, "File is empty or not a valid FASTA file."

        # Additional checks (optional)
        for record in records:
            if not record.id:
                return False, f"Record {record} has no ID."
            if not record.seq:
                return False, f"Record {record} has no sequence."

        return True, "File is a valid FASTA file."

    except Exception as e:
        return False, f"Error: {str(e)}"

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
    mapping = pd.read_csv(os.path.join(project_root,'app','references data','proteinmapping_2.csv'))
        
    odict = OrderedDict(accessionIds)

    for key, value in odict.items():
        
        if (len (mapping[mapping['Accession'] == key])>0):
            gn = mapping.loc[mapping['Accession'] == key, "GN"].iloc[0]
            species = mapping.loc[mapping['Accession'] == key, "Species"].iloc[0]
            protien = mapping.loc[mapping['Accession'] == key, "Protein"].iloc[0]
        else:
            # Fetch from UniProt if not found in local file
            uniprot_url = f"https://www.uniprot.org/uniprot/{key}.txt"
            response = requests.get(uniprot_url)
            
            if response.status_code == 200:
                data = response.text
                gn = ""
                species = ""
                protien = ""

                # Parsing the UniProt text data
                for line in data.split('\n'):
                    if line.startswith("GN   Name=") and not gn:
                        gn = line.split('=')[1].split(';')[0].strip()
                    elif line.startswith("OS   ") and not species:
                        species = line[5:].strip()
                    elif line.startswith("DE   RecName: Full=") and not protein:
                        protein = line.split('=')[1].strip(';').strip()

                # If GN and Species not found in UniProt data
                if not gn:
                    gn = "Unknown"
                if not species:
                    species = "Unknown"
                if not protien:
                    protien = "Unknown"
            else:
                gn = "Unknown"
                species = "Unknown"
                protien = "Unknown"

        odict[key]= [value, gn, species, protien]

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

def cross_check_the_allele(items, allele_dict):
    """
    Check if each item is present in the 'Allele name standardised' column of the given DataFrame.

    Args:
        items (str): A comma-separated string of alleles to check.
        allele_dict (pd.DataFrame): A DataFrame containing the "Allele name standardised" column.

    Returns:
        tuple: (bool, str) indicating if all alleles are found and a relevant message.
    """
    # Ensure the required column is present
    if 'Allele name standardised' not in allele_dict.columns:
        return False, "'Allele name standardised' column not found in the DataFrame."

    # Get the list of alleles from the DataFrame
    valid_alleles = allele_dict['Allele name standardised'].tolist()

    # Split the input items into a list
    input_items = [item.strip() for item in items.split(',')]

    # Check if each item is present in the valid_alleles list
    for item in input_items:
        if item not in valid_alleles:
            return False, f"Allele '{item}' not found in the DataFrame."

    return True, "All alleles are present in the DataFrame."

@app.route('/get_species', methods=['POST'])
def get_species():
    # Get unique species (assuming 'Gene' column contains species names)
    species_list = ALLELE_DICTIONARY['Gene'].unique()
    
    # Return the species list as a JSON response
    return jsonify(list(species_list))

@app.route('/get_mhc_classes', methods=['POST'])
def get_mhc_classes():
    species = request.json['species']  # Use request.json to access JSON data
    filtered_classes = ALLELE_DICTIONARY[ALLELE_DICTIONARY['Gene'] == species]['Class'].unique()
    return jsonify(list(filtered_classes))

@app.route('/get_alleles', methods=['POST'])
def get_alleles():
    data = request.json  # Access the JSON payload
    species = data['species']
    mhc_class = data['mhc_class']
    
    # Filter the DataFrame based on species and MHC class
    filtered_alleles = ALLELE_DICTIONARY[
        (ALLELE_DICTIONARY['Gene'] == species) & (ALLELE_DICTIONARY['Class'] == mhc_class)
    ]['Allele name standardised'].unique()
    
    # Return the filtered alleles as a JSON response
    return jsonify(list(filtered_alleles))

def zip_job_exports(taskId):
    export_folder = os.path.join(
        project_root, "app", "static", "images", taskId, "export"
    )
    zip_path = os.path.join(
        project_root, "app", "static", "images", taskId, f"export_files_{taskId}.zip"
    )

    with zipfile.ZipFile(zip_path, 'w') as zipf:
        for root, _, files in os.walk(export_folder):
            for file in files:
                file_path = os.path.join(root, file)
                # Arcname should be relative to export_folder
                arcname = os.path.relpath(file_path, export_folder)
                zipf.write(file_path, arcname=arcname)

    return zip_path

@app.route("/download-peptide-zip/<taskId>/<filename>")
def download_job_data_files_zip(taskId, filename):
    # Define a safe base directory where files are stored
    safe_base_dir = os.path.join(project_root, "app", "static", "images", taskId)
    
    # Prevent directory traversal by joining the path and normalizing it
    safe_filename = os.path.normpath(filename)  # Normalize the path to prevent traversal
    full_path = os.path.join(safe_base_dir, safe_filename)

    # Ensure the file is inside the base directory
    if not full_path.startswith(safe_base_dir):
        return "Forbidden", 403

    # Check if the file exists
    if os.path.exists(full_path):
        return send_file(full_path, as_attachment=True)
    else:
        return NotFound("File not found")
    
@app.route('/download_seq2logo_peptides/<taskid>/<sample>/<replicate>')
def download_seq2logo_peptides(taskid, sample, replicate):
    logger.info(f"Download request for peptides: taskid={taskid}, sample={sample}, replicate={replicate}")

    motif_length_path = os.path.join(project_root, 'app', 'static', 'images', taskid, 'motif_length.txt')
    
    if not os.path.isfile(motif_length_path):
        logger.error(f"Motif length file not found: {motif_length_path}")
        return abort(404, description="Motif length file not found.")

    try:
        with open(motif_length_path, 'r') as f:
            motif_length = f.read().strip()
        logger.debug(f"Motif length for task {taskid}: {motif_length}")
    except Exception:
        logger.exception(f"Failed to read motif length file: {motif_length_path}")
        return abort(500, description="Error reading motif length file.")

    # FIXED HERE: look inside sample folder, not replicate subfolder
    search_dir = os.path.join(data_mount, taskid, sample)
    pattern = f'*_{motif_length}mer.txt'
    matches = glob.glob(os.path.join(search_dir, pattern))

    if not matches:
        logger.warning(f"No file matching '*_{motif_length}mer.txt' found in {search_dir}")
        return abort(404, description=f"No peptide file found for motif length {motif_length}.")

    peptides_file = matches[0]
    logger.info(f"Serving peptide file: {peptides_file}")

    try:
        return send_file(peptides_file, as_attachment=True)
    except Exception:
        logger.exception(f"Failed to send file: {peptides_file}")
        return abort(500, description="Error sending peptide file.")

@app.route('/download_gibbscluster_core/<taskid>/<sample>/<replicate>/<cluster_attempt>')
def download_gibbscluster_core(taskid, sample, replicate, cluster_attempt):
    logger.info(f"Download request for Gibbs core: taskid={taskid}, sample={sample}, replicate={replicate}, cluster={cluster_attempt}")
    
    core_dir = os.path.join(project_root, 'app', 'static', 'images', taskid, sample, 'gibbscluster', replicate, 'cores')
    pattern = f'*{cluster_attempt}*'
    matches = glob.glob(os.path.join(core_dir, pattern))

    if not matches:
        logger.warning(f"No core file found for cluster {cluster_attempt} in {core_dir}")
        return abort(404, description=f"No core file found for cluster {cluster_attempt}.")

    core_file = matches[0]
    logger.info(f"Serving core file: {core_file}")

    try:
        return send_file(core_file, as_attachment=True)
    except Exception:
        logger.exception(f"Failed to send core file: {core_file}")
        return abort(500, description="Error sending core file.")
