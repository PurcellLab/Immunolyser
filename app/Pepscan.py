import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('agg')
import seaborn as sns
import io
import os

project_root = os.path.dirname(os.path.realpath(os.path.join(__file__, "..")))


class PepScan:
    def __init__(self):
        #filename of PEAKS output csv file containing peptide information
        self.filename = None

        #filename of FASTA file containing proteome to compare with peptide information
        self.proteome = None

        #input csv file converted into a dataframe and cleaned
        self.peptide_frame = None

        #dataframe containing information on each peptide and its location in its source protein
        self.protein_frame = None
        
    def load_proteome(self, proteome_file):
        """Loads a FASTA file containing a proteome to reference when assigning source protein sequences to peptides

        Keyword Arguments:
        proteome_file -- filename of the .fasta file containing the proteome

        Return:
        """
        #store the filename of proteome
        self.proteome = proteome_file
        proteome_file = open(proteome_file, "r")
        proteome = proteome_file.read()
        proteome_file.close()
        #format and extract proteome into a dataframe
        proteome = proteome.split("\n>")
        proteome[0] = proteome[0].replace(">", "", 1)
        for i in range(len(proteome)):
            proteome[i] = proteome[i].split("\n")
            proteome[i][0] = proteome[i][0].split("|")
            proteome[i][1] = "".join(proteome[i][1:len(proteome[i])])
            proteome[i][0].append(proteome[i][1])
            proteome[i] = proteome[i][0]
        proteome_frame = pd.DataFrame(data = proteome, columns = ["Database", "ID", "Protein Name", "Sequence"])
        return proteome_frame

    def clean(self):
        """Cleans the dataframe generated from the PEAKS file stored in self.peptide_frame, formats accession and peptide values, and removes entries with NA Accession values 
        """
        clean_frame = self.peptide_frame
        clean_frame["Accession"] = clean_frame["Accession"].str.replace(r'\#.+\#[^\:]+\:', '')
        clean_frame["Accession"] = clean_frame["Accession"].str.replace(r'\:\#.+\#[^\:]+$', '')
        clean_frame["Accession"] = clean_frame["Accession"].str.replace(r'\#.+\#[^\:]+$', 'NaN')
        clean_frame = clean_frame.replace("NaN", np.nan)
        #clean_frame["Peptide"] = clean_frame["Peptide"].str.replace(r'\(.*\.\d*\)', '')
        clean_frame["Peptide"] = clean_frame["Peptide"].str.replace(r'\(.{4,7}\)', '')
        clean_frame = clean_frame.dropna(subset = ["Accession"]).reset_index(drop = True)

        # Calculate the initial number of rows
        initial_count = len(clean_frame)
        clean_frame = clean_frame[~clean_frame['Accession'].str.contains('CONTAM|DECOY', na=False)].reset_index(drop=True)
        # Calculate the number of rows dropped
        rows_dropped = initial_count - len(clean_frame)
        # Print the number of rows dropped
        print(f"Number of rows dropped where Accession had 'CONTAN' or 'DECOY' marked,: {rows_dropped}")

        self.peptide_frame = clean_frame

    def get_locs(self):
        """Finds the location of all the peptides stored in self.protein frame within their corresponding protein sequence

        Return:
        """
        locs = []
        for row in range(len(self.protein_frame)):
            #find the position of the peptide within the origin protein  
            loc = str(self.protein_frame["Sequence"][row]).find(self.protein_frame["Peptide"][row])
            if loc == -1:
                locs.append("No Position Found")
            else:
                locs.append(loc)
        return locs
        
    def get_rel_locs(self):
        """Finds the relative location (in decimal) of all the peptides stored in self.protein frame within their corresponding protein sequence
        """
        rel_locs = []
        for row in range(len(self.protein_frame)):
            if self.protein_frame["Location"][row] == "No Position Found":
                rel_locs.append("NA")
            else:
                rel_locs.append(self.protein_frame["Location"][row]/len(self.protein_frame["Sequence"][row]))
        return rel_locs

    def search_proteome(self, peptide_file, proteome_file = "uniprot-proteome_UP000005640.fasta"):
        """Searches a proteome (proteome_file) and retrieves all protein sequences referenced by protein IDs corresponding to peptides stored in a PEAKS output csv (peptide_file)

        Keyword Arguments:
        peptide_file -- a csv file containing the PEAKS output of immunopeptidome sequencing
        proteome_file -- filename of the .fasta file containing the proteome
        """
        self.filename = peptide_file
        #Load proteome into a dataframe
        proteome_frame = self.load_proteome(proteome_file)
        #Load peptides into a dataframe
        self.peptide_frame = pd.read_csv(peptide_file)
        self.clean()
        ids = []
        use_id = []

        #Retrieve most confident accession for each peptide
        for all_accession in self.peptide_frame["Accession"]:
            all_accession = all_accession.split(":")
            for accession in range(len(all_accession)):
                all_accession[accession] = all_accession[accession].split("|")

                for entry in all_accession[accession]:

                    if len(entry) <4:
                        continue
                    if entry.find('_') != -1:
                        continue
                    # elif entry in accessionids:
                    else:
                        all_accession[accession] = entry.split('-')[0]
                        break

                # all_accession[accession] = all_accession[accession][1].split('-')[0]
            ids.append(all_accession)
            use_id.append(all_accession[0])
        #add the ID of the origin protein
        self.peptide_frame["ID"] = use_id
        self.peptide_frame["ID List"] = ids
        self.peptide_frame["Location"] = "No protein found"
        #Add peptide count column to proteome to keep track of peptides/protein
        proteome_frame["Peptides"] = 0
        
        # Keep only rows with valid string/int IDs
        valid_peptide_frame = self.peptide_frame[self.peptide_frame["ID"].apply(lambda x: isinstance(x, (str, int)))]

        self.protein_frame = valid_peptide_frame[["ID", "Peptide"]].join(
            proteome_frame[["ID", "Sequence"]].set_index("ID"),
            on="ID",
            how="inner",
            lsuffix="_pep",
            rsuffix="_prot"
        )

        #retrieve locations for all peptides
        self.protein_frame.reset_index(drop = True, inplace = True)
        self.protein_frame["Location"] = self.get_locs()
        self.protein_frame["Relative Location"] = self.get_rel_locs()

    def protein_dist(self, protein):
        """Generates a list of length 100 with each item representing the number of detected peptides at that relative location within the source protein

        Keyword Arguments:
        protein -- source protein to search file for peptide fragments

        Return: A list containing the number of peptides at each relative position within the protein sequence
        """
        #find all peptides assigned to source protein
        peptides = self.protein_frame[self.protein_frame["ID"] == protein]
        peptides.reset_index(drop = True, inplace = True)
        locations = [0]*100

        #calculate relative length of peptides
        for pep in range(len(peptides)):

            # Ammend by Prithvi: Skipping adding the location of the peptide if peptide not found in the protein
            if peptides["Location"][pep] == "No Position Found": 
                continue

            start = int(peptides["Relative Location"][pep]//0.01)
            end = int(((int(peptides["Location"][pep]) + len(peptides["Peptide"][pep]) - 1)/len(peptides["Sequence"][pep]))//0.01)
            #increment count for each relative position
            for position in range(start, end + 1):
                locations[position] += 1
        return locations

    def peptide_dist(self, proteins, taskId):
        """Displays a heatmap of peptide locations on selected peptides

        Keyword Arguments:
        proteins -- accession code or list of accession codes of proteins to display peptide distribution maps of  
        """
        if type(proteins) is str:
            proteins = [proteins]
        protein_list = []
        for accession in proteins:
            protein_list.append(self.protein_dist(accession))
        plt.clf()
        g = sns.heatmap(protein_list, cmap = "Reds", yticklabels = proteins, cbar_kws={'label': 'Number of unique peptides'})
        g.set_yticklabels(g.get_yticklabels(), rotation=0, horizontalalignment='right')

        # fig,ax=plt.subplots(figsize=(6,6))
        # canvas=FigureCanvas(fig)
        # img = io.BytesIO()
        # fig.savefig(img)
        # img.seek(0)
        plt.xlabel("Normalised protein length") 
        plt.savefig(os.path.join(project_root,'app','static','images',taskId,'pepscanner.png'), bbox_inches='tight')
        print("Pepscanner Heatplot saved at :"+ os.path.join(project_root,'app','static','images',taskId,'pepscanner.png'))

        # return img

if __name__ == "__main__":
    #initialise object
    scanner = PepScan()
    #search proteome file for source proteins of peptide file
    # scanner.search_proteome(peptide_file = "VMM1_1st_nil_DT9_peptide_C0303_NetMHCpan_binders.csv", proteome_file = "uniprot-proteome_UP000005640.fasta")
    scanner.search_proteome(peptide_file = "elutiondata.csv", proteome_file = "uniprot-proteome_UP000005640.fasta")

    scanner.peptide_dist(["Q96C19", "P04818", "F8W6I7", "A0A0U1RQF0"])
