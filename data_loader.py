
import pandas as pd
import numpy as np

def make_decoy_col_msgf(row):
    if row["Protein"].startswith("XXX_"):
        return True
    else:
        return False

def make_decoy_col_msfragger(row):
    if row["Protein"].startswith("rev"):
        return True
    else:
        return False

def make_decoy_col_maxquant(row):
    if row["Proteins"].startswith("REV"):
        return True
    else:
        return False

#how we're going to format the peptide modifictations:
# +15.995
# +57.021
#works for msfragger
def format_oxidation(row, column, to_replace):
    peptide = row[column]
#     print(to_replace)
    replace_with = "+15.995"
    if pd.isna(peptide):
        new_pep = peptide
    else:
        if to_replace in peptide:
            new_pep = peptide.replace(to_replace, replace_with)
        else:
            new_pep = peptide
    return new_pep

def format_carbamidomethyl(row, column, to_replace):
    peptide = row[column]
#     print(to_replace)
    replace_with = ""
    if pd.isna(peptide):
        new_pep = peptide
    else:
        if to_replace in peptide:
            new_pep = peptide.replace(to_replace, replace_with)
        else:
            new_pep = peptide
    return new_pep

#load data
def clean_msgfplus(file_name, prob_method="PEP"):
    msgfplus_files = {}
    #Single cell
    msgfplus_files["singleCell_1"] = "data/msgfplus/RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep1.gz"
    msgfplus_files["singleCell_2"] = "data/msgfplus/RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep2.gz"
    msgfplus_files["singleCell_3"] = "data/msgfplus/RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep3.gz"
    msgfplus_files["singleCell_4"] = "data/msgfplus/RC1051_DDA_SingleCell_HeLa_1-14-2021_Rep1.gz"
    #50 ng
    msgfplus_files["50ng_1"] = "data/msgfplus/RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep1.gz"
    msgfplus_files["50ng_2"] = "data/msgfplus/RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep2.gz"
    msgfplus_files["50ng_3"] = "data/msgfplus/RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep3.gz"
    msgfplus_files["50ng_4"] = "data/msgfplus/RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep4.gz"
    msgfplus_files["50ng_5"] = "data/msgfplus/RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep5.gz"
    #2 ng
    msgfplus_files["2ng"] = "data/msgfplus/Ex_Auto_DrM3_30umT4_2ngQC_60m_half.tsv.gz"
    msgfplus_files[".2ng"] = "data/msgfplus/Ex_Auto_DrM3_30umT4_02ngQC_60m_half.tsv.gz"


    file_path = msgfplus_files.get(file_name)

    df = pd.read_csv(file_path, sep='\t')#, sep='\t', header=0, index_col=0)

    #make a new col that includes modifide peptides
    #it's already formatted correly for oxidation
    df["new_pep"] = df.apply(lambda row: format_carbamidomethyl(row, "Peptide", "+57.021"), axis=1)

    df['decoy'] = df.apply (lambda row: make_decoy_col_msgf(row), axis=1)
    if prob_method=="PEP":
        df = df.rename({'ScanNum': 'scan', 'new_pep': 'peptide', 'PepQValue': 'probability'}, axis=1)
    else:
        df = df.rename({'ScanNum': 'scan', 'new_pep': 'peptide', 'QValue': 'probability'}, axis=1)

    # df = df.filter(['decoy', 'scan', 'peptide', 'probability'])

    return df

def clean_spectromine(file_name):
        combined_df = pd.read_csv("data/spectromine/20210129_140856_SingleCell_PSM Report_20210201_171706.csv.gz")

        #we're going to have to spearate out the files based on file name.
        spectro_files = {}
        spectro_files["singleCell_1"] = "RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep1.raw"
        spectro_files["singleCell_2"] = "RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep2.raw"
        spectro_files["singleCell_3"] = "RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep3.raw"
        spectro_files["singleCell_4"] = "RC1051_DDA_SingleCell_HeLa_1-14-2021_Rep1.raw"
        spectro_files["50ng_1"] = "RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep1.raw"
        spectro_files["50ng_2"] = "RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep2.raw"
        spectro_files["50ng_3"] = "RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep3.raw"
        spectro_files["50ng_4"] = "RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep4.raw"
        spectro_files["50ng_5"] = "RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep5.raw"

        spectro_files["2ng"] = "Ex_Auto_DrM3_30umT4_2ngQC_60m_half.raw"
        spectro_files[".2ng"] = "Ex_Auto_DrM3_30umT4_02ngQC_60m_half.raw"

        #then use the file name to select
        file_path = spectro_files.get(file_name)
        df = combined_df[combined_df["R.FileName"]==file_path]

        df = df.rename({"PEP.IsDecoy": "decoy", "PSM.MS2ScanNumber": "scan", "PEP.StrippedSequence": "peptide", "PEP.QValue": "probability"}, axis=1)
        # df = df.filter(['decoy', 'scan', 'peptide', 'probability'])

        return df

def clean_msfragger(file_name):
    msfragger_files = {}

    #2 ng
    msfragger_files["2ng"] = "data/msfragger/Ex_Auto_DrM3_30umT4_2ngQC_60m_halfpsm.tsv.gz"
    msfragger_files[".2ng"] = "data/msfragger/Ex_Auto_DrM3_30umT4_02ngQC_60m_halfpsm.tsv.gz"

    file_path = msfragger_files.get(file_name)
    df = pd.read_csv(file_path, sep='\t')

    #make a new col that includes modifide peptides
    df['temp_peptide'] = df.apply(lambda row: format_oxidation(row, "Modified Peptide", "[147]"), axis=1)
    df['temp2'] = np.where(pd.isna(df['temp_peptide']), df['Peptide'], df['temp_peptide'])


    df = df.rename({"temp2":"peptide", "Spectrum":"scan", 'PeptideProphet Probability': 'probability'}, axis=1)

    df["decoy"] = df.apply(lambda row: make_decoy_col_msfragger(row), axis=1)
    # df = df.filter(['decoy', 'scan', 'peptide', 'probability'])

    return df

def clean_metamorph(file_name, prob_method="PEP"):
    mm_files = {}

    mm_files[".2ng"] = 'data/MetaMorpheus/Ex_Auto_DrM3_30umT4_02ngQC_60m_half_PSMs.psmtsv.gz'
    mm_files["2ng"] = 'data/MetaMorpheus/Ex_Auto_DrM3_30umT4_2ngQC_60m_half_PSMs.psmtsv.gz'

    input_file = mm_files.get(file_name)
    data = pd.read_csv(input_file, sep = '\t')

    #make a new col that includes modifide peptides
    data['temp_peptide'] = data.apply(lambda row: format_oxidation(row, "Full Sequence", "[Common Variable:Oxidation on M]"), axis=1)
    data["temp2"] = data.apply(lambda row: format_carbamidomethyl(row, "temp_peptide", "[Common Fixed:Carbamidomethyl on C]"), axis=1)

    data = data.replace({"Decoy": {'Y': True, 'N': False}})
    #uniform naming
    if prob_method=="PEP":
        data_new = data.rename({"Decoy": "decoy", "Scan Number": "scan", "temp2": "peptide", 'PEP_QValue': 'probability'}, axis=1)
    else:
        data_new = data.rename({"Decoy": "decoy", "Scan Number": "scan", "temp2": "peptide", 'QValue': 'probability'}, axis=1)

    # data_new = data_new.filter((['decoy', 'scan', 'peptide','probability' ]))

    return data_new

def clean_maxquant(file_name, prob_method="PEP"):
    combined_df = pd.read_csv("data/maxquant/msmsScans.txt.gz", sep="\t")

    #get the 2 ng file
    maxq_files = {}
    maxq_files["2ng"] = "Ex_Auto_DrM3_30umT4_2ngQC_60m_half"
    maxq_files[".2ng"] = "Ex_Auto_DrM3_30umT4_02ngQC_60m_half"
    file_path = maxq_files.get(file_name)
    df = combined_df[combined_df["Raw file"]==file_path]

    #make a new col that includes modifide
    df['temp_peptide'] = df.apply(lambda row: format_oxidation(row, "Modified sequence", "(Oxidation (M))"), axis=1)
    df["temp_peptide"] = df["temp_peptide"].str[1:-1]

    df["decoy"] = df.apply(lambda row: make_decoy_col_maxquant(row), axis=1)

    if prob_method=="PEP":
        df = df.rename({"Scan number": "scan", "temp_peptide": "peptide", 'PEP':'probability'}, axis=1)
    else:
        df = df.rename({"Scan number": "scan", "temp_peptide": "peptide", 'Score':'probability'}, axis=1)

    # df = df.filter(['decoy', 'scan', 'peptide', 'probability'])

    return df


#set cutoff
#can choose column
#Use the cuttoff of 1 pweacent ) (.01)
#how many unique peptides

# only report scan numbers
# remove decoy

# number of unique peptides

# unique scan (spectra)
# unique peptide

#report data
