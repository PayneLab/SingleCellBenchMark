
import pandas as pd
import numpy as np
import os

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
def clean_msgfplus(file_name):
    msgfplus_files = {}
    #2ng
    msgfplus_files["2ng_rep1"] = "data/msgfplus/Ex_Auto_J3_30umTB_2ngQC_60m_1.tsv.gz"
    msgfplus_files["2ng_rep2"] = "data/msgfplus/Ex_Auto_J3_30umTB_2ngQC_60m_2.tsv.gz"
    msgfplus_files["2ng_rep3"] = "data/msgfplus/Ex_Auto_K13_30umTA_2ngQC_60m_1.tsv.gz"
    msgfplus_files["2ng_rep4"] = "data/msgfplus/Ex_Auto_K13_30umTA_2ngQC_60m_2.tsv.gz"
    msgfplus_files["2ng_rep5"] = "data/msgfplus/Ex_Auto_W17_30umTB_2ngQC_60m_1.tsv.gz"
    msgfplus_files["2ng_rep6"] = "data/msgfplus/Ex_Auto_W17_30umTB_2ngQC_60m_2.tsv.gz"

    #0.2 ng
    msgfplus_files["0.2ng_rep1"] = "data/msgfplus/Ex_Auto_J3_30umTB_02ngQC_60m_1.tsv.gz"
    msgfplus_files["0.2ng_rep2"] = "data/msgfplus/Ex_Auto_J3_30umTB_02ngQC_60m_2.tsv.gz"
    msgfplus_files["0.2ng_rep3"] = "data/msgfplus/Ex_Auto_K13_30umTA_02ngQC_60m_1.tsv.gz"
    msgfplus_files["0.2ng_rep4"] = "data/msgfplus/Ex_Auto_K13_30umTA_02ngQC_60m_2.tsv.gz"
    msgfplus_files["0.2ng_rep5"] = "data/msgfplus/Ex_Auto_W17_30umTA_02ngQC_60m_3.tsv.gz"
    msgfplus_files["0.2ng_rep6"] = "data/msgfplus/Ex_Auto_W17_30umTA_02ngQC_60m_4.tsv.gz"

    path_to_data_loader = os.path.abspath(os.path.dirname(__file__)) # This gets the absolute path to the location of the data_loader.py file
    complete_path_to_data = os.path.join(path_to_data_loader, msgfplus_files.get(file_name)) # We then append the relative path to the data files

    df = pd.read_csv(complete_path_to_data, sep = '\t')

    #make a new col that includes modifide peptides
    #it's already formatted correly for oxidation
    df["new_pep"] = df.apply(lambda row: format_carbamidomethyl(row, "Peptide", "+57.021"), axis=1)

    df['decoy'] = df.apply (lambda row: make_decoy_col_msgf(row), axis=1)

    df = df.rename({'ScanNum': 'scan', 'new_pep': 'peptide'}, axis=1)


    return df

def clean_msfragger(file_name):
    #it came as a combined file, so we need to parse out individual mm_files
    path_to_data_loader = os.path.abspath(os.path.dirname(__file__)) # This gets the absolute path to the location of the data_loader.py file
    complete_path_to_data = os.path.join(path_to_data_loader, 'data/msfragger/psm.tsv.gz') # We then append the relative path to the data files

    combined_df = pd.read_csv(complete_path_to_data, sep = '\t') #combined file

    msfragger_files = {}
    msfragger_files["2ng_rep1"] = "Ex_Auto_J3_30umTB_2ngQC_60m_1"
    msfragger_files["2ng_rep2"] = "Ex_Auto_J3_30umTB_2ngQC_60m_2"
    msfragger_files["2ng_rep3"] = "Ex_Auto_K13_30umTA_2ngQC_60m_1"
    msfragger_files["2ng_rep4"] = "Ex_Auto_K13_30umTA_2ngQC_60m_2"
    msfragger_files["2ng_rep5"] = "Ex_Auto_W17_30umTB_2ngQC_60m_1"
    msfragger_files["2ng_rep6"] = "Ex_Auto_W17_30umTB_2ngQC_60m_2"

    #0.2 ng
    msfragger_files["0.2ng_rep1"] = "Ex_Auto_J3_30umTB_02ngQC_60m_1"
    msfragger_files["0.2ng_rep2"] = "Ex_Auto_J3_30umTB_02ngQC_60m_2"
    msfragger_files["0.2ng_rep3"] = "Ex_Auto_K13_30umTA_02ngQC_60m_1"
    msfragger_files["0.2ng_rep4"] = "Ex_Auto_K13_30umTA_02ngQC_60m_2"
    msfragger_files["0.2ng_rep5"] = "Ex_Auto_W17_30umTA_02ngQC_60m_3"
    msfragger_files["0.2ng_rep6"] = "Ex_Auto_W17_30umTA_02ngQC_60m_4"

    file_path = msfragger_files.get(file_name)
    df = combined_df[combined_df["Spectrum File"].str.contains(file_path)]

    #make a new col that includes modifide peptides
    # df['temp_peptide'] = df.apply(lambda row: format_oxidation(row, "Modified Peptide", "[147]"), axis=1)
    # df['temp2'] = np.where(pd.isna(df['temp_peptide']), df['Peptide'], df['temp_peptide'])

    df = df.assign(temp_peptide=df.apply(lambda row: format_oxidation(row, "Modified Peptide", "[147]"), axis=1))
    df = df.assign(temp2=np.where(pd.isna(df['temp_peptide']), df['Peptide'], df['temp_peptide']))

    df = df.rename({"temp2":"peptide", "Spectrum":"scan"}, axis=1)

    df["decoy"] = df.apply(lambda row: make_decoy_col_msfragger(row), axis=1)
    # df = df.filter(['decoy', 'scan', 'peptide', 'probability'])

    return df

def clean_metamorph(file_name):
    mm_files = {}

    #2ng
    mm_files["2ng_rep1"] = "data/MetaMorpheus/Ex_Auto_J3_30umTB_2ngQC_60m_1-calib_PSMs.psmtsv.gz"
    mm_files["2ng_rep2"] = "data/MetaMorpheus/Ex_Auto_J3_30umTB_2ngQC_60m_2-calib_PSMs.psmtsv.gz"
    mm_files["2ng_rep3"] = "data/MetaMorpheus/Ex_Auto_K13_30umTA_2ngQC_60m_1-calib_PSMs.psmtsv.gz"
    mm_files["2ng_rep4"] = "data/MetaMorpheus/Ex_Auto_K13_30umTA_2ngQC_60m_2-calib_PSMs.psmtsv.gz"
    mm_files["2ng_rep5"] = "data/MetaMorpheus/Ex_Auto_W17_30umTB_2ngQC_60m_1-calib_PSMs.psmtsv.gz"
    mm_files["2ng_rep6"] = "data/MetaMorpheus/Ex_Auto_W17_30umTB_2ngQC_60m_2-calib_PSMs.psmtsv.gz"

    #0.2 ng
    mm_files["0.2ng_rep1"] = "data/MetaMorpheus/Ex_Auto_J3_30umTB_02ngQC_60m_1-calib_PSMs.psmtsv.gz"
    mm_files["0.2ng_rep2"] = "data/MetaMorpheus/Ex_Auto_J3_30umTB_02ngQC_60m_2-calib_PSMs.psmtsv.gz"
    mm_files["0.2ng_rep3"] = "data/MetaMorpheus/Ex_Auto_K13_30umTA_02ngQC_60m_1-calib_PSMs.psmtsv.gz"
    mm_files["0.2ng_rep4"] = "data/MetaMorpheus/Ex_Auto_K13_30umTA_02ngQC_60m_2-calib_PSMs.psmtsv.gz"
    mm_files["0.2ng_rep5"] = "data/MetaMorpheus/Ex_Auto_W17_30umTA_02ngQC_60m_3-calib_PSMs.psmtsv.gz"
    mm_files["0.2ng_rep6"] = "data/MetaMorpheus/Ex_Auto_W17_30umTA_02ngQC_60m_4-calib_PSMs.psmtsv.gz"


    path_to_data_loader = os.path.abspath(os.path.dirname(__file__)) # This gets the absolute path to the location of the data_loader.py file
    complete_path_to_data = os.path.join(path_to_data_loader, mm_files.get(file_name)) # We then append the relative path to the data files

    data = pd.read_csv(complete_path_to_data, sep = '\t')

    #make a new col that includes modifide peptides
    data['temp_peptide'] = data.apply(lambda row: format_oxidation(row, "Full Sequence", "[Common Variable:Oxidation on M]"), axis=1)
    data["temp2"] = data.apply(lambda row: format_carbamidomethyl(row, "temp_peptide", "[Common Fixed:Carbamidomethyl on C]"), axis=1)

    data = data.replace({"Decoy": {'Y': True, 'N': False}})
    #uniform naming
    data_new = data.rename({"Decoy": "decoy", "Scan Number": "scan", "temp2": "peptide"}, axis=1)

    # data_new = data_new.filter((['decoy', 'scan', 'peptide','probability' ]))

    return data_new

def clean_maxquant(file_name):
    path_to_data_loader = os.path.abspath(os.path.dirname(__file__)) # This gets the absolute path to the location of the data_loader.py file

    complete_path_to_2ngdata = os.path.join(path_to_data_loader, "data/maxquant/msmsScans2ng.txt.gz") # We then append the relative path to the data files
    combined_df_2ng = pd.read_csv(complete_path_to_2ngdata, sep = '\t')

    complete_path_to_02ngdata = os.path.join(path_to_data_loader, "data/maxquant/msmsScans02ng.txt.gz") # We then append the relative path to the data files
    combined_df_02ng = pd.read_csv(complete_path_to_02ngdata, sep = '\t')

    #get the 2 ng file
    maxq_files = {}
    maxq_files["2ng_rep1"] = "Ex_Auto_J3_30umTB_2ngQC_60m_1"
    maxq_files["2ng_rep2"] = 'Ex_Auto_J3_30umTB_2ngQC_60m_2'
    maxq_files["2ng_rep3"] = 'Ex_Auto_K13_30umTA_2ngQC_60m_1'
    maxq_files["2ng_rep4"] = 'Ex_Auto_K13_30umTA_2ngQC_60m_2'
    maxq_files["2ng_rep5"] = 'Ex_Auto_W17_30umTB_2ngQC_60m_1'
    maxq_files["2ng_rep6"] = 'Ex_Auto_W17_30umTB_2ngQC_60m_2'

    maxq_files["0.2ng_rep1"] = "Ex_Auto_J3_30umTB_02ngQC_60m_1"
    maxq_files["0.2ng_rep2"] = 'Ex_Auto_J3_30umTB_02ngQC_60m_2'
    maxq_files["0.2ng_rep3"] = 'Ex_Auto_K13_30umTA_02ngQC_60m_1'
    maxq_files["0.2ng_rep4"] = 'Ex_Auto_K13_30umTA_02ngQC_60m_2'
    maxq_files["0.2ng_rep5"] = 'Ex_Auto_W17_30umTA_02ngQC_60m_3'
    maxq_files["0.2ng_rep6"] = 'Ex_Auto_W17_30umTA_02ngQC_60m_4'

    if "0.2" in file_name:
        file_path = maxq_files.get(file_name)
        df = combined_df_02ng[combined_df_02ng["Raw file"]==file_path]
    else:
        file_path = maxq_files.get(file_name)
        df = combined_df_2ng[combined_df_2ng["Raw file"]==file_path]

    # file_path = maxq_files.get(file_name)
    # df = combined_df[combined_df["Raw file"]==file_path]

    #make a new col that includes modifide
    df['temp_peptide'] = df.apply(lambda row: format_oxidation(row, "Modified sequence", "(Oxidation (M))"), axis=1)
    df["temp_peptide"] = df["temp_peptide"].str[1:-1]

    df["decoy"] = df.apply(lambda row: make_decoy_col_maxquant(row), axis=1)

    df = df.rename({"Scan number": "scan", "temp_peptide": "peptide"}, axis=1)
    # else:
    #     df = df.rename({"Scan number": "scan", "temp_peptide": "peptide", 'Score':'probability'}, axis=1)

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

        df = df.rename({"PEP.IsDecoy": "decoy", "PSM.MS2ScanNumber": "scan", "PEP.StrippedSequence": "peptide", "PEP.QValue": "PEP_QValue"}, axis=1)
        # df = df.filter(['decoy', 'scan', 'peptide', 'probability'])

        return df
    
def get_pin_file(file_name):
    mm_pin_files = {}

    #2ng
    mm_pin_files["2ng_rep1"] = "data/MetaMorpheus/pin_files/Ex_Auto_J3_30umTB_2ngQC_60m_1-calib_PSMsFormattedForPercolator.tab.gz"
    mm_pin_files["2ng_rep2"] = "data/MetaMorpheus/pin_files/Ex_Auto_J3_30umTB_2ngQC_60m_2-calib_PSMsFormattedForPercolator.tab.gz"
    mm_pin_files["2ng_rep3"] = "data/MetaMorpheus/pin_files/Ex_Auto_K13_30umTA_2ngQC_60m_1-calib_PSMsFormattedForPercolator.tab.gz"
    mm_pin_files["2ng_rep4"] = "data/MetaMorpheus/pin_files/Ex_Auto_K13_30umTA_2ngQC_60m_2-calib_PSMsFormattedForPercolator.tab.gz"
    mm_pin_files["2ng_rep5"] = "data/MetaMorpheus/pin_files/Ex_Auto_W17_30umTB_2ngQC_60m_1-calib_PSMsFormattedForPercolator.tab.gz"
    mm_pin_files["2ng_rep6"] = "data/MetaMorpheus/pin_files/Ex_Auto_W17_30umTB_2ngQC_60m_2-calib_PSMsFormattedForPercolator.tab.gz"

    #0.2 ng
    mm_pin_files["0.2ng_rep1"] = "data/MetaMorpheus/pin_files/Ex_Auto_J3_30umTB_02ngQC_60m_1-calib_PSMsFormattedForPercolator.tab.gz"
    mm_pin_files["0.2ng_rep2"] = "data/MetaMorpheus/pin_files/Ex_Auto_J3_30umTB_02ngQC_60m_2-calib_PSMsFormattedForPercolator.tab.gz"
    mm_pin_files["0.2ng_rep3"] = "data/MetaMorpheus/pin_files/Ex_Auto_K13_30umTA_02ngQC_60m_1-calib_PSMsFormattedForPercolator.tab.gz"
    mm_pin_files["0.2ng_rep4"] = "data/MetaMorpheus/pin_files/Ex_Auto_K13_30umTA_02ngQC_60m_2-calib_PSMsFormattedForPercolator.tab.gz"
    mm_pin_files["0.2ng_rep5"] = "data/MetaMorpheus/pin_files/Ex_Auto_W17_30umTA_02ngQC_60m_3-calib_PSMsFormattedForPercolator.tab.gz"
    mm_pin_files["0.2ng_rep6"] = "data/MetaMorpheus/pin_files/Ex_Auto_W17_30umTA_02ngQC_60m_4-calib_PSMsFormattedForPercolator.tab.gz"
    
    path_to_data_loader = os.path.abspath(os.path.dirname(__file__)) # This gets the absolute path to the location of the data_loader.py file
    complete_path_to_data = os.path.join(path_to_data_loader, mm_pin_files.get(file_name)) # We then append the relative path to the data files

    df = pd.read_csv(complete_path_to_data, sep = '\t')

    return df
