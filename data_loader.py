
import pandas as pd

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



#load data

def clean_msgfplus(file_name):
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
    msgfplus_files[".2ng"] = "data/msgfplus/Ex_Auto_DrM3_30umT4_2ngQC_60m_half.tsv.gz"
    msgfplus_files["2ng"] = "data/msgfplus/Ex_Auto_DrM3_30umT4_02ngQC_60m_half.tsv.gz"


    file_path = msgfplus_files.get(file_name)

    df = pd.read_csv(file_path, sep='\t')#, sep='\t', header=0, index_col=0)

    df['decoy'] = df.apply (lambda row: make_decoy_col_msgf(row), axis=1)
    df = df.rename({'ScanNum': 'scan', 'Peptide': 'peptide'}, axis=1)
    df = df.filter(['decoy', 'scan', 'peptide', 'QValue'])

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

        df = df.rename({"PEP.IsDecoy": "decoy", "PSM.MS2ScanNumber": "scan", "PEP.StrippedSequence": "peptide"}, axis=1)
        df = df.filter(['decoy', 'scan', 'peptide', 'PEP.QValue'])

        return df

def clean_msfragger(file_name):
    msfragger_files = {}

    #2 ng
    msfragger_files["2ng"] = "data/msfragger/Ex_Auto_DrM3_30umT4_2ngQC_60m_halfpsm.tsv.gz"
    msfragger_files[".2ng"] = "data/msfragger/Ex_Auto_DrM3_30umT4_02ngQC_60m_halfpsm.tsv.gz"


    file_path = msfragger_files.get(file_name)
    df = pd.read_csv(file_path, sep='\t')#, sep='\t', header=0, index_col=0)

    df = df.rename({"Peptide":"peptide", "Spectrum":"scan"}, axis=1)
    df["decoy"] = df.apply(lambda row: make_decoy_col_msfragger(row), axis=1)
    df = df.filter(['decoy', 'scan', 'peptide', 'PeptideProphet Probability'])

    return df


def parse_meta(file_name, cutoff):
    meta_files = {}
    #Single cell
    meta_files["singleCell_1"] = "data/MetaMorpheus/RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep1_Peptides.psmtsv.gz"
    meta_files["singleCell_2"] = "data/MetaMorpheus/RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep2_Peptides.psmtsv.gz"
    meta_files["singleCell_3"] = "data/MetaMorpheus/RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep3_Peptides.psmtsv.gz"
    meta_files["singleCell_4"] = "data/MetaMorpheus/RC1051_DDA_SingleCell_HeLa_1-14-2021_Rep1_Peptides.psmtsv.gz"
    #50 ng
    meta_files["50ng_1"] = "data/MetaMorpheus/RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep1_Peptides.psmtsv.gz"
    meta_files["50ng_2"] = "data/MetaMorpheus/RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep2_Peptides.psmtsv.gz"
    meta_files["50ng_3"] = "data/MetaMorpheus/RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep3_Peptides.psmtsv.gz"
    meta_files["50ng_4"] = "data/MetaMorpheus/RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep4_Peptides.psmtsv.gz"
    meta_files["50ng_5"] = "data/MetaMorpheus/RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep5_Peptides.psmtsv.gz"

    #2ng run
    meta_files[".2ng"] = ""

    data = pd.read_csv(meta_files.get(file_name), sep="\t")

    data.insert(len(data.columns), 'Decoy_Id', 1)

    decoy_count = 0;
    for index, row in data['Cumulative Decoy'].iteritems():
        decoy_num = int(data['Cumulative Decoy'][index])
        if(decoy_num > decoy_count):
            decoy_count +=1
            data.iloc[index,len(data.columns)-1] = 'NaN'


    decoys_dropped_df = data[data.Decoy_Id != "NaN"]

    no_duplicates_df = decoys_dropped_df.drop_duplicates(subset=['Scan Number'])
    only_unique_df = no_duplicates_df.drop_duplicates(subset=['Base Sequence'])
    cut_off_df = only_unique_df[only_unique_df.PEP_QValue <= 0.01]

    return cut_off_df






def print_stats(df):
    #Find how many unique peptides there are
    unique_pep = df['new_peptide'].nunique()
    number_uniq_pep = df['new_peptide'].nunique()

    #fine how many unique spectra there are
    number_spectra = df["SpecID"].nunique()

    print("Number of unique peptides:", number_uniq_pep,"   " "Number of unique spectra:", number_spectra)


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
