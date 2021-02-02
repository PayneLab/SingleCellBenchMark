
import pandas as pd

#load data
def load_and_clean(file_num):

    # file_dict[1] = "peptide_ID/RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep1.gz"
    # file = "peptide_ID/RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep1.gz"

    fileD = {}
    #Single cell
    fileD["singleCell_1"] = "peptide_ID/RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep1.gz"
    fileD["singleCell_2"] = "peptide_ID/RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep2.gz"
    fileD["singleCell_3"] = "peptide_ID/RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep3.gz"
    fileD["singleCell_4"] = "peptide_ID/RC1051_DDA_SingleCell_HeLa_1-14-2021_Rep1.gz"
    #50 ng
    fileD["50ng_1"] = "peptide_ID/RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep1.gz"
    fileD["50ng_2"] = "peptide_ID/RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep2.gz"
    fileD["50ng_3"] = "peptide_ID/RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep3.gz"
    fileD["50ng_4"] = "peptide_ID/RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep4.gz"
    fileD["50ng_5"] = "peptide_ID/RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep5.gz"
    #2 ng
    fileD["2ng_1"] = "peptide_ID/RC1051_DDA_QC_HeLa_2ng_12-28-2020.gz"
    fileD["2ng_2_after_outage"] = "peptide_ID/RC1051_QC_After_outage_2ng_QC_HeLa_12-29-2020.gz"

    file = fileD.get(file_num)

    df = pd.read_csv(file, sep='\t')#, sep='\t', header=0, index_col=0)

    #drop decoy
    df = df[~df["Protein"].str.startswith("XXX_")]

    #drop duplicate scans
    df = df.drop_duplicates(subset=["ScanNum"], keep="first")

    #Find how many unique peptides there are
    df["new_peptide"] = df["Peptide"].str[2:-2]


    return df

def set_cutoff(cutoff, df):
        #set cutoff
        df = df[df.QValue <= cutoff]
        return df

def print_stats(df):
    #Find how many unique peptides there are
    unique_pep = df['new_peptide'].nunique()
    number_uniq_pep = df['new_peptide'].nunique()

    #fine how many unique spectra there are
    number_spectra = df["SpecID"].nunique()

    print("Number of unique peptides: ", number_uniq_pep, "Number of unique spectra: ", number_spectra)


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
