
import pandas as pd

#load data
def parse_msgfplus(file_name, cutoff):

    # file_dict[1] = "peptide_ID/RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep1.gz"
    # file = "peptide_ID/RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep1.gz"

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
    msgfplus_files["2ng_1"] = "data/msgfplus/RC1051_DDA_QC_HeLa_2ng_12-28-2020.gz"
    msgfplus_files["2ng_2_after_outage"] = "data/msgfplus/RC1051_QC_After_outage_2ng_QC_HeLa_12-29-2020.gz"


    file_path = msgfplus_files.get(file_name)

    df = pd.read_csv(file_path, sep='\t')#, sep='\t', header=0, index_col=0)

    #drop decoy
    df = df[~df["Protein"].str.startswith("XXX_")]

    #drop duplicate scans
    df = df.drop_duplicates(subset=["ScanNum"], keep="first") #keep highest coring

    #drop not unique peptides (keep highest scoring one
    #Add col to help find how many unique peptides there are
    df["new_peptide"] = df["Peptide"].str[2:-2]

    #set cutoff
    df = df[df.QValue <= cutoff]

    return df
def parse_spectromine(file_name, cutoff):
        combined_df = pd.read_csv("data/spectromine/20210129_140856_SingleCell_PSM Report_20210201_171706.csv")

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

        #then use the file name to select
        file_path = spectro_files.get(file_name)
        df = combined_df[combined_df["R.FileName"]==file_path]

        #drop decoy
        df = df[~df["PEP.IsDecoy"] == True]

        #drop duplicate scans
        df = df.drop_duplicates(subset=["PSM.MS2ScanNumber"], keep="first")

        #set_cutoff
        df = df[df["PEP.QValue"] <= cutoff]

        return df
def parse_msfragger(file_name, cutoff):
        msfragger_files = {}
        #Single cell
        msfragger_files["singleCell_1"] = "data/msfragger/RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep1psm.tsv.gz"
        msfragger_files["singleCell_2"] = "data/msfragger/RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep2psm.tsv.gz"
        msfragger_files["singleCell_3"] = "data/msfragger/RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep3psm.tsv.gz"
        msfragger_files["singleCell_4"] = "data/msfragger/RC1051_DDA_SingleCell_HeLa_1-14-2021_Rep1psm.tsv.gz"
        #50 ng
        msfragger_files["50ng_1"] = "data/msfragger/RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep1psm.tsv.gz"
        msfragger_files["50ng_2"] = "data/msfragger/RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep2psm.tsv.gz"
        msfragger_files["50ng_3"] = "data/msfragger/RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep3psm.tsv.gz"
        msfragger_files["50ng_4"] = "data/msfragger/RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep4psm.tsv.gz"
        msfragger_files["50ng_5"] = "data/msfragger/RC1051_Library_DDA_QC_HeLa_50ng_12-28-2020_Rep5psm.tsv.gz"
        #2 ng
        msfragger_files["2ng"] = "data/msfragger/Ex_Auto_DrM3_30umT4_2ngQC_60m_halfpsm.tsv.gz"
        msfragger_files[".2"] = "data/msfragger/Ex_Auto_DrM3_30umT4_02ngQC_60m_halfpsm.tsv.gz"


        file_path = msfragger_files.get(file_name)

        df = pd.read_csv(file_path, sep='\t')#, sep='\t', header=0, index_col=0)

        return df


#These functions will probably be deleted
def parse_msgfplus_remove_decoy(file_name):

    # file_dict[1] = "peptide_ID/RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep1.gz"
    # file = "peptide_ID/RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep1.gz"

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
    msgfplus_files["2ng_1"] = "data/msgfplus/RC1051_DDA_QC_HeLa_2ng_12-28-2020.gz"
    msgfplus_files["2ng_2_after_outage"] = "data/msgfplus/RC1051_QC_After_outage_2ng_QC_HeLa_12-29-2020.gz"


    file_path = msgfplus_files.get(file_name)

    df = pd.read_csv(file_path, sep='\t')#, sep='\t', header=0, index_col=0)

    #drop decoy
    df = df[~df["Protein"].str.startswith("XXX_")]

    #Add col to help find how many unique peptides there are
    df["new_peptide"] = df["Peptide"].str[2:-2]


    return df
def parse_spectromine_remove_decoy(file_name):
        combined_df = pd.read_csv("data/spectromine/20210129_140856_SingleCell_PSM Report_20210201_171706.csv")

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



        #then use the file name to select
        file_path = spectro_files.get(file_name)
        df = combined_df[combined_df["R.FileName"]==file_path]

        #drop decoy
        df = df[~df["PEP.IsDecoy"] == True]

        return df

def parse_msgfplus_no_cutoff(file_name):

    # file_dict[1] = "peptide_ID/RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep1.gz"
    # file = "peptide_ID/RC1051_DDA_SingleCell_HeLa_1-1-2021_Rep1.gz"

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
    msgfplus_files["2ng_1"] = "data/msgfplus/RC1051_DDA_QC_HeLa_2ng_12-28-2020.gz"
    msgfplus_files["2ng_2_after_outage"] = "data/msgfplus/RC1051_QC_After_outage_2ng_QC_HeLa_12-29-2020.gz"


    file_path = msgfplus_files.get(file_name)

    df = pd.read_csv(file_path, sep='\t')#, sep='\t', header=0, index_col=0)

    #drop decoy
    df = df[~df["Protein"].str.startswith("XXX_")]

    #drop duplicate scans
    df = df.drop_duplicates(subset=["ScanNum"], keep="first")

    #Add col to help find how many unique peptides there are
    df["new_peptide"] = df["Peptide"].str[2:-2]

    return df
def parse_spectromine_no_cutoff(file_name):
        combined_df = pd.read_csv("data/spectromine/20210129_140856_SingleCell_PSM Report_20210201_171706.csv")

        #we're going to have to spearate out the files based on file name.
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

        #then use the file name to select
        file_path = spectro_files.get(file_name)
        df = combined_df[combined_df["R.FileName"]==file_path]

        #drop decoy
        df = df[~df["PEP.IsDecoy"] == True]

        #drop duplicate scans
        df = df.drop_duplicates(subset=["PSM.MS2ScanNumber"], keep="first")

        return df


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
