#Script that will run the background functions for the template notebook
import os
import mokapot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append("..")
sys.path
import data_loader as dl


#This takes the data and filters it by dropping decoys and duplicates.
def filter_data(df, prob_column):
    #drop decoys
    df = df[df["decoy"] == False]
    # sort by qvalue
    df = df.sort_values(prob_column)
    # drop duplicate scans
    df = df.drop_duplicates(subset=["scan"], keep="first")  # keep highest scoring

    return df

def filter_input(input_files, probability):
    #code to parse their file into a trimmed and ready df
    formatted_files = {} 
    
    file_names = ["2ng_rep1", "2ng_rep2", "2ng_rep3", "2ng_rep4", "2ng_rep5", "2ng_rep6",
             "0.2ng_rep1", "0.2ng_rep2", "0.2ng_rep3", "0.2ng_rep4", "0.2ng_rep5", "0.2ng_rep6"]

    for file in file_names: 
        df = pd.read_csv(input_files[file], low_memory=False) 

         # make sure decoy column is a boolean column and then drop decoys
        df = df.assign(decoy=df['decoy'].astype(bool))
        df = df[df["decoy"] == False]
        # sort by qvalue
        df = df.sort_values(probability)
        # drop duplicate scans
        df = df.drop_duplicates(subset=["scan"], keep="first")  # keep highest scoring

        #Filtering out just the columns that we want
        df = df.filter(['scan', 'peptide', probability])

        formatted_files[file] = df
    return formatted_files
        
        
#Reading in all of our data from the parser and formatting it to be combined in a megascript.

# Reading in and formatting the original MetaMorpheus data
def get_orgininal_mm_data(file):
    mm_original_df = dl.clean_metamorph(file)
    mm_original_df = filter_data(mm_original_df, "QValue")

    # filter out just the ScanNr and QValue and/or PEP columns
    mm_original_df = mm_original_df.filter(items=['scan', 'peptide', 'QValue', 'PEP'])

    return mm_original_df


# Reading in and formatting the orginial MsFragger data
def set_probablility(row):
    new_prob = 1 - row["PeptideProphet Probability"]
    return new_prob

# pulling only scan numbers out
def extractScanNum(row):
    string = row
    spot = string.find('.')
    new_st = string[spot + 1:]
    spot = new_st.find('.')
    final_st = new_st[:spot]

    if final_st[0] == "0":
        final_st = final_st[1:]
    return final_st

def get_original_msf_data(file):
    msf_original_df = dl.clean_msfragger(file)

    # Extracting scan number from file number
    msf_original_df['scan'] = msf_original_df['scan'].apply(extractScanNum)

    msf_original_df = filter_data(msf_original_df, 'PeptideProphet Probability')

    # Changing the probabilities to the same scale all the other tools use
    msf_original_df["Updated_probability"] = msf_original_df.apply(set_probablility, 1)

    # filter out just the ScanNr and QValue and/or PEP columns
    msf_original_df = msf_original_df.filter(['scan', 'peptide', 'Updated_probability'])

    return msf_original_df


# Reading in and formatting the orginial MsgfPlus data
def get_original_msg_data(file):
    msg_original_df = dl.clean_msgfplus(file)
    msg_original_df = filter_data(msg_original_df, "QValue")

    # filter out just the ScanNr and QValue and/or PEP columns
    msg_original_df = msg_original_df.filter(['scan', 'peptide', "QValue"])

    return msg_original_df


# Reading in and formatting the orginial MaxQuant data
def get_original_mq_data(file):
    mq_original_df = dl.clean_maxquant(file)

    # Formatting and dropping any rows that are missing the sequence
    mq_original_df['Sequence'].replace(' ', np.nan, inplace=True)
    mq_original_df.dropna(subset=['Sequence'], inplace=True)

    mq_original_df = filter_data(mq_original_df, 'PEP')

    # filter out the ScanNr, peptide, and PEP columns
    mq_original_df = mq_original_df.filter(['scan', 'peptide', 'PEP'])

    return mq_original_df


# Reading in and formatting the data to be compared to the benchmarked data
def get_input_data(file, probability):
    df = pd.read_csv(file)

    df = filter_data(df, probability)

    # filter out the ScanNr, peptide, and probability columns
    df = df.filter(['scan', probability])

    return df

#Here we will read in the megascript that contains the output data from a single raw file that was ran through each tool.

# read in the megaScript and reformat it
def clean_meagScript(file):
    df = pd.read_csv(file, low_memory=False, header=[0, 1])
    df = df.drop("Unnamed: 0_level_0", axis = 1, level = 0)

    return df

#Slicing out the Peptide Prophet Probability values for MsFragger. There is no qvalue or PEP, so this is the row we are
# using. Counting how many are at or under the cutoff
def get_msf_prob_len(df, cutoff):
    msf_probability = df["MsFragger"]['Updated_probability']
    msf_probability =  msf_probability.dropna()
    msf_under_cutoff = len(msf_probability.loc[msf_probability <= cutoff])
    return msf_under_cutoff

#Slicing out the qvalues from MetaMorpheus and counting how many are at or under the cutoff
def get_mm_Qval_len(df, cutoff):
    mm_qval = df["MetaMorpheus"]["QValue"]
    mm_qval =  mm_qval.dropna()
    mm_under_cutoff = len(mm_qval.loc[mm_qval <= cutoff])
    return mm_under_cutoff

#Slicing out the PEP values from MetaMorpheus and counting how many are at or under the cutoff **Are we keeping this?
def get_mm_PEP_len(df, cutoff):
    mm_PEP = df["MetaMorpheus"]["PEP"]
    mm_PEP =  mm_PEP.dropna()
    value_under_cutoff = len(mm_PEP.loc[mm_PEP <= cutoff])
    return value_under_cutoff

#Slicing out the qvalues from MsgfPlus and counting how many are at or under the cutoff
def get_msg_Qval_len(df, cutoff):
    msg_qval = df["MsgfPlus"]["QValue"]
    msg_qval =  msg_qval.dropna()
    msg_under_cutoff = len(msg_qval.loc[msg_qval <= cutoff])
    return msg_under_cutoff

#Slicing out the PEP from MaxQuant. Counting how many are at or under the cutoff
def get_mq_PEP_len(df, cutoff):
    mq_PEP = df["MaxQuant"]["PEP"]
    mq_PEP =  mq_PEP.dropna()
    mq_under_cutoff = len(mq_PEP.loc[mq_PEP <= cutoff])
    return mq_under_cutoff

#Slicing out the probability column from the inputted data and counting how many scans are at or under the cutoff
def get_input_len(df, cutoff, probability, tool_name):
    df = df[tool_name][probability]
    df =  df.dropna()
    df_under_cutoff = len(df.loc[df <= cutoff])
    return df_under_cutoff

#This function gets the number of scan values that were at or below the cutoff for each tool and returns them.
def get_file_values(file, cutoff, inputs_probability, tool_name):
    df = clean_meagScript(file)
    msf = get_msf_prob_len(df, cutoff)
    MM_QVal = get_mm_Qval_len(df, cutoff)
    MM_PEP = get_mm_PEP_len(df, cutoff)
    msg_QVal = get_msg_Qval_len(df, cutoff)
    MQ_PEP = get_mq_PEP_len(df, cutoff)
    input_data = get_input_len(df, cutoff, inputs_probability, tool_name)
    values_list = {"msf" : msf, "MM_QVal" : MM_QVal, "MM_PEP" : MM_PEP, "msg_QVal" : msg_QVal, "MQ_PEP" : MQ_PEP, 'input_data':                      input_data}
    return values_list

#Reading in the data and making the graph for the 2ng data at a certain cutoff
def make_2ng_graph(input_probability, tool_name, cutoff = 0.01):
    File1 = get_file_values("benchmark_MegaScript_2ng_rep1.csv", cutoff, input_probability, tool_name)
    File2 = get_file_values("benchmark_MegaScript_2ng_rep2.csv", cutoff, input_probability, tool_name)
    File3 = get_file_values("benchmark_MegaScript_2ng_rep3.csv", cutoff, input_probability, tool_name)
    File4 = get_file_values("benchmark_MegaScript_2ng_rep4.csv", cutoff, input_probability, tool_name)
    File5 = get_file_values("benchmark_MegaScript_2ng_rep5.csv", cutoff, input_probability, tool_name)
    File6 = get_file_values("benchmark_MegaScript_2ng_rep6.csv", cutoff, input_probability, tool_name)

    # set width of bars
    barWidth = 0.14

    # set heights of bars
    msf_prob = [File1['msf'], File2['msf'], File3['msf'], File4['msf'], File5['msf'], File6['msf']]
    MM_PEP = [File1['MM_PEP'], File2['MM_PEP'], File3['MM_PEP'], File4['MM_PEP'], File5['MM_PEP'], File6['MM_PEP']]
    MM_qval = [File1['MM_QVal'], File2['MM_QVal'], File3['MM_QVal'], File4['MM_QVal'], File5['MM_QVal'],
               File6['MM_QVal']]
    msg_qval = [File1['msg_QVal'], File2['msg_QVal'], File3['msg_QVal'], File4['msg_QVal'], File5['msg_QVal'],
                File6['msg_QVal']]
    mq_PEP = [File1['MQ_PEP'], File2['MQ_PEP'], File3['MQ_PEP'], File4['MQ_PEP'], File5['MQ_PEP'], File6['MQ_PEP']]
    input_prob = [File1['input_data'], File2['input_data'], File3['input_data'], File4['input_data'],
                  File5['input_data'], File6['input_data']]

    # Set position of bar on X axis
    r1 = np.arange(len(msf_prob))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]
    r4 = [x + barWidth for x in r3]
    r5 = [x + barWidth for x in r4]
    r6 = [x + barWidth for x in r5]

    # Make the plot
    plt.bar(r1, msf_prob, width=barWidth, edgecolor='white', label='MsFragger Peptide Prophet Probability')
    plt.bar(r2, MM_qval, width=barWidth, edgecolor='white', label='MetaMorpheus Q-Value')
    plt.bar(r3, msg_qval, width=barWidth, edgecolor='white', label='MsgfPlus Q-Value')
    plt.bar(r4, mq_PEP, width=barWidth, edgecolor='white', label='MaxQuant PEP')
    plt.bar(r5, MM_PEP, width=barWidth, edgecolor='white', label='MetaMorpheus PEP')
    plt.bar(r6, input_prob, width=barWidth, edgecolor='white', label= tool_name + " " + input_probability)

    # Add xticks on the middle of the group bars
    plt.ylabel('Number of PSMs found')
    plt.xlabel(tool_name + " compared to the benchmarked 2ng data")
    plt.title('2ng')
    plt.xticks([r + barWidth for r in range(len(msf_prob))], ['File1', 'File2', 'File3', 'File4', 'File5', 'File6'])

    # Create legend & Show graph
    plt.legend(loc="upper right", bbox_to_anchor=(1.73, 1))
    plt.show()
    #plt.savefig('2ng_native_score.png')

#Reading in the data and making the graph for the 0.2ng data at a certain cutoff
def make_02ng_graph(input_probability, tool_name, cutoff = 0.01):
    File1 = get_file_values("benchmark_MegaScript_0.2ng_rep1.csv", cutoff, input_probability, tool_name)
    File2 = get_file_values("benchmark_MegaScript_0.2ng_rep2.csv", cutoff, input_probability, tool_name)
    File3 = get_file_values("benchmark_MegaScript_0.2ng_rep3.csv", cutoff, input_probability, tool_name)
    File4 = get_file_values("benchmark_MegaScript_0.2ng_rep4.csv", cutoff, input_probability, tool_name)
    File5 = get_file_values("benchmark_MegaScript_0.2ng_rep5.csv", cutoff, input_probability, tool_name)
    File6 = get_file_values("benchmark_MegaScript_0.2ng_rep6.csv", cutoff, input_probability, tool_name)

    # set width of bars
    barWidth = 0.14

    # set heights of bars
    msf_prob = [File1['msf'], File2['msf'], File3['msf'], File4['msf'], File5['msf'], File6['msf']]
    MM_PEP = [File1['MM_PEP'], File2['MM_PEP'], File3['MM_PEP'], File4['MM_PEP'], File5['MM_PEP'], File6['MM_PEP']]
    MM_qval = [File1['MM_QVal'], File2['MM_QVal'], File3['MM_QVal'], File4['MM_QVal'], File5['MM_QVal'],
               File6['MM_QVal']]
    msg_qval = [File1['msg_QVal'], File2['msg_QVal'], File3['msg_QVal'], File4['msg_QVal'], File5['msg_QVal'],
                File6['msg_QVal']]
    mq_PEP = [File1['MQ_PEP'], File2['MQ_PEP'], File3['MQ_PEP'], File4['MQ_PEP'], File5['MQ_PEP'], File6['MQ_PEP']]
    input_prob = [File1['input_data'], File2['input_data'], File3['input_data'], File4['input_data'],
                  File5['input_data'], File6['input_data']]

    # Set position of bar on X axis
    r1 = np.arange(len(msf_prob))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]
    r4 = [x + barWidth for x in r3]
    r5 = [x + barWidth for x in r4]
    r6 = [x + barWidth for x in r5]

    # Make the plot
    plt.bar(r1, msf_prob, width=barWidth, edgecolor='white', label='MsFragger Peptide Prophet Probability')
    plt.bar(r2, MM_qval, width=barWidth, edgecolor='white', label='MetaMorpheus Q-Value')
    plt.bar(r3, msg_qval, width=barWidth, edgecolor='white', label='MsgfPlus Q-Value')
    plt.bar(r4, mq_PEP, width=barWidth, edgecolor='white', label='MaxQuant PEP')
    plt.bar(r5, MM_PEP, width=barWidth, edgecolor='white', label='MetaMorpheus PEP')
    plt.bar(r6, input_prob, width=barWidth, edgecolor='white', label= tool_name + " " + input_probability)

    # Add xticks on the middle of the group bars
    plt.ylabel('Number of PSMs found')
    plt.xlabel(tool_name + " compared to the benchmarked 0.2ng data")
    plt.title('0.2ng')
    plt.xticks([r + barWidth for r in range(len(msf_prob))], ['File1', 'File2', 'File3', 'File4', 'File5', 'File6'])

    # Create legend & Show graph
    plt.legend(loc="upper right", bbox_to_anchor=(1.73, 1))
    plt.show()
    #plt.savefig('0.2ng_native_score.png')


#names of all the files we are going to read in to upload our data
data_files = ["2ng_rep1", "2ng_rep2", "2ng_rep3", "2ng_rep4", "2ng_rep5", "2ng_rep6",
             "0.2ng_rep1", "0.2ng_rep2", "0.2ng_rep3", "0.2ng_rep4", "0.2ng_rep5", "0.2ng_rep6"]

#This is here to help us with naming out output files
input_names = {} 
    #2ng files
input_names["2ng_rep1"] = "Ex_Auto_J3_30umTB_2ngQC_60m_1"
input_names["2ng_rep2"] = "Ex_Auto_J3_30umTB_2ngQC_60m_2"
input_names["2ng_rep3"] = "Ex_Auto_K13_30umTA_2ngQC_60m_1"
input_names["2ng_rep4"] = "Ex_Auto_K13_30umTA_2ngQC_60m_2"
input_names["2ng_rep5"] = "Ex_Auto_W17_30umTA_2ngQC_60m_3"
input_names["2ng_rep6"] = "Ex_Auto_W17_30umTA_2ngQC_60m_4"

    #0.2ng files
input_names["0.2ng_rep1"] = "Ex_Auto_J3_30umTB_02ngQC_60m_1"
input_names["0.2ng_rep2"] = "Ex_Auto_J3_30umTB_02ngQC_60m_2"
input_names["0.2ng_rep3"] = "Ex_Auto_K13_30umTA_02ngQC_60m_1"
input_names["0.2ng_rep4"] = "Ex_Auto_K13_30umTA_02ngQC_60m_2"
input_names["0.2ng_rep5"] = "Ex_Auto_W17_30umTA_02ngQC_60m_3"
input_names["0.2ng_rep6"] = "Ex_Auto_W17_30umTA_02ngQC_60m_4"

#Our data is split into 12 output files. We will read in the data one output file at a time for each tool and make a
# megascript. The megascript allows us to look at the output for each tool from that specific raw input file.
#We begin by reading in the output file from each tool. These are the output files that each tool gives us from running
# a specific raw file. We will then set the index to the scan column for each of the individual dataframes.

def make_megascript(input_files, tool_name):
    # This gets the absolute path to the location of the data_loader.py file
    path_to_data_loader = os.path.abspath(os.path.dirname(__file__))  
    
    
    for file in data_files:
        mm_df = get_orgininal_mm_data(file)
        msf_df = get_original_msf_data(file)
        msg_df = get_original_msg_data(file)
        mq_df = get_original_mq_data(file)
        input_data = input_files[file]

        # Switching index to ScanNr to join dataframes based on their scan numbers.
        MsgfPlus = msg_df.set_index('scan')
        MsFragger = msf_df.set_index('scan')
        MetaMorpheus = mm_df.set_index('scan')
        MaxQuant = mq_df.set_index('scan')
        InputData = input_data.set_index('scan')

        # concating all the individual joined dataframes together to make one megascript
        megaScript = pd.concat(dict(MsFragger=MsFragger, MsgfPlus=MsgfPlus, MetaMorpheus=MetaMorpheus,
                                    MaxQuant=MaxQuant, InputData = InputData), axis=1)
        megaScript.reset_index(inplace=True)
        
        megaScript = megaScript.rename(columns={'InputData' : tool_name})

        # saving the megascript
        megaScript.to_csv("benchmark_MegaScript_" + file + ".csv")