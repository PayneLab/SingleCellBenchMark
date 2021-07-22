# SingleCellBenchMark

## Introduction

All data from the peptide identification tools are sorted in the data folder, which are parsed and loaded into pandas dataframes by data_loader.py. To make Figures 1 and 2, data is aggregated across all tools and stored in a table called megaScript. 

## Running algorithms
### pepetide identifcation
The 12 files were run through each algorithm and the output files for each tool are stored in the /data folder  “~/data/”. The parameters for each tool are stored here. 
### Mokapot
The psm file from each tool was run through mokapot and the results were sotred in '~/MokaPot/MokaPot_output'.

### Making the megascript
First, all PSM files are read in via data_loader.py that formats the data to make it joinable. For each tool, the scan number, peptide, and statistical score columns were kept. The mokapot output files are also joined to megaScript keeping scan number, mokapot score, and mokapot q-value columns. Because Proteome Discoverer automatically runs Percolator, it has a percolator q-value that is included.  12 megascripts are created, one for each replicate. 

The code to make the megascripts is also found in the repository at “~MokaPot/make_MegaScript.ipynb”.

## Making figure 1 and 2
After the mega-scripts are made, we use them to produce Figure1 and Figure2. 
To make Figure1, for each mega-script, we slice out the different tools and their native statistical columns. The number of PSMs with a statistical score at or below 0.01 were then counted and totalled. The exact code for making Figure 1 can be found in the file “~/makeFigure1.ipynb”.

Figure2 is calculated in the same fashion by slicing and counting the number of PSMs with a mokapot q-value percolator q-value under the 0.01 cutoff. Because Proteome Discoverer automatically runs Percolator, it is only comparble to mokapot results, and as such is included in Figure 2 and not Figure 1. 

For both figures, the scores from the 2ng and the 0.2ng data are combined and compared in separate graphs, leading to a part A and B for both Figure1 and Figure2. 
The exact code for making Figure 2 can be found in the file "~/makeFigure2.ipynb". 

## Optimizations
Three new columns were added onto the PSM identifications of MetaMorpheus as additional features for mokapot to use in its rescoring: number of consecutive y peaks, the percent of annotated peaks in a consecutive series, and the absolute difference between predicted retention time and actual retention time. These new features were joined with the MetaMorpheus features (pin file). Mokapot was run in the same manner as above, with the exception that we ran the program 25 times and recorded the output each time to look for variability in the non-deterministic optimization. The code to add the new features and calculations can be found in the repository at “~/Optimization” and "~/updated_features.ipynb”. 

Feature 1. The longest consecutive set of annotated y peaks was found parsing the MetaMorpheus output column “Matched Ion Series” and reporting the longest ladder of consecutive y peaks. For example if the y peaks in a ladder are [1,2,3,6,7] then the longest consecutive peak count would be 3. See ~/Optimization/calculate_consecutive_peaks.ipynb for the full code. 

Feature 2. The second feature is the percent of annotated y peaks which are part of a consecutive series. This is calculated by taking the number of consecutive peaks in a ladder and dividing it by the length of the ladder. For example, if the annotated y peaks are [1,5,6] the number of consecutive peaks is 2 and total annotated peaks is 3. This feature would be 66%. The code to calculate this feature can be found in the GitHub repository at “~/Optimization/calculate_consecutive_peaks.ipynb”. 

Feature 3. The final new feature was the absolute difference between observed and predicted retention time from AutoRT16. The training file can be found in the GitHub repository at “/Optimization/RT_training”. The code to produce this file is found at "Optimization/retention_feature_prep.ipynb". AutoRT was trained by taking all 10,138 scans from one file with a q-value = 0.0, after decoys and duplicated scans were dropped (see “/Optimization/RT_training.tsv”). The model was then trained using recommended settings. The training model given back from AutoRT is found at “/PXD006109_prediction/test.tsv”. Predictions using the model can be found in the repository at “~/Optimization/predicted_data”.

Once the new features were calculated they were joined with the original pin file from the MetaMorphues output and run through mokapot. The code to add the new features, rerun the file through mokapot, and the calculations, can be found in the repository at “~/updated_features.ipynb”. 
