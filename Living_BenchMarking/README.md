# Running the Living Benchmarking Program:

In order to run your program's output to our benchmarked data you will use the Living_Benchmarking_Template.ipynb notebook.
You will first need to run the same data through your program. Begin by downloading six 2ng and six 0.2ng files from MassIVE (MSV000087524). The benchmarking program will be working with the peptides, so you will need find the output file that gives the peptides. The instructions in the notebook will direct you to create the .csv files that have a "scan", "decoy", "peptide", and a probability column. Follow these instructions to create files in the correct format to run through the benchmarking program. Once you have the files in the correct format, you will just need to continue to follow the instructions to run the files through the benchmarking program. 

When you are done, two graphs will save to your computer as '2ng_benchmark_graph.png' and '0.2ng_benchmark_graph.png'. This will show how many PSMs at or below a 1% FDR cutoff each tool, including yours, is finding for the relative raw files.