# Protein-Schema
The programs and data files in this repository are for the article "A percolation theory analysis of continuous functional paths in protein sequence space affirms previous insights on the optimization of proteins for adaptability": https://doi.org/10.1371/journal.pone.0314929.


The program **schema_analysis.py** reads the FireProtDB data file **Data Files/fireprotdb_results.csv**, and it imports functions from **schema_functions.py**. The included output file used FireProtDB data, last updated on Feb 8, 2022. The latest data can be downloaded from https://loschmidt.chemi.muni.cz/fireprotdb/.The program then cleans and analyzes the data, saving the output to the Excel file **Data Files/mut_out 30.xlsx**. If the variable **SAMPLE** is changed from 30 to another number, the required number of entries for a protein in FireProtDB to be analyzed will change to the new value. The name of the output file will also reflect the change. 

The output file contains the **Samples** worksheet, which lists the percentage of destabilizing mutations for each protein for both the 0.5 kcal/mol and 1.0 kcal/mol cutoffs. The file also contains the **Combined** worksheet, which lists the percentage of discoverable proteins for the different cutoffs. The program also reads the Tawfik protein data from the file **Data Files/mut_taw.xls** and saves the analysis output in the worksheet **Tawfik Data**. 

In addition, **schema_analysis.py** creates several figures. The program **schema_figs.py** creates additional figures, including the FSH plots. The program **seq_space.py** converts target c-values in nucleotide sequence space to amino acid sequence space. 
