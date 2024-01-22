# R-Shiny-Application
The following contains code to create an R Shiny application intended for use in analyzing mRNA-seq data. Data was derived from GEO accession GSE64810 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810). 

The application allows users to upload three csv files, sample information/metadata, counts matrix, and deseq2 results. The RShiny application contains four main tabs for viewing and analyzing the uploaded data: Samples, Counts, Differential Gene Expression, and Visualize Gene Expression. 

The samples tab provides an overview of the types of information available and some summary statistics for the samples. It allows the user to view the file in a table format and to visualize graphically the distribution of some sample characteristics across sample groups.

The counts tab allows the user to filter the counts matrix by percent variance and number of non-zero samples in each row. A summary table is provided to indicate how the filtering affects the counts matrix and diagnostic scatterplots are also available to visualize this graphically. A heatmap is also shown to display results for individual genes after filtering. Finally, a Principle Component Analysis is available and the user can adjust what components are displayed on the x and y axes. 

The differential gene expression tab allows the user to visualize the deseq2 results in a table format with the ability to search for specific values and sort column values. The results can also be viewed in a volcano plot format where the user can select the y and x axes and color results based on their selected p-adjusted magnitude. 

The visualize gene expression tab gives the user the ability to select a gene of interest and visualize the expression results by sample categories. The user can select from several different plot type options to visualize the expression data. 
