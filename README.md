# R-Shiny-Application
The following contains code to create an R Shiny application intended for use in analyzing mRNA-seq data. Data was derived from GEO accession GSE64810 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810). 

The application allows users to upload three csv files, sample information/metadata, counts matrix, and deseq2 results. The RShiny application contains four main tabs for viewing and analyzing the uploaded data: Samples, Counts, Differential Gene Expression, and Visualize Gene Expression. 

The samples tab provides an overview of the types of information available and some summary statistics for the samples. It allows the user to view the file in a table format and to visualize graphically the distribution of some sample characteristics across sample groups.

The counts tab allows the user to filter the counts matrix by percent variance and number of non-zero samples in each row. A summary table is provided to indicate how the filtering affects the counts matrix and diagnostic scatterplots are also available to visualize this graphically. A heatmap is also shown to display results for individual genes after filtering. Finally, a Principle Component Analysis is available and the user can adjust what components are displayed on the x and y axes. 


