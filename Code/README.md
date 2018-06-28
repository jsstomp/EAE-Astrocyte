# Code #

* This part of the repository consists of scripts used during this study.

### Scripts used are described here: ###

|Script                  |Description                                                                                                                                                               |
|------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|__prepare_expression.R__|R script called from the commandline, it will initiate *DAFS.R* and *PCA_plotter.R* if required by command line arguments                                                 |
|DAFS.R                  |R script that filter lowly expressed genes out of the count data                                                                                                          |
|PCA_plotter.R           |R script that plots PCA results based on gene expression levels of filtered count data                                                                                    |
|__DEA_DESeq2.R__        |R script called from the commandline, that starts differential expression analysis, also functions as main for the folowing scripts which can be activated as arguments   |
|DEA_analysis.R          |R script that creates visualization of DEA results in the form of heatmaps for all results as well as for different regions and populations                               |
|venn_DEA.R              |R script that creates Venn diagrams portraying the amount of signature genes for the different regions and populations                                                    |
|GOI.R                   |R script that enriches DEA results for GO terms and makes a top 20 barplot for each comparison done in the DEA                                                            |
|GO_grouper.R            |R script that creates groups of GO terms based on semantic similarity, uses *circos_plotter.R* for visualizing the results compared to comparisons of DEA                 |
|circos_plotter.R        |R script that visualizes the results compared to comparisons of DEA in cicular plots that are easy to read                                                                |
|WGCNA.R                 |R script that performs WGCNA, finds gene clusters based on gene expression values and enriches them for different pathways                                                |

* Scripts in bold can be run via command line, use --help for more options.