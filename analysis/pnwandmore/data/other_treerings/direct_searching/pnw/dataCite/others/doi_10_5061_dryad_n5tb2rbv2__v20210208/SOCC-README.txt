For each of the directories listed below input and/or output files are provided for each analysis described in the paper.

Admixture
Input/ - The input file is called 'populations_r20.haplotypes.filtered_m70_randomSNP_recoded.ped'. It was used to run admixture using k values 1-10.


DAPC
Input/ - The input file for running DAPC using the R package ADEGENET is called 'populations_r20.haplotypes.filtered_m70_randomSNP.str'. It was used to run the `find.clusters` function with max.n.clust=10, n.pca = 100, choose.n.clust = FALSE.


Dadi
The population key file is called 'Population-Key-Full.txt'. Population abbreviations are as follows: West Sierra Nevada (WSN); Pacific Northwest (PNW); East Sierra Nevada (ESN); Great Basin (GB); Southern California (SCA) 
Inputs/ - The input files for running dadi are in .tsv format. The input file for the full dataset is called 'populations_r70.haplotypes.filtered_m50_randomSNP.tsv'. The inputs for pairwise population comparisons are:
    SCA_ESNGB-populations_r50.haplotypes.filtered_m70_randomSNP.tsv
    WSN_PNW-populations_r50.haplotypes.filtered_m70_randomSNP.tsv
    WSN_ESN-populations_r50.haplotypes.filtered_m70_randomSNP.tsv
    PNW_ESN-populations_r50.haplotypes.filtered_m70_randomSNP.tsv
    ESN_GB-populations_r50.haplotypes.filtered_m70_randomSNP.tsv
Logs/ - Log files from running the analyses.


EEMS
Inputs/ - The input files are called 'Scelop.diffs', 'Scelop.outer', 'Scelop.coord'. They were used to run EEMS with the included .ini files, called '700-Scelop-chain*.ini'.


SNAPP
Input/ - The alignment file used to create the SNAPP input is called 'populations_r20.haplotypes.filtered_m70_randomSNP_SNAPP_downsampled.nex'. The .xml file used at the SNAPP input is called 'socc_snapp_input.xml'.


SplitsTree
Input/ - The formatted alignment file used as an input for SplitsTree is called 'populations_r20.haplotypes.filtered_m70_randomSNP.splitstree.nex'. 


Stacks_pipeline
The input files for Stacks are hosted on the NCBI SRA (SRA PRJNA616381).
Logs/ - The log files from running the analyses.
Outputs/ - Formatted output files from running Stacks_pipeline that were used as inputs for downstream analyses.


TreeMix
Input/ - The input file is called 'SOCC_TreeMix_input_filtered.txt.gz'. TreeMix was run as described in the Methods.