This FERRENBERG_ETAL__DATA_README.txt file was generated on 2022-09-11 by Scott Ferrenberg

GENERAL INFORMATION

1. Title of Dataset: Data from: Divergent growth-differentiation balance strategies and resource competition shape mortality patterns in ponderosa pine

2. Author Information
	Corresponding Investigator 
		Name: Dr Scott Ferrenberg
		Institution: University of Montana, Missoula, MT, USA; New Mexico State University, Las Cruces, NM, USA
		Email: ferrenbe@nmsu.edu

3. Date of data collection: 2018-2020

4. Geographic location of data collection: Gila National Forest/surrounding area, New Mexico, USA

5. Funding sources that supported the collection of the data: NONE

6. Recommended citation for this dataset: Ferrenberg, Scott, et al. (2022), Data from: Divergent growth-differentiation balance strategies and resource competition shape mortality patterns in ponderosa pine, Dryad, Dataset


DATA & FILE OVERVIEW

1. Description of dataset

These data were generated to quantify the growth-defense relationships of Pinus ponderosa that were killed by bark beetles vs. those that survived the insect epidemic. The data are in 4 .csv files and include measures of tree growth, resin duct production, resin duct sizes, and tree attributes: measures of tree age (years), radius covered by the increment core (mm), diameter at 1.4 m above the ground (cm), height (m), and competition neighborhood density (count) and basal area (cm) at the time of sampling. Competition measures are subset into intra- and inter-specific values. 
See associated manuscript in Ecosphere for complete methods and data collection details.


2. File List: 
	File 1 Name: Tree_attributes.csv
	File 1 Description: Tree ID and status, tree age (years), radius (mm), diameter at 1.4 m above the ground (cm), height (m), and competition neighborhood density (count) and basal area (cm) at the time of sampling.

	File 2 Name: Gila_PIPO_growth.1.csv
	File 2 Description: Annual increment growth (mm) for paired rings of all trees. This file also contains an annual "lag" value for helping with data alignment. Column headers indicate the ring pair, time lag, and tree IDs (L or D, for live/dead at the time of sampling, and a pair number).
	
	File 3 Name:  Gila_PIPO_RDcount.1.csv
	File 3 Description: Counts of resin ducts in annual increment growth rings. Column headers include ring pair (1 to 29 years indicating the most recent to historical rings in each pair of trees) and the tree IDs.

	File 4 Name: Gila_PIPO_RDarea.1.csv
	File 4 Description: Cross-sectional area (mm^2) of all resin ducts in annual increment growth rings for all trees. Column headers include ring pair (1 to 29 years indicating the most recent to historical rings in each pair of trees)and the tree IDs.

	
METHODOLOGICAL INFORMATION

See associated manuscript in Ecosphere for complete methods and data collection details. Tree ring data were measured from 12 mm increment cores collected at 1.4 m above the ground surface. 
Tree attribute data were collected in the field setting; tree height was measured with a hypsometer, diameter with a standard DBH tape, and competitor density was measured from plot established around focal trees.


DATA-SPECIFIC INFORMATION FOR: Tree_attributes.csv

1. Number of variables: 12

2. Number of cases/rows: 77 (this file is in a "long" format)

3. Variable List: 
	tree: a tree ID value for each entry that corresponds to columns in other data files.
	age: a numerical interger value indicating the age of the tree in years at the time (calendar year) when bark beetles killed the "dead" tree in each pair of live/dead trees.
	rad_end: a numerical continuous value that indicates the radius of xylem from the pith to the last annual growth ring in the time series derived from increment cores.
	DBH: a numerical continuous value that indicates the diameter of the tree at 1.4 m above the ground surface; note, this value includes all annual growth increments and the bark so in cases where the "live" tree out-lived the "dead" tree in a pair it could include many additional growth rings.
	HT: a numerical continous value that describes each trees' height in meters.
	comp_tot: a numerical integer value indicating the total number of "competitor" trees greater than 1.4 m in height within the local neighborhood of focal trees.
	comp_ba: a numerical continuous value indicating the cummulative basal area of "competitor" trees greater than 1.4 m in height within the local neighborhood of focal trees.
	comp_intra: a numerical integer value indicating the total number of intraspecific "competitor" trees greater than 1.4 m in height within the local neighborhood of focal trees.
	intra_ba: a numerical continuous value indicating the cummulative basal area of intraspecific "competitor" trees greater than 1.4 m in height within the local neighborhood of focal trees.
	comp_inter: a numerical integer value indicating the total number of interspecific "competitor" trees greater than 1.4 m in height within the local neighborhood of focal trees.
	inter_ba: a numerical continuous value indicating the cummulative basal area of interspecific "competitor" trees greater than 1.4 m in height within the local neighborhood of focal trees.
	mean_intra_dbh: a numerical continuous value indicating the mean diameter at 1.4 m above the ground surface intraspecific "competitor" trees greater than 1.4 m in height within the local neighborhood of focal trees.

4. Missing data codes: 

5. Abbreviations used: 
	N/A; not applicable

6. Other relevant information: 
	

DATA-SPECIFIC INFORMATION FOR: Gila_PIPO_growth.1.csv

1. Number of variables: 3

2. Number of cases/rows: 29 (this file is in a "wide" format)

3. Variable List:
	Ring_pair: a numerical integer value that aligns pairs of rings among live "L" and dead "D" trees
	Lag: a numerical integer value that describes the time lag within the ring series.
	Growth values: Column headers indicate tree ID; columns contain continous numerical values in wide-data format that indicate the annual increment growth of each tree for an aligned pair of 29 years ending with the last year of growth in the tree killed by bark beetles in each pair, also known as the dead "D" tree for simplicity.
	
	
4. Missing data codes: 'NA' for not available, meaning no growth was measured because of damage or other issues.
	
5. Abbreviations used: "L" for live tree and "D" for dead tree.
	
6. Other relevant information: Numbers in headers indicate tree pairs.
	

DATA-SPECIFIC INFORMATION FOR: Gila_PIPO_RDcount.1.csv

1. Number of variables: 2

2. Number of cases/rows: 29

3. Variable List: 
	Ring_pair: a numerical interger value that aligns pairs of rings among live "L" and dead "D" trees
	Resin_duct_counts: Column headers indicate tree ID; columns contain numerical integer values in wide-data format that indicate the number of resin ducts present in increment growth rings of each tree for an aligned pair of 29 years ending with the last year of growth in the tree killed by bark beetles in each pair, also known as the dead "D" tree for simplicity.
	
4. Missing data codes: 'NA' for not available, meaning no resin ducts were counted because of damage or other issues.

5. Abbreviations used: "L" for live tree and "D" for dead tree.

6. Other relevant information: Numbers in headers indicate tree pairs, the 'rd' suffix indicates "resin ducts" and is included to identify the variable type for those combining data using code.
	

DATA-SPECIFIC INFORMATION FOR: Gila_PIPO_RDarea.1.csv
1. Number of variables: 2

2. Number of cases/rows: 29

3. Variable List: 
	Ring_pair: a numerical interger value that aligns pairs of rings among live "L" and dead "D" trees
	Resin_duct_area: Column headers indicate tree ID; columns contain numerical continous values in wide-data format that indicate the cumulative cross sectional area of resin ducts present in increment growth rings of each tree for an aligned pair of 29 years ending with the last year of growth in the tree killed by bark beetles in each pair, also known as the dead "D" tree for simplicity.
	
4. Missing data codes: 'NA' for not available, meaning no resin ducts were counted because of damage or other issues.

5. Abbreviations used: "L" for live tree and "D" for dead tree.
	
6. Other relevant information: Numbers in headers indicate tree pairs, the 'da' suffix indicates "duct area" and is included to identify the variable type for those combining data using code.
	
	
