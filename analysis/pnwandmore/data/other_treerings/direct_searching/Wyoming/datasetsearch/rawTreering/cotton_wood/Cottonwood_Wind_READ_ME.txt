Contents of File Cottonwood_Wind_Tree_Data.csv and Cottonwood_Wind_Ring_Widths


Explanation of fields in File Cottonwood_Wind_Tree_Data.csv

This is a csv file with labelled fields:
	
tree		The name of each tree included in the analysis. The first two letters designate the study reach. "BO" are trees near Boysen Reservoir, "LD" are trees downstream of Diversion Dam, and OW are trees near Owl Creek that were not included in the analysis. The third letter refers to any special status of the tree ("C" refers to a cored, living tree that was part of the random sample, "D" refers to a sampled dead tree,  "Q" refers to a cored live tree that was not part of the random sample, and “L” refers to nonrandom living trees cored at Owl Creek. The last two to three numbers are the specific ID of the tree in its reach. 

date		Date tree was sampled, MM/DD/YYYY. Unknown dates reported as "UNK".

reach		Location where tree is found. Values include: "boysen" (trees near Boysen Reservoir), "boysen east" (dead trees collected in the eastern part of the Boysen study area, "boysen west" (dead trees collected in the western part of the Boysen study area", and "lower diversion dam" (trees downstram of Diversion Dam).

random		If tree was part of the random sample it is given the value "r", otherwise "nr".

species		When possible to determine, the species of cottonwood was listed. "PODE" is Populus deltoides, "POAN" is Populus angustifolia, "Hybrid" is a hybrid cottonwood (deltoides X angustifolia) and "UNK" means species unknown.

center year	Estimate of the center year of the tree. Some trees included in the analysis could not have a center year estimated with confidence; center years of these trees are reported as "UNK".

dbh		Diameter at 1.2 m above ground, measured in centimeters. Unknown dbh reported as “UNK”.

height		Height of tree in meters. Unknown height reported as "UNK".

vigor		Percentage of live tree canopy, from 0 (completely dead) to 100 (full living canopy) in increments of 5. Unknown vigor reported as "UNK".

status		Notes if tree was living ("L"), standing dead ("S"), fallen dead ("F"), or previously cut down ("C"; not included in death dates). 
	

Explanation of File Cottonwood_Wind_Ring_Widths

This is a text file in Tucson format with measurements in microns (thousandths of millimeters).

Column 1 (7 characters total) includes the tree name (5 to 6 characters) and the core name (A or B).

Column 2 (4 characters) is the year of the ring width given in Column 3.

Columns 3-12 give widths (in thousandths of millimeters) of rings formed in successive years beginning with the year given in column 2. End of ring width series is delimited by -9999.  
