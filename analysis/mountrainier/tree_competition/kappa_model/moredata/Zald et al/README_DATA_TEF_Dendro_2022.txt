This README_DATA_TEF_Dendro_2022.txt file was generated on March 15, 2022 by Harold Zald


GENERAL INFORMATION

1. Title of Dataset: Data from: Tree growth responses to extreme drought after mechanical thinning and prescribed fire in a Sierra Nevada mixed-conifer forest, USA. 

2. Author Information
	Corresponding Investigator 
		Name: Dr Harold Zald
		Institution: USDA Forest Service, Corvallis, Oregon, USA
		Email: harold.zald@usda.gov

	Co-investigator 1
		Name: Chance Callahan 
		Institution: Cal Poly Humboldt, Arcata, California, USA

	Co-investigator 2
		Name: Dr Matthew Hurteau
		Institution: University of New Mexico, Albuquerque, New Mexico, USA

	Co-investigator 3
		Name: Marissa Goodwin
		Institution: University of New Mexico, Albuquerque, New Mexico, USA

	Co-investigator 4
		Name: Dr Malcolm North
		Institution: USDA Forest Service, Mammoth Lakes, California, USA
	
3. Date of data collection: 2011-2017

4. Geographic location of data collection: Teakettle Experimental Forest, Sierra National Forest, California, USA.

5. Funding sources that supported the collection of the data: USDI Joint Fire Science Program (Project ID 15-1-07-6), California Department of Forestry and Fire Protection as part of the California Climate Investments Program (Grant #8GG14803).

6. Recommended citation for this dataset: TEF_Dendro, Zald et al. (2022), Data from: Tree growth responses to extreme drought after mechanical thinning and prescribed fire in a Sierra Nevada mixed-conifer forest, USA, Dryad, Dataset.

DATA & FILE OVERVIEW

1. Description of dataset

These data were generated to examined tree growth responses after thinning, prescribed burning, and extreme drought at the Teakettle Experimental Forest (TEF), a historically frequent fire mixed-conifer forest in the southern Sierra Nevada of California, USA. Growth responses to thinning and prescribed burning (individually and in combination) were quantified using annual ring width data from 1440 increment cores collected from 720 trees (two cores collected per sampled tree). Annual ring width series (tk.rwl) were cross-dated to ensure correct calendar year assignment of ring-widths. Series that could not be cross-dated due to rotten or fragmented cores were discarded, resulting in cross-dated series for 1,401 of the 1,440 cores collected (713 of 720 trees sampled). Individual tree data (tk.lookup) includes: the tree tag (numbered aluminum tag nailed to tree), experimental unit (plot), four letter species code, secondary tree identifier, date of increment core collection, tree diameter in centimeters at breast height (DBH, 1.31 m) collected in 2017, tree height in meters, live crown ratio as a percentage, categorical variable of increment core type, sampling strata each tree is in, and unique increment core identifiers. Competition and topographic data (tk.comp.topo) are derived from the 2011 remeasurement of the stem map of 15,508 trees at TEF and GIS analysis of lidar detection and ranging (lidar) data collected in 2010. The tk.comp.topo data includes: the tree identifier, experimental unit (plot), four letter species code, categorical identifier of tree cored status, DBH in 2017 and 2011, basal area of each tree in 2011, heights of cored trees, live crown ratio of cored trees, number of trees within 10 m of each tree, number of trees less than 25 cm DBH within 10 m of each tree, number of trees greater or equal to 25 cm DBH within 10 m of each tree, total basal area of all trees within 10 m of each tree, total basal area of all trees less than 25 cm DBH within 10 m of each tree, total basal area of all trees greater or equal to 25 cm DBH within 10 m of each tree, potential solar radiation and topographic wetness index were derived from lidar data collected in 2011, and are the average values from 1 m rasters extracted in a 10 m radius around each tree. 

lidar data and derived GIS rasters have not been included in the dataset due to their large size. Requests to access these data should be made to the Harold Zald (harold.zald@usda.gov).

2. File List: 
	File 1 Name: tk.rwl
	File 1 Description: R dataframe of cross-dated annual ring width series from increment cores

	File 2 Name: tk.lookup
	File 2 Description: R dataframe of attributes for cored trees

	File 3 Name: tk.comp.topo
	File 3 Description: R dataframe of tree attributes, competition metrics, and topographic indcies for all stem mapped trees measured in 2011

METHODOLOGICAL INFORMATION

The study was conducted within the 1300 ha Teakettle Experimental Forest (TEF), located in the Sierra National Forest, approximately 80 km east of Fresno, California, USA. In 1998, 18 permanent 4 ha treatment units were established in a factorial design, with two levels of prescribed burning (no burn and fall burn) and three levels of thinning (no thin, understory thin, and overstory thin), for a total of six treatment combinations. Three replicate treatment units were assigned to each treatment combination, with thinning randomly assigned, while burn units were assigned with restricted randomization due to fire line and containment considerations. Thin and burn treatments were thinned in 2000 and burned in 2001, while thin-only treatments were thinned in 2001. Understory thinning removed trees 25-76 cm in diameter while retaining at least 40% of pre-treatment tree canopy cover. Overstory thinning removed trees greater than 25 cm in diameter, while retaining approximately 22 regularly spaced large diameter trees (generally >100 cm) per hectare.

Prior to treatment implementation (1998-2000 for treated plots, 2001-2002 for control plots) a complete census was conducted of all trees and snags > 5 cm diameter at breast height (DBH) within the treatment units. Trees and snags were permanently tagged, identified to species, diameters measured, and geographic coordinates mapped using a surveyor’s total station. All trees were measured in 1999, 2004, 2011, and 2017. In 2017, tree cores and additional tree measurements were collected from a stratified random sample based on the 2011 census data. Sampling strata included all six treatment combinations, the four dominant species (white fir, incense-cedar, sugar pine, and Jeffrey pine), three diameter classes (10-25 cm, 25-55 cm, and > 55 cm), and two local competition classes (high versus low competition). Local competition was quantified by generating Thiessen polygons derived from live tree geographic coordinates in the 2011 census, with Thiessen polygon area (m^2) around each tree as the metric of competition. Low and high competition classes were based on Theissen polygon areas greater than or less than the median polygon area within a given treatment combination. Five replicate trees were selected randomly for each combination of six treatments, four species, three diameter classes, and two competition classes, resulting in 740 (720 plus a few extra) trees sampled across gradients of tree size and localized competition in each treatment combination. For each sampled tree we recorded the species, diameter, height, live crown ratio, and canopy class (dominant, co-dominant, intermediate, overtopped), and two increment cores were collected at breast height on the uphill and parallel to slope sides of the tree. Two cores were extracted from each tree with a standard 5.15 mm increment borer, except one out of every five replicate trees had its second core collected with a 12 mm diameter increment borer for a companion study of carbon stable isotopes. Cores were taped onto wooden mounting sticks until they dried, then glued and sanded with progressively finer grit sandpaper to visualize tree-ring boundaries. Ring-widths were measured to the nearest 0.001 mm using either a high resolution flatbed scanner with WinDENDRO software (Regent Instruments, Quebec, Canada) or a stereo zoom microscope and Velmex Unislide TA tree-ring measuring system (Velmex, Bloomfield, New York). Tree ring series were cross-dated to ensure correct calendar year assignment of ring-widths using the dplR package (Bunn et al. 2021). Series that could not be cross-dated due to rotten or fragmented cores were discarded, resulting in cross-dated series for 1,401 of the 1,480 cores collected (731 of 740 trees sampled).

Tree competition was calculated from the 2011 census data. Following other studies in Sierra Nevada mixed-conifer forests, we defined local competition as a 10 m radius around the target tree (Das et al., 2011, 2008; Steel et al., 2021). For the 10 m radius around each tree, the basal area of all trees, trees less than 25 cm diameter, and trees greater than 25 cm diameter was calculated. Mean annual solar radiation and topographic wetness index were also calculated within a 10 m radius around each tree.

Gridded metrics of potential solar radiation and topographic wetness were generated for TEF using a digital terrain model (dtm) derived from airborne discrete return light detection and ranging (lidar) data. Lidar data was collected in October 2010 by Watershed Sciences Inc. (Portland, OR USA) as part of a larger acquisition for the USDA Forest Service. Lidar was collected using dual Leica ALS50 Phase II sensors mounted on a Cessna Caravan 208B flown at 1,100 and 1,500 m above ground level. Lidar survey specifications included a pulse rate of 83 kHz, mirror scan rate of 54 Hz, field of view ± 14⁰ from nadir, and opposing flight line swath overlap of 50%. Total pulse and ground pulse densities across the entire acquisition area were 8.8 pts/m2 and 0.89 pts/m2, respectively. We clipped the contractor provided 1 m resolution dtm to the geographic extent of TEF, then used the clipped dtm to calculate potential solar radiation and topographic wetness. Potential solar radiation on a sloping surface was calculated using the Areal Solar Radiation Model in ArcMap 10.7.1 (ESRI, 2019), an insolation model that accounts for atmospheric conditions, elevation, surface orientation, and surrounding topography (Fu and Rich, 2002). Topographic wetness index was calculated using the physically-based basin contribution model (Beven and Kirkby, 1979) using the following Equation:
Topographic wetness index=ln(α)/(tanβ+c)
Where α is the upslope contributing basin area (Moore et al., 1991) calculated with the watershed function in ArcMap, β is the slope at that cell, and c is a small constant (c = 0.01) to avoid division by zero in cells with flat terrain. Lower values of topographic wetness indicate greater topographic wetness.

Beven, K.J., Kirkby, M.J., 1979. A physically based, variable contributing area model of basin hydrology. Hydrol. Sci. J. 24, 43–69. https://doi.org/https://doi.org/10.1080/02626667909491834

Bunn, A., Korpela, M., Biondi, F., Campelo, F., Mérian, P., Qeadan, F., Zang, C., Buras, A., Cecile, J., Mudelsee, M., 2021. Package ‘dplR.’

Das, A., Battles, J., Stephenson, N.L., van Mantgem, P.J., 2011. The contribution of competition to tree mortality in old-growth coniferous forests. For. Ecol. Manage. 261, 1203–1213. https://doi.org/10.1016/j.foreco.2010.12.035

Das, A., Battles, J., van Mantgem, P.J., Stephenson, N.L., 2008. Spatial elements of mortality risk in old‐growth forests. Ecology 89, 1744–1756. https://doi.org/https://doi.org/10.1890/07-0524.1

ESRI, 2019. ArcGIS Desktop.

Fu, P., Rich, P.M., 2002. A geometric solar radiation model with applications in agriculture and forestry. Comput. Electron. Agric. 37, 25–35. https://doi.org/https://doi.org/10.1016/S0168-1699(02)00115-1

Moore, I.D., Grayson, R.B., Ladson, A.R., 1991. Digital terrain modelling: A review of hydrological, geomorphological, and biological applications. Hydrol. Process. 5, 3–30. https://doi.org/https://doi.org/10.1002/hyp.3360050103

Steel, Z.L., Goodwin, M.J., Meyer, M.D., Fricker, G.A., Zald, H.S.J., Hurteau, M.D., North, M.P., 2021. Do forest fuel reduction treatments confer resistance to beetle infestation and drought mortality? Ecosphere 12, e03344. https://doi.org/https://doi.org/10.1002/ecs2.3344


DATA-SPECIFIC INFORMATION FOR: tk.rwl

1. Number of variables: 1401

2. Number of cases/rows: 458

3. Variable List: 
	1401 unique alphanumeric identifiers for each crossed annual ring width series 

4. Missing data codes: 
	NA

5. Abbreviations used: 
	NA; not applicable

6. Other relevant information: 
	Data are cross dated ring widths to the nearest 0.001 mm. The longest series contained 458 annual years going back to 1560. However, the majority of the trees have fewer than 200 rings, resulting in a large amount of NA values.


DATA-SPECIFIC INFORMATION FOR: tk.lookup

1. Number of variables: 13

2. Number of cases/rows: 731

3. Variable List: 
	PLOT: Treatment unit ID
	SPECIES: 4 letter species identification code
	TAG: Tree numeric ID 
	Tree: Tree alphanumeric ID combining PLOT and TAG
	DATE: Date of tree core collection
	DBH: Diameter in centimeters at breast height (1.37 m above ground) collected in 2017
	HT: Height of tree in meters
	CRNR: Live crown ratio as a percentage
	CANCLASS: Categorical canopy position
	ISOTOPE: Did tree have a 12 mm diameter isotope core collected (Y) or only standard 5 mm diameter cores collected (N)
	STRATA: Cored tree sampling strata
	COREID_A: Unique ID for the A core collected from trees
	COREID_B: Unique ID for the B core collected from trees
4. Missing data codes: 
	NA

5. Abbreviations used: 
	NA; not applicable

6. Other relevant information: 
	None


DATA-SPECIFIC INFORMATION FOR: tk.comp.topo

1. Number of variables: 17

2. Number of cases/rows: 17

3. Variable List: 
	plot: Treatment unit ID
	species 4 letter species identification code
	tree: Tree alphanumeric ID combining PLOT and TAG
	cored: Was mapped tree cored (Y/N)
	dbh: Diameter in centimeters at breast height (1.37 m above ground) in 2017
	dbh.2011: Diameter in centimeters at breast height (1.37 m above ground) in 2011
	height: Height of tree in meters
	crown.ratio: Live crown ratio as a percentage
	ba.2011: basal area of tree in square meters, in 2011
	n.trees: number of trees within 10 meters of a tree
	n.trees: number of trees less than 25 cm DBH within 10 meters of a tree
	n.trees.ge25: number of trees greater or equal to 25 cm DBH within 10 meters of a tree
	ba: basal area in square meters of all trees within 10 m of a given tree
	ba.lt25: basal area in square meters of all trees less than 25 cm DBH within 10 m of a given tree
	ba.ge25: basal area in square meters of all trees greater or equal to 25 cm DBH within 10 m of a given tree
	psr: potential solar radiation
	twi: topographic wetness index
4. Missing data codes: 
	NA
5. Abbreviations used: 
	NA; not applicable
6. Other relevant information: 
	None
