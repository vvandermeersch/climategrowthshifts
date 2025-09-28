# Data from: Contemporary fires are less frequent but more severe in dry conifer forests of the southwestern United States

[https://doi.org/10.5061/dryad.98sf7m0sn](https://doi.org/10.5061/dryad.98sf7m0sn)

## Description of the data and file structure

Field data were collected in Survey123, then adapted into .csv files for import and analysis in R. Sampling was performed in June - August 2021, and involved visits to 74 previously sampled tree-ring fire-scar sites. At each site center or found sampled tree, a 10-meter radius plot was installed, in which we recorded tree diameter at breast height (dbh), species, and status (live or dead). We measured diameter and assigned species for downed logs and recorded an overall count of trees both live and dead, standing and down. We also completed a modified (simplified) CBI protocol to assess severity, and counted seedlings by species across the plot. Sampling was focused on 6 geographic areas where networks of fire history sites had been established prior to wildfires occurring over the past ten years (2011-2020): the Jemez, Rincon, Santa Catalina, Pinaleño, and Chiricahua Mountains, and the Kaibab Plateau.

**Tree species codes explained**

| *ABCO* | *Abies concolor*        | White fir               |
| :----- | :---------------------- | :---------------------- |
| *JUDE* | *Juniperus deppeana*    | Alligator juniper       |
| *JUSC* | *Juniperus scopulorum*  | Rocky Mountain juniper  |
| *PIAZ* | *Pinus arizonica*       | Arizona pine            |
| *PIED* | *Pinus edulis*          | Piñon pine              |
| *PIEN* | *Picea engelmannii*     | Engelmann spruce        |
| *PIFL* | *Pinus flexilis*        | Limber pine             |
| *PIPO* | *Pinus ponderosa*       | Ponderosa pine          |
| *PISF* | *Pinus strobiformis*    | Southwestern white pine |
| *PPTR* | *Populus tremuloides*   | Quaking aspen           |
| *PSME* | *Pseudotsuga menziesii* | Douglas-fir             |
| *QUGA* | *Quercus gambellii*     | Gambel oak              |

### Files and variables

#### File: field\_data\_all\_trees.csv

**Description:** all live and dead trees, standing and down, were identified to species and assessed as being killed by fire (if dead).

##### Variables

* Site Name: identification code for field plot
* ID: identification code for original fire history site
* species: tree species code (spelled out in metadata)
* DBH: diameter at breast height
* Status: live or dead at time of field sampling
* Killed by fire: qualitative assessment by field data collectors as to whether tree was killed by most recent fire
* NA values derived from original metadata (signify missing data)

#### File: field\_data\_seed.csv

**Description:** all seedlings were counted and identified to species in a 5-meter subplot within the larger plot area, centered at plot center

##### Variables

* Site Name: identification code for field plot
* ID: identification code for original fire history site
* mtnrange: geographic area of interest (one of six identified for field sampling)
* PIPO: plot seedling count for ponderosa pine
* PIST: plot seedling count for southwestern white pine
* PIAR5: plot seedling count for Arizona pine
* PSME: plot seedling count for Douglas fir
* ABCO: plot seedling count for white fir
* PIEN: plot seedling count for Engelmann spruce
* JUOS: plot seedling count for juniper
* POTR: plot seedling count for quaking aspen
* count_conifer: plot seedling count summed for all conifer species
* count_total: plot seedling count summed for all species

#### File: field\_data\_sev.csv

**Description:** a modified CBI (Composite Burn Index, typically conducted one-year post-fire) plot protocol was utilized to assess severity, with particular focus on trees and heavy fuels. all measurements performed by ocular estimation.

##### Variables

* Site Name: identification code for field plot
* ID: identification code for original fire history site
* mtnrange: geographic area of interest (one of six identified for field sampling)
* Latitude: latitudinal coordinate, in decimal degrees
* Longitude: longitudinal coordinate, in decimal degrees
* Aspect: direction of slope of plot
* Slope: average slope of plot, by ocular estimate, in percent
* % torch: proportion of tree canopy consumed by fire, averaged across all trees in plot area
* % scorch: proportion of tree canopy damaged, but not consumed by fire, averaged across all trees in plot area
* % mortality: proportion of trees in plot killed by fire
* average char height (m): maximum height along tree bole where bark was damaged by fire, averaged across all trees in plot area
* % consumption of heavy fuels: presumed percent reduction in 1000 hour Time Lag Fuel Moisture fuels on the ground, averaged across plot area
* torch_cbi: CBI value (ranging from 0 to 3) corresponding to % torch value assessed above
* scorch_cbi: CBI value (ranging from 0 to 3) corresponding to % scorch value assessed above
* mort_cbi: CBI value (ranging from 0 to 3) corresponding to % mortality value assessed above
* char_cbi: CBI value (ranging from 0 to 3) corresponding to char height value assessed above
* fuel_cbi: CBI value (ranging from 0 to 3) corresponding to % consumption of heavy fuels value assessed above
* cbi_field: total plot CBI value calculated by averaging five previous individual CBI values
* NA values signify attributes we couldn't assess in the field (eg. no scorch % if no live trees present pre-fire)

#### File: field\_data\_site.csv

**Description:** general site characteristics of each plot installed

##### Variables

* Site Name: identification code for field plot
* Date Visited: date of data collection
* Collected By: field sampling team
* Fire History: most recent fire which burned plot location, name and year
* Latitude: latitudinal coordinate, in decimal degrees
* Longitude: longitudinal coordinate, in decimal degrees
* x: longitudinal coordinate, in decimal degrees
* y: latitudinal coordinate, in decimal degrees
* Aspect: direction of slope of plot
* Slope: average slope of plot, by ocular estimate, in percent
* Vegetation Type: qualitative assessment of dominant vegetation in plot area
* Observations: notes on plot conditions and parameters, to add descriptive / qualitative data to dataset

#### File: field\_data\_treecount.csv

**Description:** summary of tree data for each plot

##### Variables

* Site Name: identification code for field plot
* Coordinates: latitude and longitude, in decimal degrees
* Aspect: direction of slope of plot
* Slope: average slope of plot, by ocular estimate, in percent
* Number Live Overstory: count of live trees in plot area, with DBH >= 12.7 cm
* Number Live Sapling: count of live trees in plot area, with DBH < 12.7 cm
* Total Live Tree: count of all live trees in plot area
* Number Standing Dead: count of standing dead trees in plot area, with DBH >= 12.7 cm 
* Number Dead & Down: count of logs