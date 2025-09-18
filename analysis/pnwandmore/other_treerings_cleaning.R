rm(list = ls())
wd <- '/home/victor/projects/climategrowthshifts/analysis/pnwandmore'

# ----------------------------------------------------------- #
# Clean and merge treering data from other sources than ITRDB #
# ----------------------------------------------------------- #

# Just an example below
folder <- 'data/other_treerings/direct_searching/Utah/datasearch/rawTreeRing/1'
rawdata <- read.csv(file.path(wd, folder, 'UT_rw.csv'))

# The metadata is quite clear in this example:
# UT_rw.csv:          Comma separated file containing ring width (RW) in millimeters for all available species (SPCD) and Years (Year)
# Ring widths (RW) can be linked to the inventory data with plot control number (PLT_CN) and tree control number (TRE_CN)

processed_data <- data.frame(
  dataset = '10.5849/jof.15-097', # this is the doi of the paper
  site_id = as.character(rawdata$PLT_CN),
  tree_id = as.character(rawdata$TRE_CN),
  year = rawdata$Year,
  rw_mm = rawdata$RW
)

# But we also want to now the coordinates of the sites
# It seems these data come from FIA (/!\ beware of duplicates with the FIA folder then! These data are maybe already in your folder!)
# So you can likely find the coordinates of the sites thanks to the plot control number