# README for EvansDeyHeilman_PinusEdulis_ringwidth_2024

# EvansDeyHeilman_PinusEdulis_ringwidth_2024
readme.txt (this file)


# EvansDeyHeilman_PinusEdulis_ringwidth_2024: 
This dataset and metadata was prepared for publication with the manuscript titled for PNAS in 2024 by coauthors: Margaret E. K. Evans, Sharmila M. N. Dey, Kelly A. Heilman, John R. Tipton, R. Justin DeRose, Stefan Klesse, Emily L. Schultz, and John D. Shaw

This directory contains 2 files with tree-ring data collected as part of the US Forest Service Forest Inventory and Analysis (FIA) program in the Interior Western US. The samples (increment cores) were collected from Pinus edulis trees associated with FIA plots in two different inventory designs: periodic inventories before the year 1999, and after 1999, annual inventories, in which a temporally and spatially stratified ~10% of plots are sampled each year, such that each plot is visited every ~10 years (DeRose et al., 2017; USDA Forest Service, 2015, 1999). Data stored in this repository includes increment cores from 1558 trees, most of which were sampled during periodic inventories, from tally trees, site trees, or representative trees adjacent to FIA plots. Of these trees with increment cores, many were in plots carried over to the annual design, but others were orphaned after the implementation of the annual design in 1999. 

Tree-ring widths were measured by researchers at the University of Arizona (in the Laboratory of Tree Ring Research) and at Utah State University.

We note that publicly available FIA data can be found here: https://www.fia.fs.fed.us/tools-data/

  1. PIED_Meta.txt
     This file contains the metadata associated with each tree from which an increment core was collected. The columns on this file are: "Core_CN","PLT_CN","STATECD","COUNTYCD","PLOT","SUBP","DIA","SPCD". The Core Control number (Core_CN) can be linked to the CN column in PIED_RW_Publish.txt. This file includes other variables and metadata as defined in the FIA users manual (see https://www.fia.fs.fed.us/library/database-documentation/ for more information). Additional information that connects to the FIA database include the  Plot control number (PLT_CN), State FIPS code (STATE_CD), county (COUNTYCD), plot number (PLOT), subplot number (SUBP), diameter in inches at the time of sampling (DIA), and species code (SPCD), which can be used to link these data to the publicly available dataset. In this dataset, longer Core_CN numbers indicate persistent unique identifiers for the tree cores within the FIA database, while shorter Core_CN numbers are a concatenation of COUNTYCD, PLOT, and SUBP. To match cores up to publicly available FIADB plot information, PLT_CN, STATECD, COUNTYCD, PLOT, and SUBP columns can be used.

  2. PIED_RW.txt
     This file includes the columns "Core_CN","Year", and "Growth". The radial tree-ring widths are in millimeters (Growth), the year in which that tree ring growth occurred (Year), and the core control number (Core_CN), which can be linked to the CORE_CN column to join the information in PIED_Meta_Publish.txt and PIED_RW_Publish.txt. Note that for both of these files, in order to preserve the data in the CN columns, they should be treated as text in most software programs, not as numbers. We note that outermost radial tree ring widths may not match the measure year (MEASYEAR) or inventory year (INVYR) contained in the publicly available FIA database if increment core samples were collected before/during the growing season, and thus the last full year of growth would occur before the measure year/inventory year. Additional reasons for mismatches between MEASYEAR and the last measured growth ring include loss of the bark and outmost years of growth during the transport, remounting, and sanding of cores. Year assignments were made based on alignment of marker years and cross dating.

For questions regarding the tree ring and forest inventory datasets contact:
  John Shaw (USFS) john.d.shaw@usda.gov

For questions about the associated manuscript prepared for PNAS contact:
      Margaret E. K. Evans margaretekevans@gmail.com
