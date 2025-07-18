# README for FIA_PIPO_AZ_TREE_RING_DEC2021

# FIA_PIPO_AZ_TREE_RING_DEC2021
readme.txt (this file)

# FIA_PIPO_AZ_TREE_RING_DEC2021: This directory contains 2 files with tree ring data collected from the US Forest Service Forest Inventory and Analysis (FIA) program in the state of Arizona. These tree ring data are from Pinus ponderosa trees and were collected from FIA plots were sampled under two different inventory designs: periodic inventories before the year 1999, and after 1999, annual inventories, in which a temporally and spatially stratified ~10% of plots are sampled each year, such that each plot is visited every ~10 years (DeRose et al., 2017; USDA Forest Service, 2015, 1999). Data stored in this repository included increment cores from 518 individual trees sampled during periodic inventories (in 1995 or 1996). Of these trees with increment cores, many were in plots carried over to the annual design, but others were “orphaned” after the implementation of the annual design in 1999  Hence, these 518 trees have 1 or 2 measurements of diameter at breast height (DBH) taken between 1995-2010.

Tree ring widths were measured by researchers at University of Arizona, in the Laboratory of Tree Ring Research.

We note that publicly available FIA data can be found here: https://www.fia.fs.fed.us/tools-data/

  1. PIPOCores518Meta.txt
     This file contains the metadata associated with the tree cores on the plot. The Core Control number (CORE_CN) can be linked to the CN column in PIPOCores518RingWidths.txt. This file includes other variables and metadata as defined in the FIA users manual (see https://www.fia.fs.fed.us/library/database-documentation/ for more information). In addition to standard FIA metadata, this data includes some calculated variables, including Stand Density Index (SDI), calculated by the summation method. DateFirst and DateEnd, which refer to the start and end years for the tree ring data for each tree. RINGCOUNT is the number of growth rings in the record, and METHOD indicates how the tree ring widths were measured ("scope" indicates measured under a microscope and measuring stand). T1_DIA is the tree diameter measured during the first inventory (in inches), and T2_DIA is the tree diameter measured during the second inventory. Additional Control Numbers include the Plot control number (PLT_CN), tree control number (TRE_CN), and Site control number (SIT_CN), which can be used to link these data to the publicly available dataset.

  2. PIPOCores518RingWidths.txt
     This file includes the tree ring widths in millimeters (RW), the year in which that tree ring growth occurred (Year), and the core control number (CN), which can be linked to the CORE_CN column to join the information in PIPOCores518RingWidths.txt and PIPOCores518Meta.txt. Note that for both of these files, in order to preserve the CN columns, these CN columns need to be treated as text in most software programs, not as numbers.




For questions regarding the Arizona PIPO tree ring and forest inventory datasets contact:
  John Shaw (USFS)


For questions about the associated manuscript prepared for Global Change Biology contact:
      Kelly Heilman  kellyannheilman@gmail.com

