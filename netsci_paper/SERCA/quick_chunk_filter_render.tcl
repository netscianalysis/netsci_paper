## Simple script to be used for rendering the different beta values 
## Based on the high or low threshold use of the column percentile. 
## Meant to look at chunks of protins that is highly correlated with a 
## lot of the rest of the protein. 

set name "apo"

#set text_in [gets stdin]

source ${name}_delta_25_coloring.tcl
render TachyonLOptiXInternal ${name}_filtered_front.ppm

set name "ATP"

source ${name}_delta_25_coloring.tcl
render TachyonLOptiXInternal ${name}_filtered_front.ppm

set name "dATP"

source ${name}_delta_25_coloring.tcl
render TachyonLOptiXInternal ${name}_filtered_front.ppm


