set name "apo"

#set text_in [gets stdin]

source ${name}_median_coloring.tcl
render TachyonLOptiXInternal ${name}_Median_back.ppm

set name "ATP"

source ${name}_median_coloring.tcl
render TachyonLOptiXInternal ${name}_Median_back.ppm

set name "dATP"

source ${name}_median_coloring.tcl
render TachyonLOptiXInternal ${name}_Median_back.ppm


