"
This script will convert a csv format file to a xlsx format file.

SYNTAX:
   Rscript Convert2XL.R inputCSVfileName outputXLSXfileName

Author: Firoj Mahmud and Soumyadeep Nandi
Date: November 18, 2020
"
library(openxlsx)
args = commandArgs(trailingOnly=TRUE)
data=read.csv(args[1])
write.xlsx(data, file= args[2], asTable = FALSE)
