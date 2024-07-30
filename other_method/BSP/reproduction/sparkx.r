##########################################################################################################################
# Analysis using SPARKX in the comparison
#
# Ref: https://xzhoulab.github.io/SPARK/
# https://xzhoulab.github.io/SPARK/02_SPARK_Example/
#
# Usage:
# Rscript sparkx.r /home/wangjue/svgGranularity/By_Cell_Numbers/Sim_data_Cells_500.csv
##########################################################################################################################

library('SPARK')
library(pryr)
gc()
args = commandArgs(trailingOnly=TRUE)
da<-read.csv("D:/A_study/A_study/ST_cluster/code_and_data/method_BSP/data/3DRA2/RA2.csv")
# passing from input
# da<-read.csv(args[1])
expr <- t(da[,4:dim(da)[2]])
info<-data.frame(x=da[,"x"],
                 y=da[,"y"], 
                 z=da[,'z'])
rownames(info)<-rownames(da)
colnames(expr)<-rownames(da)
location        <- as.matrix(info)

## debug
start_time <- Sys.time()
sparkX <- sparkx(expr,location,numCores=1,option="mixture")
write.csv(sparkX$res_mtest, 'D:/A_study/A_study/ST_cluster/code_and_data/result/SPARK-X/RA2.txt')
## Results
# head(sparkX$res_mtest)
## debug
end_time <- Sys.time()
print(mem_used())
print(end_time - start_time)

