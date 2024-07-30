#############################################################
# Generate 3D simulations for 3D vs 2D meta-analysis comparision
#############################################################

InputDir <- "..\\Data\\3D_Standard_Sim\\"
OutputDir <- "..\\Data\\3D_Standard_Sim_2DSlice\\"

InputFiles <- list.files(InputDir)
for(InputFile in InputFiles){
  InputData <- read.csv(paste0(InputDir, InputFile))
  for(i in 1:10){
    InputData_Slice <- subset(InputData, z == i)#根据z轴的数值分成10个切片
    InputData_Slice <- InputData_Slice[,-3]
    #分别储存十个切片的基因表达数据
    write.csv(InputData_Slice, row.names = FALSE, 
              paste0(OutputDir, gsub('\\.csv', paste0('_Slice', i, '\\.csv'), InputFile)))
  }
}

