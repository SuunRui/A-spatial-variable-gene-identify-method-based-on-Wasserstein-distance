#############################################################
# Generate simulation data
#############################################################

InputData <- openxlsx::read.xlsx("mmc6.xlsx",
                                 colNames =  FALSE)
#将基因表达数据提取出来
InputData_exp <- as.numeric(unlist(InputData[,-1]))

Num_genes <- c(2000)
Num_cells <- c(1000)
for(Num_cell in Num_cells){#这里的循环实际上只执行了一次，相当于从变量Num_cells提取元素出来进行循环
  for(Num_gene in Num_genes){
    SimData <- data.frame()#初始化模拟数据的dataframe
    #随机Num_cell对2D坐标
    Coords <- spatstat.random::rpoispp(Num_cell, 
                                       win = spatstat.geom::owin(c(0, 1), c(0, 1)))
    N <- Coords$n
    #将坐标信息放到存储模拟数据的dataframe里面
    SimData <- rbind(SimData, 
                     data.frame(x = Coords$x * (sqrt(N) - 1), 
                                y = Coords$y * (sqrt(N) - 1)))
    #对function(i)执行Num_gene次
    Sim_exp_mat <- sapply(1:Num_gene, function(i){
      #每次生成单个基因的表达量的模拟数据
      Sim_exp <- sample(as.numeric(InputData_exp), 
                        N, #细胞数量
                        replace=TRUE)#从输入的基因表达数据做N次可重复抽样
      return(Sim_exp)
    })
    #对基因命名
    colnames(Sim_exp_mat) <- paste0("Gene_", 1:Num_gene)
    #拼接dataframe和基因表达数据
    SimData <- cbind(SimData, Sim_exp_mat)
    
    
    write.csv(SimData,row.names = FALSE,
              paste0("Gene",
                     Num_gene,
                     "_Cell",
                     as.character(as.numeric(Num_cell)),
                     ".csv"))
  }
}