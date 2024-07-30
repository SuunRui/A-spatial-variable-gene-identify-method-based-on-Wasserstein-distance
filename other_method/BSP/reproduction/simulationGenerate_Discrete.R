library(dplyr)
InputExp <- openxlsx::read.xlsx("mmc6.xlsx", 
                                colNames = FALSE)[,-1] %>% unlist() %>% as.numeric()
OutputDir <- "../Simulation_3D_Test/"


sim_discrete <- function(InputExp, 
                         OutputDir,
                         NGenes = 500, 
                         NX = 15, 
                         NY = 15, 
                         numSlices = 5, # 切片数量
                         Rad = 1, # 用于控制局部和全局基因表达模式的半径参数
                         Breaks = 8, # 用于控制局部基因表达模式中的中心点密度
                         Quant = 0.80, 
                         Noise_lv = 2, # 噪声水平
                         SEED = 1)
  {
  set.seed(SEED)
  SimData <- data.frame() # 初始化一个dataframe存储模拟数据
  N <- NX * NY
  for(i in 1:numSlices)
    {
    # 随机生成N个坐标点
    Coords <- spatstat.random::rpoispp(N, 
                                       win = spatstat.geom::owin(c(0, 1), 
                                                                 c(0, 1)))
    # 将坐标信息放到存储模拟数据的dataframe里面
    SimData <- rbind(SimData, 
                     data.frame(x = Coords$x * (sqrt(N) - 1), 
                                y = Coords$y * (sqrt(N) - 1), 
                                z = i))
    }
  N_Obs <- nrow(SimData) # 细胞的数量
  
  set.seed(SEED)
  # 对function(i)执行Num_gene次
  Sim_SVGs <- sapply(1:NGenes, 
                     function(Gene)
                       {
    return(sample(InputExp, # 通过采样每次生成单个基因的表达量的模拟数据
                  N_Obs, 
                  replace = TRUE)) # 从输入的基因表达数据做N次可重复抽样
                       })
  
  # local patterns / discrete patterns
  spikedcells <- function(coordDF,
                          Rad,
                          Breaks)
    {
    # 生成局部表达中心的坐标
    Center_pts <- as.data.frame(expand.grid(seq(from = 3,
                                                to = NX,
                                                by = Breaks), # 从 3 到 NX，以 Breaks 为步长的序列
                                           seq(from = 3,
                                               to = NY,
                                               by = Breaks) # 从 3 到 NX，以 Breaks 为步长的序列
                                           ))
    Center_pts$Z <- (numSlices + 1) / 2
    # 生成局部表达中心数量相同的在 [-2, 2] 范围内均匀分布的随机数
    Center_pts <- Center_pts + runif(length(Center_pts),
                                     -2,
                                     2)

    spikedCoords <- sapply(1:nrow(Center_pts), # nrow(Center_pts)为中心点的数量
                           function(Center_pt)
                             {
      # get distance from every cell to the center of the coordinate set;
      # assuming hot spot is centered here
      distFromCenter <- sp::spDists(
        x = coordDF %>% dplyr::select(x, y, z) %>% as.matrix(),
        y = matrix(as.numeric(Center_pts[Center_pt, ]),
                   ncol = 3),
        longlat = FALSE
      ) %>% as.vector()

      # identify coordinates in hot spot
      spikedCoord <- coordDF[distFromCenter <= (Rad),] %>% row.names() %>% as.numeric()
      return(spikedCoord)
      })
    spikedCoords <- unique(unlist(spikedCoords))
    return(spikedCoords)
    }

  # ==================global patterns / continuous patterns=====================
  spikedcells2 <- function(coordDF,
                           Rad, Breaks)
    {
    NCenter <- length(seq(from = 3,
                          to = NX,
                          by = Breaks)) * length(seq(from = 3,
                                                     to = NY,
                                                     by = Breaks))

    Phi <- 2 * pi * runif(NCenter) # 球坐标系中的角度
    Theta <- 2 * pi * runif(NCenter)# 球坐标系中的角度

    # 通过在球坐标系上运动的累积和来创建的三维坐标点
    Center_pts <- cbind((NX + 1)/2 + cumsum(2 * cos(Phi) * sin(Theta)),
                        (NY + 1)/2 + cumsum(2 * cos(Phi) * cos(Theta)),
                        cumsum(2 * sin(Phi)) + numSlices/2)
    #plot(cumsum(2*sin(Theta)), cumsum(2*cos(Theta)),type = "b", xlim = c(-10,10), ylim =c(-10,10) )

    spikedCoords <- sapply(1:length(Phi),
                           function(Center_pt)
                             {
      # get distance from every cell to the center of the coordinate set;
      # assuming hot spot is centered here
      distFromCenter <- sp::spDists(
        x = coordDF %>% dplyr::select(x, y, z) %>% as.matrix(),
        y = matrix(Center_pts[Center_pt, ],
                   ncol = 3),
        longlat = FALSE
      ) %>% as.vector()

      # identify coordinates in hot spot
      spikedCoord <- coordDF[distFromCenter <= Rad,] %>% row.names() %>% as.numeric()
      })
    spikedCoords <- unique(unlist(spikedCoords))

    # return coordinates
    return(spikedCoords)
    }
  # ==================global patterns / continuous patterns=====================
  # 
  Sim_SVGs1 <- Sim_SVGs[,1:floor(NGenes/2)]
  Sim_SVGs2 <- Sim_SVGs[,(1:ceiling(NGenes/2) + floor(NGenes/2))]
  
  # ===============================local pattern================================
  SpikedCellIndex <- spikedcells(SimData, 
                                 Rad, 
                                 Breaks)
  SpikedCells <- length(SpikedCellIndex)
  # 选取右分位点
  SpikedValues_upper <- InputExp[InputExp>=quantile(InputExp, 
                                                    Quant)]
  SpikedCellValues_upper <- sapply(1:floor(NGenes/2), 
                                   function(Gene)
                                     {
    return(sample(SpikedValues_upper, 
                  SpikedCells, # 要抽取的样本数量
                  replace = TRUE))
    })
  SpikedCellValues <- SpikedCellValues_upper
  SpikedCellValues <- as.matrix(SpikedCellValues)
  Sim_SVGs1[SpikedCellIndex, ] <- SpikedCellValues
  # =================================local pattern==============================
  # ================================global pattern==============================
  SpikedCellIndex <- spikedcells2(SimData, 
                                  Rad, 
                                  Breaks)
  SpikedCells <- length(SpikedCellIndex)
  SpikedValues_upper <- InputExp[InputExp>=quantile(InputExp, 
                                                    Quant)]
  SpikedCellValues_upper <- sapply(1:ceiling(NGenes/2), 
                                   function(Gene)
                                     {
    return(sample(SpikedValues_upper, 
                  SpikedCells, 
                  replace = TRUE))
                                     })
  SpikedCellValues <- SpikedCellValues_upper
  SpikedCellValues <- as.matrix(SpikedCellValues)
  Sim_SVGs2[SpikedCellIndex, ] <- SpikedCellValues
  # ================================global pattern==============================
  
  Sim_SVGs <- cbind(Sim_SVGs1, 
                    Sim_SVGs2)
  colnames(Sim_SVGs) <- paste0("SVG_", 
                               1:NGenes)
  SimData <- cbind(SimData, 
                   Sim_SVGs)
  
  perm_func <- function(Input_Data, 
                        Rep = 9)
    {
    NumGenes <- ncol(Input_Data) - 3
    NumCells <- nrow(Input_Data)
    set.seed(SEED)
    Output_Data <- sapply(1:NumGenes, 
                          function(Gene_Index)
                            {
      return(as.data.frame(sapply(1:Rep, 
                                  function(j) 
                                    return(Input_Data[sample(1:NumCells, 
                                                             replace = FALSE), 
                                                      Gene_Index + 3])
      )))})
    Output_Data <- do.call("cbind", 
                           Output_Data)
    return(Output_Data)
  }
  # Null基因相当于nSVGs
  Sim_NullGenes <- perm_func(SimData)
  # 9倍于SVG数量的NULL基因
  colnames(Sim_NullGenes) <- paste0("NULL_", 
                                    1:(NGenes * 9))
  
  SimData <- cbind(SimData, 
                   Sim_NullGenes)
  
  SD_mean <- mean(apply(SimData[,4:ncol(SimData)], 
                        2, 
                        sd))
  Noise_Matrix<- matrix(round(1*rnorm((nrow(SimData) * (ncol(SimData)-3)), 
                                      mean = 0, sd = SD_mean)),
                        nrow = nrow(SimData), 
                        ncol = (ncol(SimData)-3))
  SimData_AllExp <- SimData[,4:ncol(SimData)]
  SimData_AllExp <- SimData_AllExp + Noise_Matrix
  SimData_AllExp[SimData_AllExp<0] <- 0
  SimData[,4:ncol(SimData)] <- SimData_AllExp
  
  
  
  write.csv(SimData, 
            paste0(OutputDir,
                   "Discrete_",
                   N,
                   "by",
                   numSlices,
                   "_width",
                   Rad,
                   "_qt",
                   Quant*100,
                   "_Noise",
                   Noise_lv,
                   "_pw",
                   SEED,
                   ".csv"), row.names = FALSE, quote = FALSE)
  }






