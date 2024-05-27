#######################################
# Simulations: build 3D simulated datasets
# Yin Yihao, updated Tanuary 2024
#######################################

# this script builds simulated datasets based on empirical parameters from the
# Visium human DLPFC dataset


library(SpatialExperiment)
library(scran)
library(scater)
library(here)

# ---------
# load data
# ---------

spe <-read.csv("DLPFC.xlsx")

spe_full <- spe
dim(spe_full)

# --------------------
# empirical parameters
# --------------------
# MOBP gene is the gene expressed in white matter, data from fig.1 of the article
# calculate empirical simulation parameters based on expression of gene MOBP in
# white matter in human DLPFC dataset
# white matter is a special layer of DLPFC data set, see Appendix S1 of the nnSVG article

# identify gene and white matter
ix_gene <- which(names(spe) == "MOBP") # Index variable
is_gene <- names(spe) == "MOBP" # Logical vector
is_wm <- spe['WM'] == "1"
is_wm[is.na(is_wm)] <- FALSE# Replace all missing values with FALSE

# Real data set of gene "MOBP" expression in white matter and non-white matter regions
# ------------------------------------------------------------------------------
# This parameter can be changed to sample multiple genes or mean of multiple genes
# ------------------------------------------------------------------------------
# parameters: expressed vs. not expressed

# The logarithmic values of the MOBP gene in the white matter region in the real data set were selected
# This value will be used as the average expression level of genes in the expressed region in the simulation
par_meanLogcountsExpressed <- mean(spe[is_wm, is_gene])
par_varLogcountsExpressed <- var(spe[is_wm, is_gene])

# The logarithmic values of the MOBP gene in the non-white matter region of the real data set were selected
# This value will be used to model the average expression level of genes in non-expressed regions
par_meanLogcountsNotExpressed <- mean(spe[!is_wm, is_gene])
par_varLogcountsNotExpressed <- var(spe[!is_wm, is_gene])

# parameters: sparsity in expressed vs. not expressed regions
# Calculate the proportion of the logarithmic count of the gene "MOBP" to zero in the white matter region (sparsity)
par_sparsityExpressedRegion <-
  mean(spe[is_wm, is_gene] == 0)
# Calculate the proportion of the logarithmic count of the gene "MOBP" to zero in the non-white matter region (sparsity)
par_sparsityNotExpressedRegion <-
  mean(spe[!is_wm, is_gene] == 0)


# ---------------------
# additional parameters
# ---------------------

# additional simulation parameters

# number of non-expressed and expressed genes
n_genes <- 1000
n_genes_nonExpressed <- 900
n_genes_expressed <- 100

# strength of expression (relative to MOBP in WM vs. non-WM region)

lowExpression <-
  1 / 3 * (par_meanLogcountsExpressed - par_meanLogcountsNotExpressed) + par_meanLogcountsNotExpressed
highExpression <- 1 * par_meanLogcountsExpressed


# -------------------
# spatial coordinates
# -------------------


n_spots <- 10000


x <- runif(n_spots, min = 10, max = 100)
y <- runif(n_spots, min = 10, max = 100)
z <- runif(n_spots, min = 10, max = 100)


coords_3D <- data.frame(x = x, y = y, z = z)


# scale to range 0 to 1 in each dimension
coords_3D[, 1] <-
  (coords_3D[, 1] - min(coords_3D[, 1])) / (max(coords_3D[, 1]) - min(coords_3D[, 1]))
coords_3D[, 2] <-
  (coords_3D[, 2] - min(coords_3D[, 2])) / (max(coords_3D[, 2]) - min(coords_3D[, 2]))
coords_3D[, 3] <-
  (coords_3D[, 3] - min(coords_3D[, 3])) / (max(coords_3D[, 3]) - min(coords_3D[, 3]))

rownames(coords_3D) <- NULL

head(coords_3D)
dim(coords_3D)
summary(coords_3D)

# number of spots
n_spots <- nrow(coords_3D)


# ----------------------------------
# create masks for expressed regions
# ----------------------------------

centersLargeBandwidth_x <- 0.5
centersLargeBandwidth_y <- 0.5
centersLargeBandwidth_z <- 0.5
Big_CircleRadius <- 0.3
Small_CircelRadius <- c(0.2, 0.25, 0.3)
Big_InnerCircleRadius <- 0.15
Big_InnerCircleCenter_x <- c(0.35, 0.5, 0.65)
Big_InnerCircleCenter_y <- c(0.35, 0.5, 0.65)
Big_InnerCircleCenter_z <- c(0.35, 0.5, 0.65)

########## 大圆  ################
Largecircle_Mask <- matrix(FALSE, nrow =n_spots, ncol = 1)
for (i in seq_len(n_spots)) {
  x <- coords_3D[i, "x"]
  y <- coords_3D[i, "y"]
  z <- coords_3D[i, "z"]
    if (((x - centersLargeBandwidth_x) ^ 2 + (y - centersLargeBandwidth_y) ^ 2
          + (z - centersLargeBandwidth_z) ^2) < Big_CircleRadius ^ 2) { 
      Largecircle_Mask[i] <- TRUE 
  }
}

############# ring ########
# Generate all possible coordinate combinations 3*3*3
coordinates <- expand.grid(x = Big_InnerCircleCenter_x, 
                           y = Big_InnerCircleCenter_y, 
                           z = Big_InnerCircleCenter_z,
                           r = Big_InnerCircleRadius)


# x0,y0,z0,r0: great circle parameter
# biss: Adjust the coordinates of the smaller circle according to the larger circle
# r1: Small circle radius
mask1 <- function(x0, y0, z0, r0,biss,r1 ) {
  start_time <- Sys.time()
  Mask <- matrix(FALSE, nrow =n_spots, ncol =nrow(coordinates) )
  InnerCircleCenter_x <- c(x0-biss, x0, x0+biss)
  InnerCircleCenter_y <- c(y0-biss, y0, y0+biss)
  InnerCircleCenter_z <- c(z0-biss, z0, z0+biss)
  coordinates <- expand.grid(x = InnerCircleCenter_x, 
                             y = InnerCircleCenter_y, 
                             z = InnerCircleCenter_z,
                             r = r1)
  for (j in seq_len(nrow(coordinates))) {
    for (i in seq_len(n_spots)) {
      x <- coords_3D[i, "x"]
      y <- coords_3D[i, "y"]
      z <- coords_3D[i, "z"]
      if ((((x - x0)^2 + (y - y0)^2 + (z - z0)^2) < r0^2) &
          (((x - coordinates[j,1]) ^ 2 + (y - coordinates[j,2]) ^ 2
            + (z - coordinates[j,3]) ^2) 
           > coordinates[j,4] ^ 2)) {
          Mask[i, j] <- TRUE # If the condition is met, it becomes True
      }
    }
  }
  end_time <- Sys.time()
  running_time <- end_time - start_time
  print(paste("time cost：", running_time))
  return(Mask)
}
# Big hollow ball
Largeannulusmask<-mask1(centersLargeBandwidth_x,centersLargeBandwidth_y,
           centersLargeBandwidth_z,Big_CircleRadius,0.15,0.15)
set.seed(123)
Largeannulus_Mask <- Largeannulusmask[, sample(1:27, 9)]


Small_ExterCircleRadius <- 0.25
Small_ExterCircleCenter_x <- c(0.25, 0.75)
Small_ExterCircleCenter_y <- c(0.25, 0.75)
Small_ExterCircleCenter_z <- c(0.25, 0.75)

coordinates_1 <- expand.grid(
  x = Small_ExterCircleCenter_x,
  y = Small_ExterCircleCenter_y,
  z = Small_ExterCircleCenter_z,
  r = Small_ExterCircleRadius
)

mask_all <- matrix(nrow = n_spots, ncol = 0)
for (i in seq_len(nrow(coordinates_1) )){
  mask_all <- cbind(mask_all, mask1(coordinates_1[i,1],coordinates_1[i,2],
                                    coordinates_1[i,3],coordinates_1[i,4],0.1,0.15))
}
# -------------------mask2 func----------------
mask2 <- function(coordinates_1, r) {
  start_time <- Sys.time()
  Mask <- matrix(FALSE, nrow = n_spots, ncol = 1)
  for (j in seq_len(nrow(coordinates_1))) {
    for (i in seq_len(n_spots)) {
      x <- coords_3D[i, "x"]
      y <- coords_3D[i, "y"]
      z <- coords_3D[i, "z"]
      if (((x - coordinates_1[j, 1])^2 + (y - coordinates_1[j, 2])^2 
           + (z - coordinates_1[j, 3])^2) < r^2) {
        Mask[i, 1] <- TRUE 
      }
    }
  }
  
  end_time <- Sys.time()
  running_time <- end_time - start_time
  print(paste("time cost：", running_time))
  return(Mask)
}
Smallcircle_Mask <- matrix(nrow = n_spots, ncol = 0)
for (i in seq_along(Small_CircelRadius)) {
  Smallcircle_Mask <- cbind(Smallcircle_Mask, mask2(
    coordinates_1, Small_CircelRadius[i]))
    print('1')
  
}
#-------------------mask2 func-------

set.seed(1)
samplis <- 87
Smallannulus_Mask <- matrix(FALSE, nrow = n_spots, ncol = 0)
# results <- list()
for (i in 1:samplis) {
  set.seed(i)
  selected_columns <- seq(1, 27 * 8, by = 27)
  selected_columns <- sapply(selected_columns, function(x) sample(x:(x + 26), 1))
  result <- mask_all[, selected_columns[1]]
  for (k in 2:length(selected_columns)) {
    result <- result | mask_all[, selected_columns[k]]
  }
  Smallannulus_Mask <- cbind(Smallannulus_Mask, result)
}

Mask <- cbind(cbind(cbind(Largecircle_Mask, Largeannulus_Mask), Smallcircle_Mask), Smallannulus_Mask)
# -----------------------------------------------------------------------------
# functions to create log-expression values for noise genes and expressed genes
# -----------------------------------------------------------------------------

# note: using global arguments for parameters within functions
# note: set random seed before running functions

# noise genes

fn_createNoiseGene <- function() {
  logexpr <- rnorm(n_spots,
                   mean = par_meanLogcountsNotExpressed,
                   sd = sqrt(par_varLogcountsNotExpressed))
  q <- quantile(logexpr, par_sparsityNotExpressedRegion)
  logexpr[logexpr < q] <- 0
  logexpr[logexpr < 0] <- 0
  # return values
  logexpr
}


# expressed genes
fn_createExpressedGene <- function(meanExpression) {
  logexpr <- rnorm(n_spots,
                   mean = meanExpression,
                   sd = sqrt(par_varLogcountsExpressed))
  q <- quantile(logexpr, par_sparsityExpressedRegion)
  logexpr[logexpr < q] <- 0
  logexpr[logexpr < 0] <- 0
  # return values
  logexpr
}


# -------------------------
# create simulated datasets
# -------------------------

# gene names

gene_names <- c(sprintf("nonsvgs%03d", seq_len(900)), sprintf("svgs%03d", seq_len(100)))
spot_names <- sprintf("spot%04d", seq_len(n_spots))



fn_buildSimulatedData <- function(coords, mask, expressionStrength) {
  rowdata <- DataFrame(
    gene_name = gene_names,
    expressed = c(
      rep(FALSE, n_genes_nonExpressed),
      rep(TRUE, n_genes_expressed)
    ),
    expression_strength = c(
      rep(0, n_genes_nonExpressed),
      rep(expressionStrength, n_genes_expressed)
    )
  )
  coldata <- DataFrame(mask = mask)
  coords <- cbind(coords, coldata)
  rownames(coords) <- spot_names
  logcounts <- matrix(NA, nrow = n_spots, ncol = n_genes)
  i = 1
  for (g in seq_len(n_genes)) {
    logcounts[, g] <- fn_createNoiseGene()
    if (rowdata$expressed[g]) {
      logcounts[, g][mask[, i]] <-
        (fn_createExpressedGene(expressionStrength))[mask[, i]]
      i = i + 1
    }
  }
  rownames(logcounts) <- spot_names
  colnames(logcounts) <- gene_names
  logcounts_dataframe <- as.data.frame(logcounts)
  df <- cbind(coords, logcounts_dataframe)
  return(df)
}



# build datasets for combinations of parameters

set.seed(1)
spe_sim_largeBandwidth_fullExpr <-
  fn_buildSimulatedData(coords_3D, Mask, highExpression)
write.csv(spe_sim_largeBandwidth_fullExpr,
          "3DhighExpr.csv", 
          row.names = TRUE, 
          col.names = TRUE)



spe_sim_largeBandwidth_fullExpr <-
  fn_buildSimulatedData(coords_3D, Mask, lowExpression)
write.csv(spe_sim_largeBandwidth_fullExpr,
          "3DlowExpr.csv", 
          row.names = TRUE, 
          col.names = TRUE)
# build datasets for combinations of parameters

