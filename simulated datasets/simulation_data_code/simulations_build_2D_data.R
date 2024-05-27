#######################################
# Simulations: build 3D simulated datasets
# Sun Rui, updated Tanuary 2024
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
spe <- read.csv('DLPFC.xlsx')
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
ix_gene <- which(names(spe) == "MOBP") 
is_gene <- names(spe) == "MOBP" 
is_wm <- spe['WM'] == "1"
is_wm[is.na(is_wm)] <- FALSE 


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
#Calculate the proportion of the logarithmic count of the gene "MOBP" to zero in the white matter region (sparsity)
par_sparsityExpressedRegion <-
  mean(spe[is_wm, is_gene] == 0)
# Calculate the proportion of the logarithmic count of the gene "MOBP" to zero in the non-white matter region (sparsity)
par_sparsityNotExpressedRegion <-
  mean(spe[!is_wm, is_gene] == 0)

# check
# par_meanLogcountsExpressed  # 2.646526
# par_varLogcountsExpressed  # 0.0.4128669
# par_meanLogcountsNotExpressed  # 0.5142754
# par_varLogcountsNotExpressed  # 0.5661618
# par_sparsityExpressedRegion  # 0.01949318
# par_sparsityNotExpressedRegion  # 0.6353167


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

# using spatial coordinates from Visium slide in human DLPFC dataset
n_spots <- 2000

# Generate a random number of 10-100
x <- runif(n_spots, min = 10, max = 100)
y <- runif(n_spots, min = 10, max = 100)


# Save the result to a data box
coords <- data.frame(x = x, y = y)

# scale to range 0 to 1 in each dimension
coords[, 1] <-
  (coords[, 1] - min(coords[, 1])) / (max(coords[, 1]) - min(coords[, 1]))
coords[, 2] <-
  (coords[, 2] - min(coords[, 2])) / (max(coords[, 2]) - min(coords[, 2]))

rownames(coords) <- NULL

head(coords)
dim(coords)
summary(coords)

# number of spots
n_spots <- nrow(coords)


# ----------------------------------
# create masks for expressed regions
# ----------------------------------

# create masks to identify spatial coordinates within expressed regions

# -----------------------------Radius-----------------------------
# radius of expressed regions
# ---------------big annulus---------------
Big_ExterCircleRadius <- 0.3
Big_InnerCircleRadius <- 0.15
# ---------------big circle---------------
Big_CircleRadius <- 0.3
# ---------------small annulus---------------
Small_ExterCircleRadius <- 0.15
Small_InnerCircleRadius <- 0.1
# ---------------small circle---------------
Small_CircleRadius <- c(0.2, 0.25, 0.30)
# ------------------------------Radius---------------------------

# ----------------------------centers----------------------------
# ---------------big annulus---------------
Big_ExterCircleCenter_x <- 0.5
Big_ExterCircleCenter_y <- 0.5
Big_InnerCircleCenter_x <- c(0.35, 0.5, 0.65)
Big_InnerCircleCenter_y <- 0.5
# ---------------big circle---------------
Big_CircleCenter_x <- 0.5
Big_CircleCenter_y <- 0.5
# ---------------small annulus---------------
Small_ExterCircleCenter_x <- c(0.25, 0.75)
Small_ExterCircleCenter_y <- c(0.25, 0.75)
Small_innerCircleCenter_x <- c(0.2, 0.25, 0.3, 0.7, 0.75, 0.8)
Small_innerCircleCenter_y <- c(0.25, 0.75)
# ---------------small circle---------------
Small_CircleCenter_x <- c(0.25, 0.75)
Small_CircleCenter_y <- c(0.25, 0.75)
# ----------------------------centers----------------------------


# ---------------big annulus---------------
# Initialize mask to FALSE
Big_Annulus_Mask1 <- matrix(FALSE, nrow = n_spots, ncol = 1)
for (i in seq_len(n_spots)) {
  x <- coords[i, "x"]
  y <- coords[i, "y"]
  if (((x - Big_ExterCircleCenter_x) ^ 2 + (y - Big_ExterCircleCenter_y) ^ 2)
      < Big_ExterCircleRadius ^ 2 &
      ((x - Big_InnerCircleCenter_x[1]) ^ 2 + (y - Big_InnerCircleCenter_y) ^ 2)
      > Big_InnerCircleRadius ^ 2) {
    Big_Annulus_Mask1[i] <- TRUE # If the condition is met, it becomes True
  }
}

# ---------
Big_Annulus_Mask2 <- matrix(FALSE, nrow = n_spots, ncol = 1)
for (i in seq_len(n_spots)) {
  x <- coords[i, "x"]
  y <- coords[i, "y"]
  if (((x - Big_ExterCircleCenter_x) ^ 2 + (y - Big_ExterCircleCenter_y) ^ 2)
      < Big_ExterCircleRadius ^ 2 &
      ((x - Big_InnerCircleCenter_x[2]) ^ 2 + (y - Big_InnerCircleCenter_y) ^ 2)
      > Big_InnerCircleRadius ^ 2) {
    Big_Annulus_Mask2[i] <- TRUE 
  }
}

# ---------
Big_Annulus_Mask3 <-matrix(FALSE, nrow = n_spots, ncol = 1)
for (i in seq_len(n_spots)) {
  x <- coords[i, "x"]
  y <- coords[i, "y"]
  if (((x - Big_ExterCircleCenter_y) ^ 2 + (y - Big_ExterCircleCenter_y) ^ 2)
      < Big_ExterCircleRadius ^ 2 &
      ((x - Big_InnerCircleCenter_x[3]) ^ 2 + (y - Big_InnerCircleCenter_y) ^ 2)
      > Big_InnerCircleRadius ^ 2) {
    Big_Annulus_Mask3[i] <- TRUE
  }
}



# ---------------big circle---------------
Big_Circle_Mask1 <- matrix(FALSE, nrow = n_spots, ncol = 1)

for (i in seq_len(n_spots)) {
  x <- coords[i, "x"]
  y <- coords[i, "y"]
  if (((x - Big_CircleCenter_x) ^ 2 + (y - Big_CircleCenter_y) ^ 2)
      < Big_CircleRadius ^ 2) {
    Big_Circle_Mask1[i] <- TRUE # 满足条件的都变为True
  }
}

# ---------------small annulus---------------
Small_Annulus_Mask <- matrix(FALSE, nrow = n_spots, ncol = 81)
# Initialize an empty matrix to store the horizontal coordinates of different inner circles, there are 81 cases
# The circle on the bottom left is 11, the circle on the bottom right is 21, 
#the circle on the top left is 21, and the circle on the top right is 22
Small_innerCircleCenter <-
  expand.grid(
    Small_innerCircleCenter_x[1:3],
    Small_innerCircleCenter_x[4:6],
    Small_innerCircleCenter_x[1:3],
    Small_innerCircleCenter_x[4:6]
  )
colnames(Small_innerCircleCenter) <- c('11', '21', '12', '22')
count = 1
for (i in seq_along(Small_innerCircleCenter[,1])) {
  for (j in seq_len(n_spots)) {
    for (k in seq_along(Small_ExterCircleCenter_y)) {
      for (l in seq_along(Small_ExterCircleCenter_x)) {
        x <- coords[j, "x"]
        y <- coords[j, "y"]
        if (((x - Small_ExterCircleCenter_x[l]) ^ 2 + 
             (y - Small_ExterCircleCenter_y[k]) ^2) 
            < Small_ExterCircleRadius ^ 2 &
            ((x - Small_innerCircleCenter[i, as.character(l * 10 + k)])^2 + 
             (y - Small_innerCircleCenter_y[k]) ^ 2) 
            > Small_InnerCircleRadius ^2) {
          Small_Annulus_Mask[j, i] <- TRUE
        }
      }
    }
  }
}

# ---------------small circle---------------
Small_Circle_Mask <- matrix(FALSE, nrow = n_spots, ncol = 3)
for (i in seq_len(3)) {
  for (j in seq_len(n_spots)) {
    for (k in seq_along(Small_CircleCenter_y)) {
      for (l in seq_along(Small_CircleCenter_x)) {
        x <- coords[j, "x"]
        y <- coords[j, "y"]
        if (((x - Small_CircleCenter_x[l]) ^ 2 +
             (y - Small_CircleCenter_y[k]) ^ 2) < Small_CircleRadius[i] ^
            2)
        {
          Small_Circle_Mask[j, i] <- TRUE
        }
      }
    }
  }
}
Mask = matrix(nrow = n_spots, ncol = 0)
Mask = cbind(cbind(cbind(Mask, Big_Annulus_Mask1), Big_Annulus_Mask2), Big_Annulus_Mask3)
Mask = cbind(cbind(cbind(Mask, Big_Circle_Mask1), Small_Annulus_Mask), Small_Circle_Mask)
Mask = cbind(Mask, Mask[,1:12])

# -----------------------------------------------------------------------------
# functions to create log-expression values for noise genes and expressed genes
# -----------------------------------------------------------------------------

# note: using global arguments for parameters within functions
# note: set random seed before running functions

# noise genes

fn_createNoiseGene <- function() {
  # Extract random values from a normal distribution
  logexpr <- rnorm(n_spots,
                   mean = par_meanLogcountsNotExpressed,
                   sd = sqrt(par_varLogcountsNotExpressed))
  # Truncate values to achieve sparsity
  q <- quantile(logexpr, par_sparsityNotExpressedRegion)
  logexpr[logexpr < q] <- 0
  # Remove any negative values
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

# function to build simulated dataset

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
  # coords <- cbind(coords, coldata)
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
  fn_buildSimulatedData(coords, Mask, highExpression)
write.csv(spe_sim_largeBandwidth_fullExpr,
          "2DhighExpr.csv", 
          row.names = TRUE, 
          col.names = TRUE)
set.seed(1)
spe_sim_largeBandwidth_fullExpr <-
  fn_buildSimulatedData(coords, Mask, lowExpression)
write.csv(spe_sim_largeBandwidth_fullExpr,
          "2DlowExpr.csv", 
          row.names = TRUE, 
          col.names = TRUE)

