# Spatial Transcriptomics Analysis Pipeline

This pipeline consists of a set of Python scripts designed to analyze spatial transcriptomics data. The purpose of these scripts is to preprocess the data, fit gene expression functions to spatial coordinates, calculate Wasserstein distances, and identify variable genes based on these distances.

## Prerequisites

Ensure you have the following libraries installed:
- numpy
- pandas
- scipy
- scanpy
- anndata
- fitter

Install the required packages using pip:

```sh
pip install numpy pandas scipy scanpy anndata fitter
```

## Scripts Overview

### 1. `data_Preprocessing.py`

This script handles the reading and preprocessing of spatial transcriptomics data.

#### Functions:
- **`read_data(file_way)`**: Reads the data from a CSV file and converts it to an `anndata` format.
  - **Parameters**:
    - `file_way` (str): Path to the CSV file.
  - **Returns**: An `anndata` object containing spatial coordinates and gene expression data.

- **`data_pre(adata)`**: Preprocesses the `anndata` object by filtering genes, annotating mitochondrial genes, and normalizing the data.
  - **Parameters**:
    - `adata` (anndata): The annotated data object.
  - **Returns**: A pandas DataFrame with preprocessed gene expression data and spatial coordinates.

### 2. `caculate_wass_dist.py`

This script calculates the Wasserstein distance for gene expressions fitted to spatial coordinates.

#### Functions:
- **`sc_fit(x, expr)`**: Fits a gene expression function to the spatial coordinates using linear interpolation.
  - **Parameters**:
    - `x` (array): Coordinates of spots (2D or 3D).
    - `expr` (array): Gene expression values.
  - **Returns**: Interpolated function.

- **`sc_Wasserstein_Dist(x, expr, func)`**: Computes the Wasserstein distance between the fitted gene expression function and a mean function.
  - **Parameters**:
    - `x` (array): Coordinates of spots (2D or 3D).
    - `expr` (array): Gene expression values.
    - `func` (function): Interpolated function of gene expression.
  - **Returns**: Wasserstein distance.

- **`run(data, dataset_name)`**: Executes the calculation of Wasserstein distances for all genes in the dataset.
  - **Parameters**:
    - `data` (DataFrame): DataFrame containing spatial coordinates and gene expression data.
    - `dataset_name` (str): Name of the dataset.
  - **Returns**: DataFrame with genes and their calculated Wasserstein distances.

### 3. `fitter_dist.py`

This script fits various statistical distributions to the Wasserstein distances and determines the threshold for variable genes.

#### Functions:
- **`calculate_threshold_gene_site(store_sorted)`**: Fits distributions to the sorted Wasserstein distances and calculates the threshold for identifying variable genes.
  - **Parameters**:
    - `store_sorted` (array): Array with sorted Wasserstein distances.
  - **Returns**: Tuple containing the threshold gene site index, the number of variable genes, and the subset of the sorted array with variable genes.

### 4. `main.py`

This script orchestrates the entire process, combining the functionalities of the other scripts into a complete workflow.

### 5. `full_spot.py`

This script provides utility functions used in the preprocessing and analysis steps, such as calculating spatial coordinates and modifying gene expressions.

#### Functions:
- **`spatial_coordinates_knn()`**: Calculates K-nearest neighbors for spatial coordinates.
- **`spot_neighbor()`**: Identifies neighboring spots based on spatial coordinates.
- **`Modifying_gene_low_expression()`**: Modifies gene expressions with low values.
- **`Modifying_gene_high_expression()`**: Modifies gene expressions with high values.

## Usage

1. **Preprocess the Data**:
   ```python
   from data_Preprocessing import read_data, data_pre

   file_way = 'path/to/your/data.csv'
   adata = read_data(file_way)
   preprocessed_data = data_pre(adata)
   ```

2. **Calculate Wasserstein Distances**:
   ```python
   from caculate_wass_dist import run

   dataset_name = 'your_dataset_name'
   wass_distances = run(preprocessed_data, dataset_name)
   ```

3. **Fit Distributions and Identify Variable Genes**:
   ```python
   from fitter_dist import calculate_threshold_gene_site

   sorted_wass_dist = wass_distances['wass dist'].values
   threshold_gene_site, SVG_numbers, SVGs = calculate_threshold_gene_site(sorted_wass_dist)
   ```

## Example

An example usage of the complete workflow can be found in the `main.py` script, which integrates all the steps into a single script.

