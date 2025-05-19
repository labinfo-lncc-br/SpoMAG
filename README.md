# Probability of Sporulation potential in MAGs (SpoMAG) 


## Scope

<p align="center">
<img src="SpoMAGlogo.png" alt="SpoMAG_logo" width="600"/>
</p>

SpoMAG is an R-based machine learning tool designed to predict the sporulation potential of Metagenome-Assembled Genomes (MAGs) from uncultivated Firmicutes species that could have the potential to undergo sporulation, such as those from the classes Bacilli and Clostridia. It leverages functional annotation and ensemble learning (Random Forest + Support Vector Machine) to identify likely spore-formers from genome annotation data.


The repository for SpoMAG is at GitHub on the https://github.com/labinfo-lncc/SpoMAG. In this website you can report a bug and get help.



## Citation

Paper under publication.



## Installation of the SpoMAG package

You can install the **SpoMAG** package directly from GitHub using:

```R
# Install devtools if not already installed
install.packages("devtools")

# Install SpoMAG from GitHub
devtools::install_github("labinfo-lncc-br/SpoMAG")


### Dependencies

SpoMAG depends on the following packages:

- dplyr, version XXX
- tidyr, version XXX
- tibble, version XXX
- readr, version XXX
- caret, version XXX
- randomForest, version XXX


## Quick start
This is a quick example using the included files sporulation.csv (a known spore-former) and asporogenic.csv (a known non-spore-formers).

# Load package
library(SpoMAG)

# Load example annotation tables
file_spor <- system.file("extdata", "sporulation.csv", package = "SpoMAG")
file_aspo <- system.file("extdata", "asporogenic.csv", package = "SpoMAG")

# Read files
df_spor <- readr::read_csv(file_spor)
df_aspo <- readr::read_csv(file_aspo)

# Step 1: Extract sporulation-related genes
genes_spor <- sporulation_gene_name(df_spor)
genes_aspo <- sporulation_gene_name(df_aspo)

# Step 2: Convert to binary matrix
bin_spor <- build_binary_matrix(genes_spor)
bin_aspo <- build_binary_matrix(genes_aspo)

# Step 3: Predict using ensemble model (preloaded in package)
model_path <- system.file("extdata", "modelos_stacking2.RData", package = "SpoMAG")

result_spor <- predict_sporulation(bin_spor, model_path)
result_aspo <- predict_sporulation(bin_aspo, model_path)

# View results
print(result_spor)
print(result_aspo)

Expected output includes columns such as:
- RF_Prob: probability from Random Forest
- SVM_Prob: probability from SVM
- Meta_Prediction: final ensemble prediction
- Meta_Prob_Esporulante: predicted probability of being a spore-former
