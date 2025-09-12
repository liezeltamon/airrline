#### Notes receptor feature functions

```
# Types of features as biomarkers:
# (done) Diversity and expansion metrics (include ones like D5, D10 (hill numbers) that
# can differentiate polyclonal from monoclonal and not expanded) - numeric per sample, subsample then vegan, DescTools etc.

# (done) k-mer usage - per sequence and proportion, immunarch (done per sample)
# (done) AA motif usage - per sequence so proportion, GLIPH2 and alakazam::countGenes() to get proportions
# (done) Physiochemical properties - numeric per sequence so average, alakazam
# (check if it gives allele) + custom functions to extract family and alakazam::countGenes() to get proportions

# (done) Allele, gene or family usage - per sequence so proportion, dandelion 

# Single-cell specific:
# (done) Clonal overlap between cell types - custom function
# (done) Inter/intra-clonality of cell types - custom function

# Probability of generation - per sequence so average, OLGA (python),SHazaM
# (done) Convergence - numeric per sample, ratio of sequence nt / sequence aa

# (done) Somatic mutation frequency - from dandelion

# B-cell specific:
# Class switching - proportions of isotypes and clonal overlap using isotype status,
# separate by isotype but pool IgD and IgM (in naive, these are coexpressed))
# Mutated isotypes - Mutation of certain genes e.g. mutation of IgD/IgM etc (CLL: unmutated vs mutated IGHV status = prognostic marker.)

# Tools:
# immcantation (R)
# OLGA (python)
# immunarch (R) - table is per sequence

# Functions to assign to .calculate_measure() applying to vector of clone ids (x)

# Diversity metrics to assess repertoire richness (total number of species/receptor sequence) 
# and evenness (how uniform the distribution of members for each species/receptor sequence)
# Info below from comparison of 12 metrics using TCR repertoire: 
# https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-025-02236-5
# Gini-Simpson, Pielou, and Basharin were the most robust in both simulated and experimental data

# Probability of receptor generation
# Useful for detecting disease-specific clones that might be hard to generate but are expanded in response to disease
# OLGA implemented in python packages as detailed here: https://github.com/statbiophys/OLGA/issues/14
# a. Faster numba implementation: https://github.com/dweb0/OLGA
# b. Parallel implementation in soNNia (by OLGA developer): https://github.com/statbiophys/soNNia
# "The SONIA package has additional functions but if you are only interested in pgen evaluation you can look at the 'sonia-evaluate' command with the '--pgen' option."
# c. Also faster implementation in tcrdist3: https://tcrdist3.readthedocs.io/en/latest/pGen.html#probability-of-generation

```