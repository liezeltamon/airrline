# Derive diversity metrics

# Diversity metrics to assess repertoire richness (total number of species/receptor sequence) 
# and evenness (how uniform the distribution of members for each species/receptor sequence)
# Info below from comparison of 12 metrics using TCR repertoire: 
# https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-025-02236-5
# Gini-Simpson, Pielou, and Basharin were the most robust in both simulated and experimental data

# git clone https://github.com/BorchLab/immApex.git
source("immApex/R/diversity.R")

# ----- FUNCTIONS -----

# x is a character vector of clone identifiers

# 1. Pielou, Basharin, d50, and Gini primarily describe evenness and highly correlate with one another
#.basharin <- related to shannon so read up on it
.pielou_evenness  <- function(x, ...) { pielou_evenness(as.numeric(table(x))) }
.d50_dom          <- function(x, ...) { d50_dom(as.numeric(table(x))) }
.gini_coef        <- function(x, ...) { gini_coef(as.numeric(table(x))) }

# 2. Richness is best captured by S index, next Chao1 and ACE which also consider information on evenness
#.s_index <-
.ace_richness     <- function(x, ...) { ace_richness(as.numeric(table(x))) }
.chao1_richness   <- function(x, ...) { chao1_richness(as.numeric(table(x))) }

# 3. Shannon, Inv.Simspon, D3, D4, and Gini.Simpson measure richness and increasingly more information on evenness
.shannon_entropy  <- function(x, ...) { shannon_entropy(as.numeric(table(x))) }
.norm_entropy     <- function(x, ...) { norm_entropy(as.numeric(table(x))) }
.inv_simpson      <- function(x, ...) { inv_simpson(as.numeric(table(x))) }
# With increasing Hill number (Shannon: q = 1; Inv.Simpson and Gini.Simpson: q = 2; 
# D3:  q = 3 and D4:  q = 4 see Eq. 6), the threshold for distinguishing changes
# in only Richness increases, except for Gini.Simpson, whose response is visually
# more similar to Hill’s number  q = 4 than q = 2 
.hill_d3          <- function(x, ...) { hill_q(3)(as.numeric(table(x))) }
.hill_d4          <- function(x, ...) { hill_q(4)(as.numeric(table(x))) }
.gini_simpson     <- function(x, ...) { gini_simpson(as.numeric(table(x))) }
