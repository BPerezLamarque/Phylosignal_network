
**Guidelines to measure phylogenetic signal in a bipartite interaction network:**
====




This document indicates how to measure phylogenetic signal in a bipartite ecological network using the R-package RPANDA (Morlon *et al.*, 2016), *i.e.* to measure whether closely related species interact with similar partners.


**Citation:** Benoît Perez-Lamarque, Odile Maliet, Marc-André Selosse, Florent Martos, and Hélène Morlon, (*in prep.*)
*Do closely related species interact with similar partners? Testing for phylogenetic signal in ecological networks*.


**Contact:** Benoît Perez-Lamarque, benoit.perez.lamarque@gmail.com


# Contents:
**[Installation](#installation)**\
**[Measuring phylogenetic signal](#measuring-phylogenetic-signal)**\
**[Optionnal steps](#optionnal-steps)**





# Installation:


The functions part of the R-package RPANDA and can be installed from GitHub using devtools:

```r
library(devtools)
install_github("hmorlon/PANDA",ref="Benoit_phylosignal", dependencies = TRUE)

```



# Measuring phylogenetic signal:


<p align="center">
    <img src="https://github.com/BPerezLamarque/Phylosignal_network/blob/master/example/figures.png" width="500">
</p>

<p align="center">
    <b>A toy example of an interaction network between orchids (in green) and mycorrhizal fungi (in brown) informed with the phylogenetic trees of each guild.</b>
</p>

<br> <br>

First, load the example dataset of a mycorrhizal network between orchids and mycorrhizal fungi from La Réunion island (Martos *et al.*, 2012) along with the reconstructed phylogenetic trees of the orchids and the fungal OTUs:


```r
library(RPANDA)

data(mycorrhizal_network)

network <- mycorrhizal_network[[1]] # interaction matrix with orchids in columns and fungi in rows
tree_orchids <- mycorrhizal_network[[2]] # phylogenetic tree (phylo object)
tree_fungi <- mycorrhizal_network[[3]] # phylogenetic tree (phylo object)

```

<br> <br>


##  Step 1: Testing for phylogenetic signal in the species interactions


This first step uses the function  `phylosignal_network` to compute the phylogenetic signal in species interactions (do closely related species interact with similar partners?) using a simple Mantel test or the phylogenetic signal in the degree of generalism (do closely related species interact with the same number of partners?). Mantel tests measuring phylogenetic signal in species interactions can be computed using quantified or binary networks, with the Jaccard or UniFrac ecological distances.

```r

# compute phylogenetic signals in species intercations
phylosignal_network(network, tree_A = tree_orchids, tree_B = tree_fungi, method = "GUniFrac", correlation = "Pearson", nperm=10000)

```

| Option | Description |
| --- | --- |
| `network` | a matrix representing the ecological interaction network with species from guild A in columns and species from guild B in rows. |
| `tree_A` | a phylogenetic tree of the guild A (the columns of the interaction network). |
| `tree_B` | a phylogenetic tree of the guild B (the rows of the interaction network). |
| `method` | indicates which method is used to compute the phylogenetic signal in species interactions: you can choose "Jaccard_weighted" for computing ecological distances using Jaccard dissimilarities (or "Jaccard_binary" to not take into account the abundances of the interactions), or "GUniFrac" to compute the generalized UniFrac distances (or "UniFrac_unweighted" to not take into account the interaction abundances). |
| `correlation` |indicates which correlation is used in the Mantel test, among the Pearson, Spearman, or Kendall correlations. |
| `nperm` | number of permutations to evaluate the significance of the Mantel test.  |



The output of  `phylosignal_network` is:

| nb_A | nb_B | mantel_cor_A | pvalue_high_A | pvalue_low_A | mantel_cor_B | pvalue_high_B | pvalue_low_B |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 70 | 93 | -0.03 | 0.66 | 0.34  | 0.01  | 0.22 | 0.78 |

 which corresponds the number of orchid species (**nb_A**), the number of fungal species (**nb_B**), the Mantel correlation  between phylogenetic distances and ecological distances for orchids (**mantel_cor_A**), its associated upper p-value (**pvalue_high_A**), its associated lower p-value (**pvalue_low_A**), and  the Mantel correlation between phylogenetic distances and ecological distances for fungi (**mantel_cor_B**), its associated upper p-value (**pvalue_high_B**), abd its associated lower p-value (**pvalue_low_B**),

Thus, here we do not detect any significant phylogenetic signal in species interactions between the orchids and their mycorrhizal fungi. 

<br> <br>

##  Step 2: Testing for degrees of generalism


This second step also uses the function  `phylosignal_network` to compute the phylogenetic signal in degrees of generalism (do closely related species interact with the same number of partners?) using a simple Mantel test:

```r

# compute the phylogenetic signal in degrees of generalism for orchids 
phylosignal_network(network, tree_A = tree_orchids, method = "degree", correlation = "Pearson", nperm=10000)


```

| Option | Description |
| --- | --- |
| `network` | a matrix representing the ecological interaction network with species from guild A in columns and species from guild B in rows. |
| `tree_A` | a phylogenetic tree of the guild A (the columns of the interaction network). |
| `correlation` |indicates which correlation is used in the Mantel test, among the Pearson, Spearman, or Kendall correlations. |
| `nperm` | number of permutations to evaluate the significance of the Mantel test.  |



The output of  `phylosignal_network` is then:

| nb_A | nb_B | mantel_cor_A | pvalue_high_A | pvalue_low_A | 
| --- | --- | --- | --- | --- |
| 70 | 93 | 0.02 | 0.34 | 0.66  |

which corresponds the number of orchid species (**nb_A**), the number of fungal species (**nb_B**), and the Mantel correlation between phylogenetic distances and degree difference distances for orchids (**mantel_cor_A**), its associated upper p-value (**pvalue_high_A**), its associated lower p-value (**pvalue_low_A**).

Thus, here we do not detect any significant phylogenetic signal in degrees of generalism for the orchids.


<br> <br>


# Optionnal steps:

## Option 1: investigate clade-specific phylogenetic signals (simple Mantel tests with Bonferroni correction)




This first option uses the function  `phylosignal_sub_network` to  compute the clade-specific phylogenetic signals in species interactions. For each node of the tree A having a certain number of descending species, it computes the phylogenetic signal in the resulting sub-network by performing a Mantel test between the phylogenetic distances and the ecological distances for the given sub-clade of tree A. Mantel tests can be computed using quantified or binary networks, with the Jaccard or UniFrac ecological distances. The results of the clade-specific phylogenetic signal analysis can be represented using the function  `plot_phylosignal_sub_network`.

```r
# compute clade-specific phylogenetic signals in species intercations for orchids

results_clade_A <- phylosignal_sub_network(network, tree_A = tree_orchids, tree_B = tree_fungi, method = "GUniFrac", correlation = "Pearson", nperm=100000, minimum=10)

plot_phylosignal_sub_network(tree_A = tree_orchids, results_clade_A)

```

| Option | Description |
| --- | --- |
| `network` | a matrix representing the ecological interaction network with species from guild A in columns and species from guild B in rows. |
| `tree_A` | a phylogenetic tree of the guild A (the columns of the interaction network). |
| `tree_B` | a phylogenetic tree of the guild B (the rows of the interaction network). |
| `method` | indicates which method is used to compute the phylogenetic signal in species interactions: you can choose "Jaccard_weighted" for computing ecological distances using Jaccard dissimilarities (or "Jaccard_binary" to not take into account the abundances of the interactions), or "GUniFrac" to compute the generalized UniFrac distances (or "UniFrac_unweighted" to not take into account the interaction abundances). |
| `correlation` |indicates which correlation is used in the Mantel test, among the Pearson, Spearman, or Kendall correlations. |
| `nperm` | number of permutations to evaluate the significance of the Mantel test.  |
| `minimum` | indicates the minimal number of descending species for a node in tree A to compute its clade-specific phylogenetic signal.  |



The output of  `phylosignal_sub_network` corresponds to a table where each line corresponds to a tested clade and which contains 8 columns: the name of the node (*node*), the number of species in the sub-clade A (*nb_A*), the number of species in guild B associated with the sub-clade A (*nb_B*), the Mantel correlation for guild A (*mantel_cor*), its associated upper p-value (*pvalue_high*), its associated lower p-value (*pvalue_low*), and the corresponding Bonferroni corrected p-values (*pvalue_high_corrected* and *pvalue_low_corrected*).

The representation of the results using `plot_phylosignal_sub_network` is a phylogenetic tree with nodes colored according to the clade-specific phylogenetic signals. Blue nodes are not significant (based in the Bonferonni correction), grey nodes are not tested (less than  `minimum` descending tips), and orange-red nodes represent significant phylogenetic signals and their color indicates the strenght of the correlation.


<br> <br>



## Option 2: test the robustness of the findings to phylogenetic uncertainty and/or sampling bias

This can be achieved by testing the robustness of the results to phylogenetic uncertainty (using a Bayesian posterior of tree, by investigating the influence of polytomy…) or sampling bias (by subsampling over-represented species/clades). See Perez-Lamarque *et al.* (*in prep.*) for more details.
