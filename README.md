
**Guidelines to measure the phylogenetic signal in a bipartite interaction network:**
====




This tutorial explains how to measure the phylogenetic signal in species interactions in a bipartite ecological network using the R-package RPANDA (Morlon *et al.*, 2016), *i.e.* how to test whether closely related species interact with similar partners.


**Citation:** Benoît Perez-Lamarque, Odile Maliet, Marc-André Selosse, Florent Martos, and Hélène Morlon,
*Do closely related species interact with similar partners? Testing for phylogenetic signal in interaction networks*, bioRxiv, doi: https://doi.org/10.1101/2021.08.30.458192.


**Contact:** Benoît Perez-Lamarque, benoit.perez.lamarque@gmail.com


Simulations used to investigate the performances of the different approaches to measure the phylogenetic signal in species interactions in bipartite interaction networks are available in the folder "simulations" with an associated README. Functions to simulate the BipartiteEvol interaction networks are available in RPANDA, while the script used to simulate interaction networks with a phylogenetic signal in the number of partners (and not their identity) can be found [here](https://github.com/BPerezLamarque/Phylosignal_network/blob/master/example/sim_signal_number_partners.R).


# Contents:
**[Installation](#installation)**\
**[Measuring the phylogenetic signal](#measuring-phylogenetic-signal)**\
**[Optional steps](#optional-steps)**





# Installation:


The functions are part of the R-package RPANDA and can be installed from GitHub using devtools:

```r
library(devtools)
install_github("hmorlon/PANDA",ref="Benoit_phylosignal", dependencies = TRUE)

```



# Measuring the phylogenetic signal:


<p align="center">
    <img src="https://github.com/BPerezLamarque/Phylosignal_network/blob/master/example/figures.png" width="500">
</p>

<p align="center">
    <b>A toy example of an interaction network between orchids (in green) and mycorrhizal fungi (in brown) informed with the phylogenetic trees of each guild.</b>
</p>

First, you can load the example dataset of a mycorrhizal network between orchids and mycorrhizal fungi from La Réunion island (Martos *et al.*, 2012) along with the reconstructed phylogenetic trees of the orchids and the fungi:


```r
library(RPANDA)

data(mycorrhizal_network)

network <- mycorrhizal_network[[1]] # interaction matrix with orchids in columns and fungi in rows
tree_orchids <- mycorrhizal_network[[2]] # phylogenetic tree (phylo object)
tree_fungi <- mycorrhizal_network[[3]] # phylogenetic tree (phylo object)

```

<br> <br>


##  Step 1: Testing for the phylogenetic signal in the species interactions


This first step uses the function  `phylosignal_network` to compute the phylogenetic signal in species interactions (do closely related species interact with similar partners?) using a simple Mantel test. Mantel tests measuring the phylogenetic signal in species interactions can be computed using quantified or binary networks, with the Jaccard or UniFrac ecological distances.

```r

# compute phylogenetic signals in species interactions
phylosignal_network(network, tree_A = tree_orchids, tree_B = tree_fungi, method = "GUniFrac", correlation = "Pearson", nperm=10000)

```

| Option | Description |
| --- | --- |
| `network` | a matrix representing the bipartite interaction network with species from guild A in columns and species from guild B in rows. |
| `tree_A` | a phylogenetic tree of the guild A (the columns of the interaction network). |
| `tree_B` | a phylogenetic tree of the guild B (the rows of the interaction network). |
| `method` | indicates which method to use to compute the phylogenetic signal in species interactions: you can choose "Jaccard_weighted" for computing ecological distances using Jaccard dissimilarities (or "Jaccard_binary" to not take into account the abundances of the interactions), or "GUniFrac" to compute the weighted (or generalized) UniFrac distances (or "UniFrac_unweighted" to not take into account the interaction abundances). |
| `correlation` |indicates which correlation to use in the Mantel test, among the Pearson, Spearman, or Kendall correlations. |
| `nperm` | indicates the number of permutations to evaluate the significance of the Mantel test.  |

<br> <br>

The output of  `phylosignal_network` is:

| nb_A | nb_B | mantel_cor_A | pvalue_high_A | pvalue_low_A | mantel_cor_B | pvalue_high_B | pvalue_low_B |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 70 | 93 | -0.03 | 0.66 | 0.34  | 0.01  | 0.22 | 0.78 |

 which corresponds to the number of orchid species (**nb_A**), the number of fungal species (**nb_B**), the Mantel correlation  between the phylogenetic distances and ecological distances for orchids (**mantel_cor_A**), its associated upper p-value (**pvalue_high_A**, *i.e.* the fraction of permutations that led to higher correlation values), its associated lower p-value (**pvalue_low_A**, *i.e.* the fraction of permutations that led to lower correlation values), and  the Mantel correlation between the phylogenetic distances and ecological distances for fungi (**mantel_cor_B**), its associated upper p-value (**pvalue_high_B**), abd its associated lower p-value (**pvalue_low_B**),

Thus, here we do not detect any significant phylogenetic signal in species interactions between the orchids and their mycorrhizal fungi (all p-values>0.05). 

<br> <br>

##  Step 2: Testing for the phylogenetic signal in the number of partners


This second step also uses the function  `phylosignal_network` to compute the phylogenetic signal in the number of partners (the degree; *i.e.* do closely related species interact with the same number of partners?) using a simple Mantel test. The goal of this step is to verify if there is a significant phylogenetic signal in species interactions (step 1), whether this phylogenetic signal in species interactions can be caused by a phylogenetic signal in the number of partners (confounding effect), rather than the identity of the interacting partners.

```r

# compute the phylogenetic signal in the number of partners for orchids 
phylosignal_network(network, tree_A = tree_orchids, method = "degree", correlation = "Pearson", nperm=10000)


```

| Option | Description |
| --- | --- |
| `network` | a matrix representing the bipartite interaction network with species from guild A in columns and species from guild B in rows. |
| `tree_A` | a phylogenetic tree of the guild A (the columns of the interaction network). |
| `correlation` | indicates which correlation to use in the Mantel test, among the Pearson, Spearman, or Kendall correlations. |
| `nperm` | indicates the number of permutations to evaluate the significance of the Mantel test.  |

<br> <br>

The output of  `phylosignal_network` is then:

| nb_A | nb_B | mantel_cor_A | pvalue_high_A | pvalue_low_A | 
| --- | --- | --- | --- | --- |
| 70 | 93 | 0.02 | 0.34 | 0.66  |

which corresponds to the number of orchid species (**nb_A**), the number of fungal species (**nb_B**), and the Mantel correlation between phylogenetic distances and degree difference distances for orchids (**mantel_cor_A**), its associated upper p-value (**pvalue_high_A**, *i.e.* the fraction of permutations that led to higher correlation values), and its associated lower p-value (**pvalue_low_A**, *i.e.* the fraction of permutations that led to lower correlation values).

Thus, here we do not detect any significant phylogenetic signal in the number of partners for the orchids (p-value>0.05). 



<br> <br>


# Optional steps:

## Option 1: Investigate clade-specific phylogenetic signals


This first option uses the function  `phylosignal_sub_network` to  compute the clade-specific phylogenetic signals in species interactions using simple Mantel tests with Bonferroni correction. For each node of the tree A having a certain number of descending species (*e.g.* 10), it computes the phylogenetic signal in the resulting sub-network by performing a Mantel test between the phylogenetic distances and the ecological distances for the given sub-clade of tree A. Mantel tests can be computed using quantified or binary networks, with the Jaccard or UniFrac ecological distances. The results of the clade-specific phylogenetic signal analysis can be represented using the function  `plot_phylosignal_sub_network`.

```r
# compute clade-specific phylogenetic signals in species interactions for orchids

results_clade_A <- phylosignal_sub_network(network, tree_A = tree_orchids, tree_B = tree_fungi, method = "GUniFrac", correlation = "Pearson", nperm=100000, minimum=10, degree=F)

plot_phylosignal_sub_network(tree_A = tree_orchids, results_clade_A)

```

| Option | Description |
| --- | --- |
| `network` | a matrix representing the bipartite interaction network with species from guild A in columns and species from guild B in rows. |
| `tree_A` | a phylogenetic tree of the guild A (the columns of the interaction network). |
| `tree_B` | a phylogenetic tree of the guild B (the rows of the interaction network). |
| `method` | indicates which method to use to compute the phylogenetic signal in species interactions: you can choose "Jaccard_weighted" for computing ecological distances using Jaccard dissimilarities (or "Jaccard_binary" to not take into account the abundances of the interactions), or "GUniFrac" to compute the weighted (or generalized) UniFrac distances (or "UniFrac_unweighted" to not take into account the interaction abundances). |
| `correlation` | indicates which correlation to use in the Mantel test, among the Pearson, Spearman, or Kendall correlations. |
| `nperm` | indicates the number of permutations to evaluate the significance of the Mantel test.  |
| `minimum` | indicates the minimal number of descending species for a node in tree A to compute its clade-specific phylogenetic signal.  |
| `degree` | if  `degree=TRUE `, Mantel tests testing for phylogenetic signal in the number of partners are additionally performed in each sub-clade.|
    



The output of  `phylosignal_sub_network` corresponds to a table where each line corresponds to a tested orchid sub-clade and which contains at least 8 columns: the name of the node (**node**), the number of species in the corresponding orchid sub-clade (**nb_A**), the number of fungal species  associated with the corresponding orchid sub-clade (**nb_B**), the Mantel correlation for the orchid sub-clade (**mantel_cor**), its associated upper p-value (**pvalue_high**), its associated lower p-value (**pvalue_low**), and the corresponding Bonferroni corrected p-values (**pvalue_high_corrected** and **pvalue_low_corrected**).

The representation of the results using `plot_phylosignal_sub_network` is a phylogenetic tree with nodes colored according to the clade-specific phylogenetic signals. Blue nodes are not significant (based in the Bonferonni correction), grey nodes are not tested (less than  `minimum` descending tips), and orange-red nodes represent significant phylogenetic signals and their color indicates the strength of the correlation.


<br> <br>



## Option 2: Test the robustness of the findings to phylogenetic uncertainty and/or sampling bias

This can be achieved by testing the robustness of the results to phylogenetic uncertainty (*e.g.* using a Bayesian tree posterior) or sampling bias (by subsampling over-represented species/clades). See Perez-Lamarque *et al.* (*in prep.*) for more details.



*PS :*  using the function  `phylosignal_network`, you can also apply the approach PBLM (Ives and Godfray, 2006) by specifying  `method="PBLM"`. However, given the high type-I error of this approach when testing for phylogenetic signal in species interactions, this is discouraged.
