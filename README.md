
**Guidelines to measure the phylogenetic signal in a bipartite interaction network:**
====




This tutorial explains how to measure the phylogenetic signal in species interactions in a bipartite ecological network using the R-package RPANDA (Morlon *et al.*, 2016), *i.e.* how to test whether closely related species interact with similar partners.


**Citation:** Benoît Perez-Lamarque, Odile Maliet, Benoît Pichon, Marc-André Selosse, Florent Martos, and Hélène Morlon,
*Do closely related species interact with similar partners? Testing for phylogenetic signal in interaction networks*, bioRxiv, [doi: https://doi.org/10.1101/2021.08.30.458192](https://www.biorxiv.org/content/10.1101/2021.08.30.458192v3).


**Contact:** Benoît Perez-Lamarque, benoit.perez.lamarque@gmail.com


**Simulations:** Simulations used to investigate the performances of the different approaches measuring the phylogenetic signal in species interactions in bipartite interaction networks are available in the folder [simulations](https://github.com/BPerezLamarque/Phylosignal_network/tree/master/simulations). Functions to simulate the BipartiteEvol interaction networks are available in RPANDA, while the script used to simulate interaction networks with a phylogenetic signal in the number of partners (and not their identity) can be found [here](https://github.com/BPerezLamarque/Phylosignal_network/blob/master/simulations/script/sim_signal_number_partners.R).


# Contents:
**[Installation](#installation)**\
**[Measuring the phylogenetic signal](#measuring-phylogenetic-signal)**\
**[Optional steps](#optional-steps)**\
**[Other possibilities](#other-possibilities)**




# Installation:


The functions are part of the R-package RPANDA and can be installed from GitHub using devtools:

```r
library(devtools)
install_github("hmorlon/PANDA",ref="master", dependencies = TRUE)

```

(NB: RPANDA requires mpfr to be installed. If not already installed on your computer, you may get an error message. If so, please visit the [mpfr help page](https://www.mpfr.org/mpfr-current/mpfr.html) for installation). 

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



<p align="center">
    <img src="https://github.com/BPerezLamarque/Phylosignal_network/blob/master/example/figures_steps.png" width="500">
</p>

<p align="center">
    <b>Recommendations for measuring phylogenetic signals in species interactions.</b>
</p>


##  Step 1: Testing for the phylogenetic signal in the species interactions


This first step uses the function  `phylosignal_network` to compute the phylogenetic signal in species interactions (do closely related species interact with similar partners?) using a simple Mantel test. Mantel tests measuring the phylogenetic signal in species interactions can be computed using quantified or binary networks, with the Jaccard or UniFrac ecological distances.

```r

# compute phylogenetic signals in species interactions
phylosignal_network(network, tree_A = tree_orchids, tree_B = tree_fungi, method = "GUniFrac", correlation = "Pearson", nperm=10000)

```

| Option | Description |
| --- | --- |
| `network` | a matrix representing the bipartite interaction network with species from guild A in columns and species from guild B in rows. |
| `tree_A` | a phylogenetic tree of guild A (the columns of the interaction network). |
| `tree_B` | a phylogenetic tree of guild B (the rows of the interaction network). |
| `method` | indicates which method to use to compute the phylogenetic signal in species interactions: you can choose "Jaccard_weighted" for computing ecological distances using Jaccard dissimilarities (or "Jaccard_binary" to not take into account the abundances of the interactions), or "GUniFrac" to compute the weighted (or generalized) UniFrac distances (or "UniFrac_unweighted" to not take into account the interaction abundances). |
| `correlation` |indicates which correlation to use in the Mantel test, among the Pearson, Spearman, or Kendall correlations. |
| `nperm` | indicates the number of permutations to evaluate the significance of the Mantel test.  |

<br> <br>

The output of  `phylosignal_network` is:

| nb_A | nb_B | mantel_cor_A | pvalue_upper_A | pvalue_lower_A | mantel_cor_B | pvalue_upper_B | pvalue_lower_B |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 70 | 93 | -0.03 | 0.66 | 0.34  | 0.01  | 0.22 | 0.78 |

 which corresponds to the number of orchid species (**nb_A**), the number of fungal species (**nb_B**), the Mantel correlation  between the phylogenetic distances and ecological distances for orchids (**mantel_cor_A**), its associated upper p-value (**pvalue_upper_A**, *i.e.* the fraction of permutations that led to higher correlation values), its associated lower p-value (**pvalue_lower_A**, *i.e.* the fraction of permutations that led to lower correlation values), and  the Mantel correlation between the phylogenetic distances and ecological distances for fungi (**mantel_cor_B**), its associated upper p-value (**pvalue_upper_B**), abd its associated lower p-value (**pvalue_lower_B**),

Here, **pvalue_upper_A>0.05** so closely related orchid species do not tend to interact with similar mycorrhizal fungi. Similarly, **pvalue_upper_B>0.05** so closely related fungi do not tend to interact with similar orchids. Thus, we do not detect any significant phylogenetic signal in species interactions between the orchids and their mycorrhizal fungi (all p-values>0.05). 

<br> <br>

##  Step 2: Testing that the phylogenetic signal is still significant when accounting for the signal in the number of partners:


If there is a significant phylogenetic signal in Step 1, this second first step uses the function  `phylosignal_network` to compute the phylogenetic signal in species interactions using a simple Mantel test with permutations that keep constant the number of partners per species. Step 2 thus tests whether the phylogenetic signal (observed in Step 1) is still significant when accounting for the confounding effect of the phylogenetic signal in the number of partners.

```r

# compute phylogenetic signals in species interactions with permutations keeping constant the number fo partners 
phylosignal_network(network, tree_A = tree_orchids, tree_B = tree_fungi, method = "GUniFrac", correlation = "Pearson", nperm=1000, permutations="nbpartners")

```

| Option | Description |
| --- | --- |
| `network` | a matrix representing the bipartite interaction network with species from guild A in columns and species from guild B in rows. |
| `tree_A` | a phylogenetic tree of guild A (the columns of the interaction network). |
| `tree_B` | a phylogenetic tree of guild B (the rows of the interaction network). |
| `method` | indicates which method to use to compute the phylogenetic signal in species interactions: you can choose "Jaccard_weighted" for computing ecological distances using Jaccard dissimilarities (or "Jaccard_binary" to not take into account the abundances of the interactions), or "GUniFrac" to compute the weighted (or generalized) UniFrac distances (or "UniFrac_unweighted" to not take into account the interaction abundances). |
| `correlation` |indicates which correlation to use in the Mantel test, among the Pearson, Spearman, or Kendall correlations. |
| `nperm` | indicates the number of permutations to evaluate the significance of the Mantel test.  |
| `permutations` | indicates which permutations to compute either "shuffle" (*i.e.* random shuffling of the distance matrix, as in a regular Mantel test) or "nbpartners" (*i.e.* keeping constant the number of partners per species and shuffling at random their identity).  |

<br> <br>

The output of  `phylosignal_network` in Step 2 has the same correlation values as in Step 1. Only the **pvalues** are affected by the changes in the permutation strategy. 
For instance, if in Step 1 **pvalue_upper_A<0.05** and in  Step 2**pvalue_upper_A<0.05**, then we can conclude that there is a significant phylogenetic signal in species interactions that can not be fully explained by the phylogenetic signal in the number of partners. 
Alternatively,  if in Step 1 **pvalue_upper_A<0.05** and in  Step 2**pvalue_upper_A>0.05**, then we cannot exclude that the  phylogenetic signal in species interactions observed in Step 1 is not explained by the phylogenetic signal in the number of partners.

<br> <br>



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
| `tree_A` | a phylogenetic tree of guild A (the columns of the interaction network). |
| `tree_B` | a phylogenetic tree of guild B (the rows of the interaction network). |
| `method` | indicates which method to use to compute the phylogenetic signal in species interactions: you can choose "Jaccard_weighted" for computing ecological distances using Jaccard dissimilarities (or "Jaccard_binary" to not take into account the abundances of the interactions), or "GUniFrac" to compute the weighted (or generalized) UniFrac distances (or "UniFrac_unweighted" to not take into account the interaction abundances). |
| `correlation` | indicates which correlation to use in the Mantel test, among the Pearson, Spearman, or Kendall correlations. |
| `nperm` | indicates the number of permutations to evaluate the significance of the Mantel test.  |
| `minimum` | indicates the minimal number of descending species for a node in tree A to compute its clade-specific phylogenetic signal.  |
| `degree` | if  `degree=TRUE `, Mantel tests testing for phylogenetic signal in the number of partners are additionally performed in each sub-clade.|
    
<br> <br>

The output of  `phylosignal_sub_network` is then:


|node | nb_A | nb_B  | mantel_cor | pvalue_upper | pvalue_lower | pvalue_upper_corrected | pvalue_lower_corrected| 
| --- | --- | --- | --- | --- | --- | --- | --- |
| 71   | 70   | 93 | -0.03224387       | 0.670     |  0.331               |   1.000                | 1.000| 
| 72   | 16   | 27 | -0.07060626       | 0.639     |  0.362                |  1.000                | 1.000| 
| 73   | 15   | 26 | -0.20381036       | 0.954      | 0.047                |  1.000                | 0.752| 
| 74   | 13   | 26 | -0.05309340       | 0.488      | 0.513                |  1.000                | 1.000| 
| 75   | 12   | 24 |  0.04577263       | 0.297      | 0.704                 | 1.000                | 1.000| 
| 79   | 54   | 79  | 0.01193616       | 0.439      | 0.562                 | 1.000                | 1.000| 
| 80   | 53   | 78 | -0.12537619       | 0.924      | 0.077                 | 1.000              |   1.000| 
| 83   | 11   | 23 | -0.07965825       | 0.808      | 0.193                 | 1.000               |  1.000| 
| 87   | 39   | 56  | 0.03310306       | 0.368      | 0.633                 | 1.000                | 1.000| 
| 90   | 36   | 53  | 0.12048399       | 0.135      | 0.866                 | 1.000                | 1.000| 
| 92   | 34   | 53  | 0.36600645       | 0.001      | 1.000                 | **0.016**             |    1.000| 
| 98   | 28   | 44  | 0.15359285       | 0.005      | 0.996                 | 0.080                | 1.000| 
| 105   | 20   | 38  | 0.05440825      |  0.241     |  0.760                |  1.000               |  1.000| 
| 108   | 17   | 34 | -0.12135866      |  0.895     |  0.106                |  1.000               |  1.000| 
| 112   | 13   | 29 | -0.04986268      |  0.590     |  0.411                |  1.000               |  1.000| 
| 114   | 11   | 27  | 0.12542181       | 0.169      | 0.832                |  1.000                | 1.000| 

which corresponds to a table where each line corresponds to a tested orchid sub-clade and which contains at least 8 columns: the name of the node (**node**), the number of species in the corresponding orchid sub-clade (**nb_A**), the number of fungal species  associated with the corresponding orchid sub-clade (**nb_B**), the Mantel correlation for the orchid sub-clade (**mantel_cor**), its associated upper p-value (**pvalue_upper**), its associated lower p-value (**pvalue_lower**), and the corresponding Bonferroni corrected p-values (**pvalue_upper_corrected** and **pvalue_lower_corrected**).

The representation of the results using `plot_phylosignal_sub_network` is a phylogenetic tree with nodes colored according to the clade-specific phylogenetic signals. Blue nodes are not significant (based on the Bonferonni correction), grey nodes are not tested (less than `minimum` descending tips), and orange-red nodes represent significant phylogenetic signals and their color indicates the strength of the correlation.

<p align="center">
    <img src="https://github.com/BPerezLamarque/Phylosignal_network/blob/master/example/example_orchid_sub_clades.png" width="400">
</p>

<p align="center">
    <b> Clade-specific phylogenetic signals in species interactions in the orchids. The only significant phylogenetic signal (orange dot; pvalue_upper_corrected<0.05  ) corresponds to the Angraecineae.</b>
</p>


<br> <br>



## Option 2: Test the robustness of the findings to phylogenetic uncertainty and/or sampling bias

This can be achieved by testing the robustness of the results to phylogenetic uncertainty (*e.g.* using a Bayesian tree posterior) or sampling bias (by subsampling over-represented species/clades). See [Perez-Lamarque *et al.* (*in prep.*)](https://www.biorxiv.org/content/10.1101/2021.08.30.458192v1) for more details.


<br> <br>





# Other possibilities:





##  Testing for the phylogenetic signal in the number of partners


This  step uses the function  `phylosignal_network` to compute the phylogenetic signal in the number of partners (the degree; *i.e.* do closely related species interact with the same number of partners?) using a simple Mantel test. The goal of this step is to verify if there is a significant phylogenetic signal in species interactions (step 1), whether this phylogenetic signal in species interactions can be caused by a phylogenetic signal in the number of partners (confounding effect), rather than the identity of the interacting partners. We refer to this as a sequential Mantel test. 

```r

# compute the phylogenetic signal in the number of partners for orchids 
phylosignal_network(network, tree_A = tree_orchids, method = "degree", correlation = "Pearson", nperm=10000)


```

| Option | Description |
| --- | --- |
| `network` | a matrix representing the bipartite interaction network with species from guild A in columns and species from guild B in rows. |
| `tree_A` | a phylogenetic tree of guild A (the columns of the interaction network). |
| `correlation` | indicates which correlation to use in the Mantel test, among the Pearson, Spearman, or Kendall correlations. |
| `nperm` | indicates the number of permutations to evaluate the significance of the Mantel test.  |

<br> <br>

The output of  `phylosignal_network` is then:

| nb_A | nb_B | mantel_cor_A | pvalue_upper_A | pvalue_lower_A | 
| --- | --- | --- | --- | --- |
| 70 | 93 | 0.02 | 0.34 | 0.66  |

which corresponds to the number of orchid species (**nb_A**), the number of fungal species (**nb_B**), and the Mantel correlation between phylogenetic distances and degree difference distances for orchids (**mantel_cor_A**), its associated upper p-value (**pvalue_upper_A**, *i.e.* the fraction of permutations that led to higher correlation values), and its associated lower p-value (**pvalue_lower_A**, *i.e.* the fraction of permutations that led to lower correlation values).

Here, **pvalue_upper_A>0.05**, so we do not detect any significant phylogenetic signal in the number of partners for the orchids. 



##  Using the PBLM approach

Using the function  `phylosignal_network`, you can also apply the approach PBLM (Ives and Godfray, 2006) by specifying  `method="PBLM"`. However, given the very frequent false positives of this approach when testing for phylogenetic signals in species interactions, this is discouraged.



