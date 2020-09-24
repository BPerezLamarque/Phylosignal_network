
**Guidelines to measure phylogenetic signal in species interactions in a bipartite ecological network:**
====




This document indicates how to measure phylogenetic signal in species interactions in a bipartite ecological network using the R-package RPANDA (Morlon *et al.*, 2016).


**Citation:** Benoît Perez-Lamarque, Odile Maliet, Marc-André Selosse, Florent Martos, and Hélène Morlon, (*in prep.*)
**Do closely related species interact with similar partners? Testing for phylogenetic signal in ecological networks**.


**Contact:** Benoît Perez-Lamarque, benoit.perez.lamarque@gmail.com


# Contents:
**[Installation](#installation)**\
**[Measuring phylogenetic signal](#measuring-phylogenetic-signal)**\
**[Optionnal validations](#optionnal-validations):**\





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




##  Step 1: Test for phylogenetic signal in the species interactions (Mantel test)


The function  `phylosignal_network` computes the phylogenetic signal in species interactions (do closely related species interact with similar partners?) using a simple Mantel test or the phylogenetic signal in the degree of generalism (do closely related species interact with the same number of partners?). Mantel tests measuring phylogenetic signal in species interactions can be computed using quantified or binary networks, with the Jaccard or UniFrac ecological distances.


| Option | Description |
| --- | --- |
| `network` | a matrix representing the ecological interaction network with species from guild A in columns and species from guild B in rows. |
| `tree_A` | a phylogenetic tree of the guild A (the columns of the interaction network). |
| `tree_B` | a phylogenetic tree of the guild B (the rows of the interaction network). |
| `method` | indicates which method is used to compute the phylogenetic signal in species interactions: you can choose "Jaccard_weighted" for computing ecological distances using Jaccard dissimilarities (or "Jaccard_binary" to not take into account the abundances of the interactions), or "GUniFrac" to compute the generalized UniFrac distances ("UniFrac_unweighted" to not take into account the interaction abundances). |
| `correlation` |indicates which correlation is used in the Mantel test, among the Pearson, Spearman, or Kendall correlations. |
| `nperm` | number of permutations to evaluate the significance of the Mantel test.  |



```r

# compute the phylogenetic signal in species intercations for orchids 
phylosignal_network(network, tree_orchids, tree_fungi, method = "GUniFrac", correlation = "Pearson", nperm=10000)

```


The output of  `phylosignal_network` is:

| nb_A | nb_B | mantel_cor_A | pvalue_high_A | pvalue_low_A | mantel_cor_B | pvalue_high_B | pvalue_low_B |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 70 | 93 | -0.03 | 0.66 | 0.34  | 0.01  | 0.22 | 0.78 |

 which corresponds the number of orchid species (**nb_A**), the number of fungal species (**nb_B**), the Mantel correlation for orchids (**mantel_cor_A**), its associated upper p-value (**pvalue_high_A**), its associated lower p-value (**pvalue_low_A**), and  the Mantel correlation for fungi (**mantel_cor_B**), its associated upper p-value (**pvalue_high_B**), abd its associated lower p-value (**pvalue_low_B**),

Thus, here we do not detect any phylogenetic signal in species interactions. 




