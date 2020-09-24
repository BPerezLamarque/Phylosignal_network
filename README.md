
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
    <figcaption><center>A toy example of an interaction network between orchids (in green) and mycorrhizal fungi (in brown) informed with the phylogenetic trees of each guild.</center></figcaption>
</p>


First, load the example dataset of a mycorrhizal network between orchids and mycorrhizal fungi from La Réunion island (Martos *et al.*, 2012) along with the reconstructed phylogenetic trees of the orchids and the fungal OTUs:


```r
library(RPANDA)

data(mycorrhizal_network)

network <- mycorrhizal_network[[1]] # interaction matrix 
tree_orchids <- mycorrhizal_network[[2]] # phylogenetic tree (phylo object)
tree_fungi <- mycorrhizal_network[[3]] # phylogenetic tree (phylo object)

```




##  Step 1: test the phylogenetic signal in the species interactions (Mantel test)


The following function computes the phylogenetic signal in species interactions (do closely related species interact with similar partners?) using a simple Mantel test or the phylogenetic signal in the degree of generalism (do closely related species interact with the same number of partners?). Mantel tests measuring phylogenetic signal in species interactions can be computed using quantified or binary networks, with the Jaccard or UniFrac ecological distances.


<dl>
  <dt>Coffee</dt>
  <dd>- black hot drink</dd>
  <dt>Milk</dt>
  <dd>- white cold drink</dd>
</dl>




```r

# compute the phylogenetic signal in species intercations for orchids 
phylosignal_network(network, tree_orchids, tree_fungi, method = "GUniFrac", correlation = "Pearson")

```




