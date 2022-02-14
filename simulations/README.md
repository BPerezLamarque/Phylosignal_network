# Scripts and simulations 

This folder contains 2,400 bipartite interaction networks simulated using BipartiteEvol (Maliet et al. 2020) with the R-package RPANDA (Morlon et al. 2016; R Core Team 2020). Functions to simulate the BipartiteEvol interaction networks are available in [RPANDA](https://github.com/hmorlon/PANDA) (see function `sim.BipartiteEvol`).


We considered a total number of 500 (simul F), 1,000 (simul E), 2,000 (simul D), 3,000 (simul A), 4,000 (simul B), or 5,000 (simul C) pairs of interacting individuals per simulation. 

For each size, we simulated the evolution of: 
- 100 neutral networks (αA=0 ; αB=0)
- 120 mutualistic networks ((1) αA=1; αB=1; (2) αA=0.1; αB=0.1; (3) αA=0.01; αB=0.01; (4) αA=1; αB=0.1; (5) αA=1; αB=0.01; and (6) αA=0.1; αB=0.01) 
- 180 antagonistic networks ((1) αA=-1; αB=1; (2) αA=-0.1; αB=0.1; (3) αA=-0.01; αB=0.01; (4) αA=-1; αB=0.1; (5) αA=-1; αB=0.01; (6) αA=-0.1; αB=1; (7) αA=-0.1; αB=0.01; (8) αA=-0.01; αB=1; (9) αA=-0.01; αB=0.1). 
The number of the seed used to start the simulation is indicated for each individual network.

We used a mutation rate μ=0.01 and followed the interacting individuals during 50^6 death events. 

At the end, we extracted for each guild (A or B) a species tree from its genealogy by randomly selecting one individual per species (two trees per simulation are available in the Newick format) 
and reconstructed the corresponding weighted interaction network (a matrix in .csv format) by counting the number of occurrences of each interspecific interaction. 




**Citation:** Benoît Perez-Lamarque, Odile Maliet, Benoît Pichon, Marc-André Selosse, Florent Martos, and Hélène Morlon,
*Do closely related species interact with similar partners? Testing for phylogenetic signal in interaction networks*, bioRxiv, [doi: https://doi.org/10.1101/2021.08.30.458192](https://www.biorxiv.org/content/10.1101/2021.08.30.458192v1).


**Contact:** Benoît Perez-Lamarque, benoit.perez.lamarque@gmail.com


All the scripts for generating such networks and for measuring the phylogenetic signals in species interactions (using the methods described in Perez-Lamarque et al. 2021) are available in the folder `script/` (the file`script_phylogenetic_signal_network.R` is the main script).

