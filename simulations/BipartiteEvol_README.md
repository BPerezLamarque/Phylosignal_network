


This folder contains 2,400 bipartite interaction networks simulated using BipartiteEvol (Maliet et al. 2020) with the R-package RPANDA (Morlon et al. 2016; R Core Team 2020). 

We considered a total number of 500 (simul F), 1,000 (simul E), 2,000 (simul D), 3,000 (simul A), 4,000 (simul B), or 5,000 (simul C) pairs of interacting individuals per simulation. 

For each size, we simulated the evolution of: 
- 100 neutral networks (αA=0 ; αB=0)
- 120 mutualistic networks ((1) αA=1; αB=1; (2) αA=0.1; αB=0.1; (3) αA=0.01; αB=0.01; (4) αA=1; αB=0.1; (5) αA=1; αB=0.01; and (6) αA=0.1; αB=0.01) 
- 180 antagonistic networks ((1) αA=-1; αB=1; (2) αA=-0.1; αB=0.1; (3) αA=-0.01; αB=0.01; (4) αA=-1; αB=0.1; (5) αA=-1; αB=0.01; (6) αA=-0.1; αB=1; (7) αA=-0.1; αB=0.01; (8) αA=-0.01; αB=1; (9) αA=-0.01; αB=0.1). 
The number of the seed used to start the simulation is indicated for each individual network.

We used a mutation rate μ=0.01 and followed the interacting individuals during 50^6 death events. 

At the end, we extracted for each guild (A or B) a species tree from its genealogy by randomly selecting one individual per species (two trees per simulations in the Newick format) 
and reconstructed the corresponding weighted interaction network (a matrix in .csv format) by counting the number of occurrences of each interspecific interaction. 


The scripts for generated such networks and for measuring the phylogenetic signals (methods as described in Perez-Lamarque et al. in prep) are available in the folder "script/". The file "script_phylogenetic_signal_network.R" is the main script.
