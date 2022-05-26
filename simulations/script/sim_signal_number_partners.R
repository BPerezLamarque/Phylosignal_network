
##### Simulate phylogenetic signal in the number of partners (but not in their identity) ####


rm(list=ls())

library(phytools)
library(mvMORPH)
library(reshape2)

seed=1

# choose a d_A value of the OU process 
# (high d_A corresponds to low phylogenetic signal, while d_A=0 corresponds to a Brownian motion, i.e. high phylogenetic signal)
d_A=5 
d_A=1
d_A=0.05
d_A=0

# simulate phylogenies

set.seed(3)
# number of species
nb_A <- floor(runif(1, 40, 150))
nb_B <- floor(runif(1, 40, 150))

# simulate trees
tree_A <- phytools::pbtree(n=nb_A)
tree_B <- phytools::pbtree(n=nb_B)

# scale branch length
tree_A$edge.length <- tree_A$edge.length/max(node.depth.edgelength(tree_A))
tree_B$edge.length <- tree_B$edge.length/max(node.depth.edgelength(tree_B))

# simulate OU with mean of 0 (theta=0), variance of 0.1 (noise of the brownian motion) and  a parameter of attraction toward 0, d_A (central tendency)
# phylogenetic signal in the traits (the number of partners) of clade A

if (d_A!="0"){
  data_A <- mvSIM(tree_A, param=list(sigma=0.1, alpha=d_A, ntraits=1, theta=0), model="OU1", nsim=1)
}else{ # brownian motion
  data_A <- mvSIM(tree_A, param=list(sigma=0.1, ntraits=1, theta=0), model="BM1", nsim=1)
}


# scale the data (from specialist to generalist according to the number of species in clade B)
data_A <- round( 1+(abs(data_A)-min(abs(data_A)))/(max(abs(data_A))-min(abs(data_A)))*(nb_B-1))

# make the network
network <- matrix(0, nrow=Ntip(tree_B), ncol= Ntip(tree_A))
colnames(network) <- tree_A$tip.label
rownames(network) <- tree_B$tip.label

# atribute the identity of the interactions at random
for (i in 1:nrow(data_A)){
  
  network[sample(1:nrow(network), size = data_A[i], replace = F), rownames(data_A)[i]] <- 1
  
}

colSums(network)==data_A
rowSums(network)

# Remove species if no interactions
network <- network[,colSums(network)>0]
network <- network[rowSums(network)>0,]
tree_B <- drop.tip(tree_B, tip=tree_B$tip.label[!tree_B$tip.label %in% rownames(network)])
tree_A <- drop.tip(tree_A, tip=tree_A$tip.label[!tree_A$tip.label %in% colnames(network)])


network <- network[tree_B$tip.label,tree_A$tip.label]


image(network) # final network

plot(tree_A) # final trees
plot(tree_B)








