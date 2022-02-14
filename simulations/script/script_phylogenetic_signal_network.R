

#### Step 1:  Simulate Bipartite network - simul 1 ####

rm(list=ls())


setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/script/")

# First: compile manually the C functions
# R CMD SHLIB fitness.c

library(coda)
library(ape)
library(fields)
library(Matrix)
library(RPANDA)
library(bipartite)
library(apTreeshape)
library(ade4)
library(parallel)

library(treebase, lib.loc="/users/biodiv/bperez/packages/")
library(adephylo, lib.loc="/users/biodiv/bperez/packages/")


dyn.load("fitness.so")
dyn.load("fitness_death.so")
source("ComputeMetricsEtablissement.R")
source("Rcode_comment_new.R") # compile c functions
source("beta.R")


setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1/")


####### Parameters do not vary 
plot=F

muP=0.01          # mutation rate clade P
muH=0.01          # mutation rate clade H

thresh=1  # species delineation

rP=Inf   # for all case 
rH=Inf

nP=50              # number of updated individuals at each time step (for each clade)
nH=50

NG=100000        # number of generation (you can add some afterwards if not long enough, see at the very end)
nb_it=9

thin = 100

seed=as.integer(1)

D=6

#### run the model ####
run_bipartite <- function(seed){
  
  if (!file.exists(paste0("network_",name,"_seed_",seed,".Rdata"))){
    print(paste0("Seed: ", seed))
    
    set.seed(seed)
    print(paste0("Iteration: ", 1))
    
    mod = model.spatial(nx,ny,NG=NG,D=D,muP=muP,muH=muH,alphaP=alphaP,alphaH=alphaH,iniP=0,iniH=0,nP=nP,nH=nH,rP=rP,rH=rH,effect=1,
                        verbose=1000,thin=thin,P=NULL,H=NULL,fitness_death=FALSE)
    
    #### build the genealogy ####
    gen = make.gen(mod, verbose=F)
    if (plot==T) plot(gen$H)
    if (plot==T) plot(gen$P)
    # not everything has coalesced yet, more iterations should be added -- see at the end
    
    #### build the phylogeny ####
    # with a species definition threshold s = 1
    phy1 = define.species(gen,threshold=thresh,verbose=F, monophyly = FALSE)
    if (plot==T) plot(phy1$Pphylo$tree)
    if (plot==T) plot(phy1$Hphylo$tree)
    
    #### build the network ####
    # still s = 1, change the phylogeny to change this parameter
    xP = mod$P       # traits of the individuals in clade P
    xH = mod$H
    
    net = build.network(xP, xH, gen, phy1, alpha = alphaP, seuil=0, rP=rP, rH=rH)
    # you can put a threshold and remove association with to low fitness with seuil, but I used 0 (the default, all associations are kept) in my case
    trait.id = 1
    if (plot==T) par(mfrow = c(1,1))
    if (plot==T) plot.network(net,phy1,1)
    
    network <- as.matrix(net)
    colnames(network) <- seq(1:ncol(network))
    rownames(network) <- seq(1:nrow(network))
    
    if (plot==T) image(log(network))
    
    write.table(network, paste0("network_",name,"_seed_",seed,".csv"), quote=F,sep=";")
    if (ncol(network)>1) write.tree(phy1$Hphylo$tree, paste0("network_tree_guild_A_",name,"_seed_",seed,".tre"))
    if (nrow(network)>1) write.tree(phy1$Pphylo$tree, paste0("network_tree_guild_B_",name,"_seed_",seed,".tre"))
    
    #### add iterations ####
    
    for (rep in 1:nb_it){
      print(paste0("Iteration: ", rep+1))
      
      # you can repeat the following
      mod = model.spatial(nx,ny,NG,D=D,muP,muH,alphaP=alphaP,alphaH=alphaH,iniP=0,iniH=0,nP=nP,nH=nH,rP=rP,rH=rH,effect=1,
                          verbose=1000,thin=thin,P=mod$P,H=mod$H,fitness_death=FALSE)
      
      # update the genealogy
      gen = make.gen(mod,treeP=gen$P, treeH=gen$H, verbose=F)
      if (plot==T) par(mfrow = c(1,1))
      if (plot==T) plot(gen$H)
      
      # update the phylogenies
      phy1 = define.species(gen,threshold=thresh, verbose=F, monophyly = FALSE)
      if (plot==T) par(mfrow = c(1,1))
      if (plot==T) plot(phy1$Pphylo$tree)
      if (plot==T) plot.model(gen,phy1, 1)
      
      xP = mod$P       # traits of the individuals in clade P
      xH = mod$H
      net = build.network(xP, xH, gen, phy1, alpha = alphaP, seuil=0, rP=rP, rH=rH)
      if (plot==T) plot.model.network(gen,phy1,trait.id, net,mod)
      
      network <- as.matrix(net)
      colnames(network) <- seq(1:ncol(network))
      rownames(network) <- seq(1:nrow(network))
      
      write.table(network, paste0("network_",name,"_seed_",seed,".csv"), quote=F,sep=";")
      if (ncol(network)>1) write.tree(phy1$Hphylo$tree, paste0("network_tree_guild_A_",name,"_seed_",seed,".tre"))
      if (nrow(network)>1) write.tree(phy1$Pphylo$tree, paste0("network_tree_guild_B_",name,"_seed_",seed,".tre"))
    }
    save(list = ls(),file = paste0("network_",name,"_seed_",seed,".Rdata"))
    print(noquote(c("End: ", seed)))
  }
}


# grid size
nx=80
ny=50

name="simul_1_A_neutral"
alphaH=0
alphaP=0
mclapply(1:20,run_bipartite,mc.cores = 20)

name="simul_1_A_mutualism_1"
alphaH=1
alphaP=1
mclapply(1:20,run_bipartite,mc.cores = 20)

name="simul_1_A_mutualism_2"
alphaH=0.1
alphaP=0.1
mclapply(1:20,run_bipartite,mc.cores = 20)

name="simul_1_A_mutualism_3"
alphaH=0.01
alphaP=0.01
mclapply(1:20,run_bipartite,mc.cores = 20)

name="simul_1_A_mutualism_4"
alphaH=1
alphaP=0.1
mclapply(1:20,run_bipartite,mc.cores = 20)

name="simul_1_A_mutualism_5"
alphaH=1
alphaP=0.01
mclapply(1:20,run_bipartite,mc.cores = 20)

name="simul_1_A_mutualism_6"
alphaH=0.1
alphaP=0.01
mclapply(1:20,run_bipartite,mc.cores = 20)

name="simul_1_A_antagonism_1"
alphaH=-1
alphaP=1
mclapply(1:20,run_bipartite,mc.cores = 20)

name="simul_1_A_antagonism_2"
alphaH=-0.1
alphaP=0.1
mclapply(1:20,run_bipartite,mc.cores = 20)

name="simul_1_A_antagonism_3"
alphaH=-0.01
alphaP=0.01
mclapply(1:20,run_bipartite,mc.cores = 20)

name="simul_1_A_antagonism_4"
alphaH=-1
alphaP=0.1
mclapply(1:20,run_bipartite,mc.cores = 20)

name="simul_1_A_antagonism_5"
alphaH=-1
alphaP=0.01
mclapply(1:20,run_bipartite,mc.cores = 20)

name="simul_1_A_antagonism_6"
alphaH=-0.1
alphaP=1
mclapply(1:20,run_bipartite,mc.cores = 20)

name="simul_1_A_antagonism_7"
alphaH=-0.1
alphaP=0.01
mclapply(1:20,run_bipartite,mc.cores = 20)

name="simul_1_A_antagonism_8"
alphaH=-0.01
alphaP=1
mclapply(1:20,run_bipartite,mc.cores = 20)

name="simul_1_A_antagonism_9"
alphaH=-0.01
alphaP=0.1
mclapply(1:20,run_bipartite,mc.cores = 20)



#####  Step 2-A: Compute phylogenetic signal - simul 1 ######


# on cluster
rm(list=ls())

library(vegan)
library(phytools)
library(ggplot2)
library(RPANDA,lib.loc="/users/biodiv/bperez/packages/")
library(GUniFrac,lib.loc="/users/biodiv/bperez/packages/")
library(parallel)
require(R.utils)

setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1")

transparent_theme <- theme(
  panel.grid = element_blank(),
  axis.line = element_line("black"),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))

source("../script/function_phylo_signal_network.R")

dyn.load("../script/permute.so")


param=c("neutral", paste0("mutualism_",1:6), paste0("antagonism_",1:9))

list_networks <- as.vector(outer(c("simul_1_A","simul_1_B","simul_1_C","simul_1_D","simul_1_E", "simul_1_F"), c(paste0("neutral_seed_",1:100), paste0("mutualism_", as.vector(outer(1:6, 1:20, paste, sep="_seed_"))), paste0("antagonism_",as.vector(outer(1:9, 1:20, paste, sep="_seed_")))), paste, sep="_"))

name = list_networks[1]


compute_phylo_signal <- function(name, method, correlation, nperm){
  
  if (!file.exists(paste0("results_phylogenetic_signal_",name,"_",method,".csv"))){
    
    network <- read.table(paste0("network_",name,".csv"), header=T,sep=";")
    
    colnames(network) <- seq(1:ncol(network))
    rownames(network) <- seq(1:nrow(network))
    
    tree_A <- read.tree(paste0("network_tree_guild_A_",name,".tre"))
    tree_B <- read.tree(paste0("network_tree_guild_B_",name,".tre"))
    
    
    # check that there are several species
    if ((nrow(network)>1)&(ncol(network)!=1)){
      
      # Not all species interactins with each others
      if (length(which(network==0))>0){
        
        output <- RPANDA::phylosignal_network(network, tree_A, tree_B, method = method, nperm = nperm, correlation = correlation)
        
        if (method %in% c("PBLM", "PBLM_binary")){
          names(output) <- c("nb_A", "nb_B", "dA", "dB", "MSETotal", "MSEFull", "MSEStar", "MSEBase")
          write.table(output,paste0("results_phylogenetic_signal_",name,"_",method,".csv"), quote=F, sep=";",row.names=F)}
        
        return(output)
        
      }else{
        return(c(Ntip(tree_A), Ntip(tree_B), 0,0,0,0,0,0))}
    }else{
      return(c(NA,NA, NA,NA,NA,NA,NA,NA))}
  }else{
    return(c(NA,NA,NA,NA,NA,NA,NA,NA))}
}


correlation="Pearson"

method="PBLM"
method="PBLM_binary"

method="Jaccard_weighted"
method="Jaccard_binary"
method="GUniFrac"
method="UniFrac_unweighted"

for (correlation in c("Pearson", "Spearman", "Kendall")){
  
  nperm=10000
  if (correlation=="Kendall") {nperm=100}
  
  results_signal <- matrix(unlist(mclapply(list_networks,compute_phylo_signal,mc.cores=20, mc.preschedule=F, method=method, correlation=correlation, nperm=nperm)),nrow=length(list_networks),byrow=T)
  results_signal <- cbind(list_networks, results_signal)
  results_signal <- data.frame(results_signal, stringsAsFactors =F)
  if (method %in% c("PBLM", "PBLM_binary")){colnames(results_signal) <- c("name","nb_A", "nb_B", "dA", "dB", "MSETotal", "MSEFull", "MSEStar", "MSEBase")
  }else{colnames(results_signal) <- c("name","nb_A", "nb_B", "mantel_cor_A","pvalue_high_A","pvalue_low_A", "mantel_cor_B", "pvalue_high_B", "pvalue_low_B")}
  
  results_signal$size <- "A"
  results_signal$size[grep(pattern = "_B_", results_signal$name)] <- "B"
  results_signal$size[grep(pattern = "_C_", results_signal$name)] <- "C"
  results_signal$size[grep(pattern = "_D_", results_signal$name)] <- "D"
  results_signal$size[grep(pattern = "_E_", results_signal$name)] <- "E"
  results_signal$size[grep(pattern = "_F_", results_signal$name)] <- "F"
  results_signal$param <- "neutral"
  results_signal$param[grep(pattern = "_antagonism_1_", results_signal$name)] <- "antagonism_1"
  results_signal$param[grep(pattern = "_antagonism_2_", results_signal$name)] <- "antagonism_2"
  results_signal$param[grep(pattern = "_antagonism_3_", results_signal$name)] <- "antagonism_3"
  results_signal$param[grep(pattern = "_antagonism_4_", results_signal$name)] <- "antagonism_4"
  results_signal$param[grep(pattern = "_antagonism_5_", results_signal$name)] <- "antagonism_5"
  results_signal$param[grep(pattern = "_antagonism_6_", results_signal$name)] <- "antagonism_6"
  results_signal$param[grep(pattern = "_antagonism_7_", results_signal$name)] <- "antagonism_7"
  results_signal$param[grep(pattern = "_antagonism_8_", results_signal$name)] <- "antagonism_8"
  results_signal$param[grep(pattern = "_antagonism_9_", results_signal$name)] <- "antagonism_9"
  results_signal$param[grep(pattern = "_mutualism_1_", results_signal$name)] <- "mutualism_1"
  results_signal$param[grep(pattern = "_mutualism_2_", results_signal$name)] <- "mutualism_2"
  results_signal$param[grep(pattern = "_mutualism_3_", results_signal$name)] <- "mutualism_3"
  results_signal$param[grep(pattern = "_mutualism_4_", results_signal$name)] <- "mutualism_4"
  results_signal$param[grep(pattern = "_mutualism_5_", results_signal$name)] <- "mutualism_5"
  results_signal$param[grep(pattern = "_mutualism_6_", results_signal$name)] <- "mutualism_6"
  
  if (!method %in% c("PBLM", "PBLM_binary")){write.table(results_signal,paste0("results_phylogenetic_signal_simul_1_",method,"_", correlation,".csv"), quote=F, sep=";",row.names=F)
  }else{write.table(results_signal,paste0("results_phylogenetic_signal_simul_1_",method,".csv"), quote=F, sep=";",row.names=F)}
  
}





#####  Step 2-B: Test compute PBLM with Aik - simul 1 ######


# on cluster
rm(list=ls())

library(vegan)
library(phytools)
library(ggplot2)
library(GUniFrac,lib.loc="/users/biodiv/bperez/packages/")
library(parallel)
require(R.utils)


setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1")

transparent_theme <- theme(
  panel.grid = element_blank(),
  axis.line = element_line("black"),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))



source("../script/function_phylo_signal_network.R")

dyn.load("../script/permute.so")


param=c("neutral", paste0("mutualism_",1:6), paste0("antagonism_",1:9))

list_networks <- as.vector(outer(c("simul_1_A","simul_1_B","simul_1_C","simul_1_D","simul_1_E", "simul_F"), c(paste0("neutral_seed_",1:100), paste0("mutualism_", as.vector(outer(1:6, 1:20, paste, sep="_seed_"))), paste0("antagonism_",as.vector(outer(1:9, 1:20, paste, sep="_seed_")))), paste, sep="_"))

name = list_networks[1]

compute_phylo_signal <- function(name, method, correlation, nperm){
  
  if (!file.exists(paste0("results_phylogenetic_signal_",name,"_",method,".csv"))){
    
    network <- read.table(paste0("network_",name,".csv"), header=T,sep=";")
    
    colnames(network) <- seq(1:ncol(network))
    rownames(network) <- seq(1:nrow(network))
    
    tree_A <- read.tree(paste0("network_tree_guild_A_",name,".tre"))
    tree_B <- read.tree(paste0("network_tree_guild_B_",name,".tre"))
    
    # check that there are several species
    if ((nrow(network)>1)&(ncol(network)!=1)){
      
      # Not all species interactins with each others
      if (length(which(network==0))>0){
        
        assoc_strength <- matrix(0, ncol=ncol(network), nrow=nrow(network) )
        rownames(assoc_strength) <- rownames(network)
        colnames(assoc_strength) <- colnames(network)
        
        for (i in 1:nrow(assoc_strength)){
          for (j in 1:ncol(assoc_strength)){
            assoc_strength[i,j] <- -log(1-network[i,j]/sum(network[i,]))
          }
        }
        
        assoc_strength[is.infinite(assoc_strength)] <- max(assoc_strength[is.finite(assoc_strength)])*10^3
        
        pstart=c(0.5,0.5)
        
        model_pblm <- R.utils::withTimeout(pblm(assocs=assoc_strength, tree1=tree_B, tree2=tree_A, bootstrap=F, nreps=0, pstart=pstart), timeout = 60*60*24, onTimeout = "silent")
        
        model_pblm$MSE
        model_pblm$signal.strength
        
        
        if (method %in% c("PBLM", "PBLM_binary")){
          names(output) <- c("nb_A", "nb_B", "dA", "dB", "MSETotal", "MSEFull", "MSEStar", "MSEBase")
          write.table(output,paste0("results_phylogenetic_signal_Aik_",name,"_",method,".csv"), quote=F, sep=";",row.names=F)}
        
        return(output)
        
      }else{
        return(c(nb_A, nb_B, 0,0,0,0,0,0))}
    }else{
      return(c(NA,NA, NA,NA,NA,NA,NA,NA))}
  }else{
    return(c(NA,NA,NA,NA,NA,NA,NA,NA))}
}


correlation="Pearson"

method="PBLM"
method="PBLM_binary"

method="Jaccard_weighted"
method="Jaccard_binary"
method="GUniFrac"
method="UniFrac_unweighted"

for (correlation in c("Pearson", "Spearman", "Kendall")){
  
  nperm=10000
  if (correlation=="Kendall") {nperm=100}
  
  results_signal <- matrix(unlist(mclapply(list_networks,compute_phylo_signal,mc.cores=20,mc.preschedule=F, method=method, correlation=correlation, nperm=nperm)),nrow=length(list_networks),byrow=T)
  results_signal <- cbind(list_networks, results_signal)
  results_signal <- data.frame(results_signal, stringsAsFactors =F)
  if (method %in% c("PBLM", "PBLM_binary")){colnames(results_signal) <- c("name","nb_A", "nb_B", "dA", "dB", "MSETotal", "MSEFull", "MSEStar", "MSEBase")
  }else{colnames(results_signal) <- c("name","nb_A", "nb_B", "mantel_cor_A","pvalue_high_A","pvalue_low_A", "mantel_cor_B", "pvalue_high_B", "pvalue_low_B")}
  
  results_signal$size <- "A"
  results_signal$size[grep(pattern = "_B_", results_signal$name)] <- "B"
  results_signal$size[grep(pattern = "_C_", results_signal$name)] <- "C"
  results_signal$size[grep(pattern = "_D_", results_signal$name)] <- "D"
  results_signal$size[grep(pattern = "_E_", results_signal$name)] <- "E"
  results_signal$size[grep(pattern = "_F_", results_signal$name)] <- "F"
  results_signal$param <- "neutral"
  results_signal$param[grep(pattern = "_antagonism_1_", results_signal$name)] <- "antagonism_1"
  results_signal$param[grep(pattern = "_antagonism_2_", results_signal$name)] <- "antagonism_2"
  results_signal$param[grep(pattern = "_antagonism_3_", results_signal$name)] <- "antagonism_3"
  results_signal$param[grep(pattern = "_antagonism_4_", results_signal$name)] <- "antagonism_4"
  results_signal$param[grep(pattern = "_antagonism_5_", results_signal$name)] <- "antagonism_5"
  results_signal$param[grep(pattern = "_antagonism_6_", results_signal$name)] <- "antagonism_6"
  results_signal$param[grep(pattern = "_antagonism_7_", results_signal$name)] <- "antagonism_7"
  results_signal$param[grep(pattern = "_antagonism_8_", results_signal$name)] <- "antagonism_8"
  results_signal$param[grep(pattern = "_antagonism_9_", results_signal$name)] <- "antagonism_9"
  results_signal$param[grep(pattern = "_mutualism_1_", results_signal$name)] <- "mutualism_1"
  results_signal$param[grep(pattern = "_mutualism_2_", results_signal$name)] <- "mutualism_2"
  results_signal$param[grep(pattern = "_mutualism_3_", results_signal$name)] <- "mutualism_3"
  results_signal$param[grep(pattern = "_mutualism_4_", results_signal$name)] <- "mutualism_4"
  results_signal$param[grep(pattern = "_mutualism_5_", results_signal$name)] <- "mutualism_5"
  results_signal$param[grep(pattern = "_mutualism_6_", results_signal$name)] <- "mutualism_6"
  
  if (!method %in% c("PBLM", "PBLM_binary")){write.table(results_signal,paste0("results_phylogenetic_signal_simul_1_",method,"_", correlation,".csv"), quote=F, sep=";",row.names=F)
  }else{write.table(results_signal,paste0("results_phylogenetic_signal_simul_1_",method,".csv"), quote=F, sep=";",row.names=F)}
  
}







#####  Step 2-C:  Compute PBLM with bootstrap - simul 1 F ######


# on cluster
rm(list=ls())

library(vegan)
library(phytools)
library(ggplot2)
library(GUniFrac,lib.loc="/users/biodiv/bperez/packages/")
library(parallel)
require(R.utils)

method="PBLM"
method="PBLM_binary"

setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1")

transparent_theme <- theme(
  panel.grid = element_blank(),
  axis.line = element_line("black"),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))


source("../script/function_phylo_signal_network.R")
dyn.load("../script/permute.so")



param=c("neutral", paste0("mutualism_",1:6), paste0("antagonism_",1:9))

list_networks <- as.vector(outer(c("simul_1_F"), c(paste0("neutral_seed_",1:100), paste0("mutualism_", as.vector(outer(1:6, 1:20, paste, sep="_seed_"))), paste0("antagonism_",as.vector(outer(1:9, 1:20, paste, sep="_seed_")))), paste, sep="_"))

name = list_networks[1]

compute_phylo_signal <- function(name, method){
  
  if (!file.exists(paste0("results_phylogenetic_signal_",name,"_",method,".csv"))){
    
    print(name)
    network <- read.table(paste0("network_",name,".csv"), header=T,sep=";")
    
    colnames(network) <- seq(1:ncol(network))
    rownames(network) <- seq(1:nrow(network))
    
    tree_A <- read.tree(paste0("network_tree_guild_A_",name,".tre"))
    tree_B <- read.tree(paste0("network_tree_guild_B_",name,".tre"))
    
    # check that there are several species
    if ((nrow(network)>1)&(ncol(network)!=1)){
      
      # Not all species interactions with each others
      if (length(which(network==0))>0){
        
        if (method=="PBLM_binary"){ network[network>0] <- 1 }
        
        # test bootstrap
        model_pblm <- R.utils::withTimeout(pblm(assocs=network, tree1=tree_B, tree2=tree_A, bootstrap=T, nreps=100), timeout = 60*60*24, onTimeout = "silent")
        model_pblm$MSE
        model_pblm$signal.strength
        
        output <- c(Ntip(tree_A), Ntip(tree_B), model_pblm$signal.strength[2,2], model_pblm$signal.strength[1,2], model_pblm$MSE ,model_pblm$signal.strength[2,1], model_pblm$signal.strength[2,3], model_pblm$signal.strength[1,1], model_pblm$signal.strength[1,3])
        
        
        if (method %in% c("PBLM", "PBLM_binary")){
          names(output) <- c("nb_A", "nb_B", "dA", "dB", "MSETotal", "MSEFull", "MSEStar", "MSEBase", "dA_min", "dA_max", "dB_min", "dB_max")
          write.table(output,paste0("results_phylogenetic_signal_PBLM_bootstrap_",name,"_",method,".csv"), quote=F, sep=";",row.names=F)}
        
        return(output)
        
      }else{
        return(c(nb_A, nb_B, 0,0,0,0,0,0,0,0,0,0))}
    }else{
      return(c(NA,NA, NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))}
  }else{
    return(c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA))}
}


nb_cores=15
set.seed(NULL)
results_signal <- matrix(unlist(mclapply(sample(list_networks),compute_phylo_signal,mc.cores=nb_cores,mc.preschedule=F, method=method)),nrow=length(list_networks),byrow=T)
results_signal <- cbind(list_networks, results_signal)
results_signal <- data.frame(results_signal, stringsAsFactors =F)
colnames(results_signal) <- c("name","nb_A", "nb_B", "dA", "dB", "MSETotal", "MSEFull", "MSEStar", "MSEBase", "dA_min", "dA_max", "dB_min", "dB_max")

results_signal$size <- "A"
results_signal$size[grep(pattern = "_B_", results_signal$name)] <- "B"
results_signal$size[grep(pattern = "_C_", results_signal$name)] <- "C"
results_signal$size[grep(pattern = "_D_", results_signal$name)] <- "D"
results_signal$size[grep(pattern = "_E_", results_signal$name)] <- "E"
results_signal$size[grep(pattern = "_F_", results_signal$name)] <- "F"
results_signal$param <- "neutral"
results_signal$param[grep(pattern = "_antagonism_1_", results_signal$name)] <- "antagonism_1"
results_signal$param[grep(pattern = "_antagonism_2_", results_signal$name)] <- "antagonism_2"
results_signal$param[grep(pattern = "_antagonism_3_", results_signal$name)] <- "antagonism_3"
results_signal$param[grep(pattern = "_antagonism_4_", results_signal$name)] <- "antagonism_4"
results_signal$param[grep(pattern = "_antagonism_5_", results_signal$name)] <- "antagonism_5"
results_signal$param[grep(pattern = "_antagonism_6_", results_signal$name)] <- "antagonism_6"
results_signal$param[grep(pattern = "_antagonism_7_", results_signal$name)] <- "antagonism_7"
results_signal$param[grep(pattern = "_antagonism_8_", results_signal$name)] <- "antagonism_8"
results_signal$param[grep(pattern = "_antagonism_9_", results_signal$name)] <- "antagonism_9"
results_signal$param[grep(pattern = "_mutualism_1_", results_signal$name)] <- "mutualism_1"
results_signal$param[grep(pattern = "_mutualism_2_", results_signal$name)] <- "mutualism_2"
results_signal$param[grep(pattern = "_mutualism_3_", results_signal$name)] <- "mutualism_3"
results_signal$param[grep(pattern = "_mutualism_4_", results_signal$name)] <- "mutualism_4"
results_signal$param[grep(pattern = "_mutualism_5_", results_signal$name)] <- "mutualism_5"
results_signal$param[grep(pattern = "_mutualism_6_", results_signal$name)] <- "mutualism_6"

write.table(results_signal,paste0("results_phylogenetic_signal_simul_1_",method,"_bootstrap.csv"), quote=F, sep=";",row.names=F)



#####  Step 3-A: Plot outputs phylogenetic signal - simul 1 - Mantel tests  ######


rm(list=ls())

library(ggplot2)

setwd("/Users/bperez/ownCloud/Recherche/These/ISYEB/signal_phylo/simul_1/")

transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))


eco_matrix="Jaccard_weighted"
eco_matrix="Jaccard_binary"
eco_matrix="GUniFrac"
eco_matrix="UniFrac_unweighted"

method="Pearson"

for (eco_matrix in c("Jaccard_weighted","Jaccard_binary","GUniFrac","UniFrac_unweighted")){
  
  for (method in c("Pearson", "Spearman", "Kendall")){
    
    results_signal <- read.table(paste0("results_phylogenetic_signal_simul_1_",eco_matrix,"_", method, ".csv"), header=T, sep=";")
    
    results_signal$mantel_cor_A <- as.numeric(results_signal$mantel_cor_A)
    results_signal$mantel_cor_B <- as.numeric(results_signal$mantel_cor_B)
    results_signal$pvalue_high_A <- as.numeric(results_signal$pvalue_high_A)
    results_signal$pvalue_high_B <- as.numeric(results_signal$pvalue_high_B)
    
    # plot size A and B 
    results_signal$size <- as.character(results_signal$size)
    results_signal$size[results_signal$size=="A"] <- 3000
    results_signal$size[results_signal$size=="B"] <- 4000
    results_signal$size[results_signal$size=="C"] <- 5000
    results_signal$size[results_signal$size=="D"] <- 2000
    results_signal$size[results_signal$size=="E"] <- 1000
    results_signal$size[grep(pattern = "_F_", results_signal$name)] <- 500
    
    
    results_signal$size <- as.factor(results_signal$size)
    results_signal$size <-  factor(results_signal$size,levels(results_signal$size)[c(5,1,2,3,4,6)])
    
    
    ggplot(results_signal,aes(x=param,y=nb_A,fill=size))+xlab("Parameters")+ylab("Number of species (A)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
      ylim(0, max(results_signal$nb_A))+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      geom_boxplot(alpha=0.9,size=0.75, color="#273746")
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_size_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 10, height = 5)
    
    ggplot(results_signal,aes(x=param,y=nb_B,fill=size))+xlab("Parameters")+ylab("Number of species (B)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
      ylim(0, max(results_signal$nb_B))+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      geom_boxplot(alpha=0.9,size=0.75, color="#273746")
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_size_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 10, height = 5)
    
    
    
    # Size plot
    
    ggplot(results_signal,aes(x=param,y=mantel_cor_A,fill=size))+xlab("Parameters")+ylab("Mantel correlation (A)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
      geom_hline(yintercept=0, linetype="dashed", color="#273746")+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      geom_boxplot(alpha=0.9,size=0.75, color="#273746")
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Mantel_tests_",eco_matrix,"_", method, "_size_correlation_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)
    
    ggplot(results_signal,aes(x=param,y=mantel_cor_B,fill=size))+xlab("Parameters")+ylab("Mantel correlation (B)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
      geom_hline(yintercept=0, linetype="dashed", color="#273746")+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      geom_boxplot(alpha=0.9,size=0.75, color="#273746")
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Mantel_tests_",eco_matrix,"_", method,"_size_correlation_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)
    
    
    ggplot(results_signal,aes(x=param,y=pvalue_high_A,fill=size))+xlab("Parameters")+ylab("pvalue - Mantel correlation (A)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
      geom_hline(yintercept=0, linetype="dashed", color="#273746")+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      geom_boxplot(alpha=0.9,size=0.75, color="#273746")
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Mantel_tests_",eco_matrix,"_", method,"_size_pvalue_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)
    
    ggplot(results_signal,aes(x=param,y=pvalue_high_B,fill=size))+xlab("Parameters")+ylab("pvalue - Mantel correlation (B)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
      geom_hline(yintercept=0, linetype="dashed", color="#273746")+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      geom_boxplot(alpha=0.9,size=0.75, color="#273746")
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Mantel_tests_",eco_matrix,"_", method,"_size_pvalue_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)
    
    
    
    ### No size
    
    ggplot(results_signal,aes(x=param,y=mantel_cor_A,fill=param))+xlab("Parameters")+ylab("Mantel correlation (A)") +labs(title=" ")+scale_fill_manual(values=c("#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#27ae60","#27ae60","#27ae60","#27ae60","#27ae60","#27ae60","#f4d03f"))+transparent_theme_no + 
      geom_hline(yintercept=0, linetype="dashed", color="#273746")+ theme(legend.position="none")+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      geom_boxplot(alpha=0.9,size=0.75, color="#273746")
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Mantel_tests_",eco_matrix,"_", method, "_correlation_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 6, height = 7)
    
    ggplot(results_signal,aes(x=param,y=mantel_cor_B,fill=param))+xlab("Parameters")+ylab("Mantel correlation (B)") +labs(title=" ")+scale_fill_manual(values=c("#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#27ae60","#27ae60","#27ae60","#27ae60","#27ae60","#27ae60","#f4d03f"))+transparent_theme_no + 
      geom_hline(yintercept=0, linetype="dashed", color="#273746")+ theme(legend.position="none")+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      geom_boxplot(alpha=0.9,size=0.75, color="#273746")
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Mantel_tests_",eco_matrix,"_", method,"_correlation_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 6, height = 7)
    
    
    ggplot(results_signal,aes(x=param,y=pvalue_high_A,fill=param))+xlab("Parameters")+ylab("pvalue - Mantel correlation (A)") +labs(title=" ")+scale_fill_manual(values=c("#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#27ae60","#27ae60","#27ae60","#27ae60","#27ae60","#27ae60","#f4d03f"))+transparent_theme_no + 
      geom_hline(yintercept=0, linetype="dashed", color="#273746")+ theme(legend.position="none")+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      geom_boxplot(alpha=0.9,size=0.75, color="#273746")
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Mantel_tests_",eco_matrix,"_", method,"_pvalue_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 6, height = 7)
    
    ggplot(results_signal,aes(x=param,y=pvalue_high_B,fill=param))+xlab("Parameters")+ylab("pvalue - Mantel correlation (B)") +labs(title=" ")+scale_fill_manual(values=c("#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#27ae60","#27ae60","#27ae60","#27ae60","#27ae60","#27ae60","#f4d03f"))+transparent_theme_no + 
      geom_hline(yintercept=0, linetype="dashed", color="#273746")+ theme(legend.position="none")+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      geom_boxplot(alpha=0.9,size=0.75, color="#273746")
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Mantel_tests_",eco_matrix,"_", method,"_pvalue_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 6, height = 7)
    
    
    ##  Resume for clade A
    
    results_signal$Inference <- "Not significant signal"
    results_signal$Inference[which(results_signal$pvalue_high_A<0.05)] <- "Significant signal"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.05))] <- "Significant signal (R>0.05)"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.15))] <- "Significant signal (R>0.15)"
    
    colors <- c()
    if ("Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
    if ("Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
    if ("Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
    if ("Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}
    
    library(dplyr)
    results_signal_all <- results_signal %>% 
      group_by(param,Inference,size) %>% 
      summarise(count=n()) %>% 
      mutate(perc=count/sum(count))
    
    brks <- c(0, 0.25, 0.5, 0.75, 1)
    ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
      xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
      scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
      scale_color_manual(values=colors)+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      facet_grid(~size)
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Mantel_tests_",eco_matrix,"_", method,"_inference_size_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 14, height = 7)
    
    
    results_signal$size_tot <- "<150"
    results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>=150)] <- "150-250"
    results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>250)] <- ">250"
    
    results_signal_all <- results_signal %>% 
      group_by(param,Inference,size_tot) %>% 
      summarise(count=n()) %>% 
      mutate(perc=count/sum(count))
    
    results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
    results_signal_all$size_tot <-  factor(results_signal_all$size_tot,levels(results_signal_all$size_tot)[c(1,3,2)])
    
    
    ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
      xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
      scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
      scale_color_manual(values=colors)+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      facet_grid(~size_tot)
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Mantel_tests_",eco_matrix,"_", method,"_inference_size_tot_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 10, height = 7)
    
    
    ### Add anticorrelation
    
    results_signal$Inference <- "d) Not significant signal"
    results_signal$Inference[which(results_signal$pvalue_high_A<0.05)] <- "e) Significant signal"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.05))] <- "f) Significant signal (R>0.05)"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.15))] <- "g) Significant signal (R>0.15)"
    results_signal$Inference[which(results_signal$pvalue_low_A<0.05)] <- "c) Significant anti-signal"
    results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.05)))] <- "b) Significant anti-signal (R<-0.05)"
    results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.15)))] <- "a) Significant anti-signal (R<-0.15)"
    
    #colors=c("#943126","#e74c3c","#f1948a","#f7dc6f","#7dcea0","#27ae60","#196f3d")
    
    colors <- c()
    if ("a) Significant anti-signal (R<-0.15)" %in% results_signal$Inference) {colors <- c(colors, "#943126")}
    if ("b) Significant anti-signal (R<-0.05)" %in% results_signal$Inference) {colors <- c(colors, "#e74c3c")}
    if ("c) Significant anti-signal" %in% results_signal$Inference) {colors <- c(colors, "#f1948a")}
    if ("d) Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
    if ("e) Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
    if ("f) Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
    if ("g) Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}
    
    results_signal_all <- results_signal %>% 
      group_by(param,Inference,size_tot) %>% 
      summarise(count=n()) %>% 
      mutate(perc=count/sum(count))
    
    results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
    results_signal_all$size_tot <-  factor(results_signal_all$size_tot,levels(results_signal_all$size_tot)[c(1,3,2)])
    
    ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
      xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
      scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
      scale_color_manual(values=colors)+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      facet_grid(~size_tot)
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Mantel_tests_",eco_matrix,"_", method,"_inference_size_tot_anticorrelation_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 11, height = 7)
    
    
    # no size
    results_signal_all <- results_signal %>% 
      group_by(param,Inference) %>% 
      summarise(count=n()) %>% 
      mutate(perc=count/sum(count))
    ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
      xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
      scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
      scale_color_manual(values=colors)+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Mantel_tests_",eco_matrix,"_", method,"_inference_new_anticorrelation_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 6, height = 5.5)
    
    
    
    ##  Resume for clade B
    
    results_signal$Inference <- "Not significant signal"
    results_signal$Inference[which(results_signal$pvalue_high_B<0.05)] <- "Significant signal"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.05))] <- "Significant signal (R>0.05)"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.15))] <- "Significant signal (R>0.15)"
    
    colors <- c()
    if ("Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
    if ("Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
    if ("Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
    if ("Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}
    
    library(dplyr)
    results_signal_all <- results_signal %>% 
      group_by(param,Inference,size) %>% 
      summarise(count=n()) %>% 
      mutate(perc=count/sum(count))
    
    brks <- c(0, 0.25, 0.5, 0.75, 1)
    ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
      xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
      scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
      scale_color_manual(values=colors)+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      facet_grid(~size)
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Mantel_tests_",eco_matrix,"_", method,"_inference_size_clade_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 14, height = 7)
    
    
    results_signal$size_tot <- "<150"
    results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>=150)] <- "150-250"
    results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>250)] <- ">250"
    
    results_signal_all <- results_signal %>% 
      group_by(param,Inference,size_tot) %>% 
      summarise(count=n()) %>% 
      mutate(perc=count/sum(count))
    
    results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
    results_signal_all$size_tot <-  factor(results_signal_all$size_tot,levels(results_signal_all$size_tot)[c(1,3,2)])
    
    
    ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
      xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
      scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
      scale_color_manual(values=colors)+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      facet_grid(~size_tot)
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Mantel_tests_",eco_matrix,"_", method,"_inference_size_tot_clade_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 10, height = 7)
    
    
    ### Add anticorrelation
    
    results_signal$Inference <- "d) Not significant signal"
    results_signal$Inference[which(results_signal$pvalue_high_B<0.05)] <- "e) Significant signal"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.05))] <- "f) Significant signal (R>0.05)"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.15))] <- "g) Significant signal (R>0.15)"
    results_signal$Inference[which(results_signal$pvalue_low_B<0.05)] <- "c) Significant anti-signal"
    results_signal$Inference[intersect(which(results_signal$pvalue_low_B<0.05), which(results_signal$mantel_cor_B<(-0.05)))] <- "b) Significant anti-signal (R<-0.05)"
    results_signal$Inference[intersect(which(results_signal$pvalue_low_B<0.05), which(results_signal$mantel_cor_B<(-0.15)))] <- "a) Significant anti-signal (R<-0.15)"
    
    #colors=c("#943126","#e74c3c","#f1948a","#f7dc6f","#7dcea0","#27ae60","#196f3d")
    colors <- c()
    if ("a) Significant anti-signal (R<-0.15)" %in% results_signal$Inference) {colors <- c(colors, "#943126")}
    if ("b) Significant anti-signal (R<-0.05)" %in% results_signal$Inference) {colors <- c(colors, "#e74c3c")}
    if ("c) Significant anti-signal" %in% results_signal$Inference) {colors <- c(colors, "#f1948a")}
    if ("d) Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
    if ("e) Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
    if ("f) Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
    if ("g) Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}
    
    results_signal_all <- results_signal %>% 
      group_by(param,Inference,size_tot) %>% 
      summarise(count=n()) %>% 
      mutate(perc=count/sum(count))
    
    results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
    results_signal_all$size_tot <-  factor(results_signal_all$size_tot,levels(results_signal_all$size_tot)[c(1,3,2)])
    
    ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
      xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
      scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
      scale_color_manual(values=colors)+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      facet_grid(~size_tot)
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Mantel_tests_",eco_matrix,"_", method,"_inference_size_tot_anticorrelation_clade_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 11, height = 7)
    
    
    # no size
    results_signal_all <- results_signal %>% 
      group_by(param,Inference) %>% 
      summarise(count=n()) %>% 
      mutate(perc=count/sum(count))
    ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
      xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
      scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
      scale_color_manual(values=colors)+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Mantel_tests_",eco_matrix,"_", method,"_inference_new_anticorrelation_clade_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 6, height = 5.5)
    
    
    
    print(eco_matrix)
    print(method)
    
    # type-I error
    print("A")
    print(length(intersect(which(results_signal$pvalue_high_A<0.05),intersect(which(results_signal$param=="neutral"),which(results_signal$size_tot=="<150"))))/length(intersect(which(results_signal$param=="neutral"),which(results_signal$size_tot=="<150"))))
    print(length(intersect(which(results_signal$pvalue_high_A<0.05),intersect(which(results_signal$param=="neutral"),which(results_signal$size_tot=="150-250"))))/length(intersect(which(results_signal$param=="neutral"),which(results_signal$size_tot=="150-250"))))
    print(length(intersect(which(results_signal$pvalue_high_A<0.05),intersect(which(results_signal$param=="neutral"),which(results_signal$size_tot==">250"))))/length(intersect(which(results_signal$param=="neutral"),which(results_signal$size_tot==">250"))))
    
    print("B")
    print(length(intersect(which(results_signal$pvalue_high_B<0.05),intersect(which(results_signal$param=="neutral"),which(results_signal$size_tot=="<150"))))/length(intersect(which(results_signal$param=="neutral"),which(results_signal$size_tot=="<150"))))
    print(length(intersect(which(results_signal$pvalue_high_B<0.05),intersect(which(results_signal$param=="neutral"),which(results_signal$size_tot=="150-250"))))/length(intersect(which(results_signal$param=="neutral"),which(results_signal$size_tot=="150-250"))))
    print(length(intersect(which(results_signal$pvalue_high_B<0.05),intersect(which(results_signal$param=="neutral"),which(results_signal$size_tot==">250"))))/length(intersect(which(results_signal$param=="neutral"),which(results_signal$size_tot==">250"))))
    
    
  }
}





#####  Step 3-B: Plot outputs phylogenetic signal - simul 1 - PBLM with MSE  ######

# on cluster

rm(list=ls())

library(ggplot2)

setwd("/Users/bperez/ownCloud/Recherche/These/ISYEB/signal_phylo/simul_1/")
setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1")


transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))


eco_matrix="PBLM"
eco_matrix="PBLM_binary"

for (eco_matrix in c( "PBLM", "PBLM_binary")){
  
  info = list.files(path = "results_PBLM/", pattern = paste0('_',eco_matrix,'.csv'))
  info <- info[grep(x=info,pattern = "seed")]
  
  results_all <- c()
  for (i in info){
    
    if (length(grep(x=i, pattern = "_neutral_")==1)>0) {param="neutral"}
    if (length(grep(x=i, pattern = "_mutualism_1_")==1)>0) {param="mutualism_1"}
    if (length(grep(x=i, pattern = "_mutualism_2_")==1)>0) {param="mutualism_2"}
    if (length(grep(x=i, pattern = "_mutualism_3_")==1)>0) {param="mutualism_3"}
    if (length(grep(x=i, pattern = "_mutualism_4_")==1)>0) {param="mutualism_4"}
    if (length(grep(x=i, pattern = "_mutualism_5_")==1)>0) {param="mutualism_5"}
    if (length(grep(x=i, pattern = "_mutualism_6_")==1)>0) {param="mutualism_6"}
    if (length(grep(x=i, pattern = "_antagonism_1_")==1)>0) {param="antagonism_1"}
    if (length(grep(x=i, pattern = "_antagonism_2_")==1)>0) {param="antagonism_2"}
    if (length(grep(x=i, pattern = "_antagonism_3_")==1)>0) {param="antagonism_3"}
    if (length(grep(x=i, pattern = "_antagonism_4_")==1)>0) {param="antagonism_4"}
    if (length(grep(x=i, pattern = "_antagonism_5_")==1)>0) {param="antagonism_5"}
    if (length(grep(x=i, pattern = "_antagonism_6_")==1)>0) {param="antagonism_6"}
    if (length(grep(x=i, pattern = "_antagonism_7_")==1)>0) {param="antagonism_7"}
    if (length(grep(x=i, pattern = "_antagonism_8_")==1)>0) {param="antagonism_8"}
    if (length(grep(x=i, pattern = "_antagonism_9_")==1)>0) {param="antagonism_9"}
    
    if (length(grep(x=i, pattern = "_A_")==1)>0) {size="3000"}
    if (length(grep(x=i, pattern = "_B_")==1)>0) {size="4000"}
    if (length(grep(x=i, pattern = "_C_")==1)>0) {size="5000"}
    if (length(grep(x=i, pattern = "_D_")==1)>0) {size="2000"}
    if (length(grep(x=i, pattern = "_E_")==1)>0) {size="1000"}
    if (length(grep(x=i, pattern = "_F_")==1)>0) {size="500"}
    
    results_seed <- read.table(paste0("results_PBLM/",i), header=T, sep=";")
    
    name_simul <- as.character(unlist(gsub("_PBLM.csv", "", gsub("_PBLM_binary.csv", "", gsub("results_phylogenetic_signal_", "", i)))))
    
    # pb with D and E
    if (ncol(results_seed)==1) {results_seed = t(results_seed)}
    
    results_all <- rbind(results_all, c(name_simul, param, size, results_seed))
    
  }
  
  results_all <- data.frame(results_all)
  
  colnames(results_all) <- c("name", "param","size", "nb_A", "nb_B", "dA", "dB", "MSETotal", "MSEFull", "MSEStar", "MSEBase")
  
  results_all$param <- as.character(results_all$param)
  results_all$size <- as.character(results_all$size)
  results_all$nb_A <- as.numeric(as.character(results_all$nb_A))
  results_all$nb_B <- as.numeric(as.character(results_all$nb_B))
  results_all$dA <- as.numeric(as.character(results_all$dA))
  results_all$dB <- as.numeric(as.character(results_all$dB))
  results_all$MSETotal <- as.numeric(as.character(results_all$MSETotal))
  results_all$MSEFull <- as.numeric(as.character(results_all$MSEFull))
  results_all$MSEStar <- as.numeric(as.character(results_all$MSEStar))
  results_all$MSEBase <- as.numeric(as.character(results_all$MSEBase))
  
  results_all$size <- as.factor(results_all$size)
  results_all$size <-  factor(results_all$size,levels(results_all$size)[c(5,1,2,3,4,6)])
  
  
  results_all$Diff_MSE <- results_all$MSEStar - results_all$MSEFull
  
  results_all$Ratio_MSE <- (results_all$MSEStar - results_all$MSEFull)/results_all$MSEStar
  
  results_all$name <- unlist(results_all$name)
  
  write.table(results_all,paste0("results_phylogenetic_signal_simul_1_",eco_matrix,".csv"), quote=F, sep=";",row.names=F)
  
  
  results_all_failed <- results_all
  results_all <- results_all[which(results_all$dA!="inference failed"),]
  
  # Plot 
  
  ggplot(results_all,aes(x=param,y=dA,fill=size))+xlab("Parameters")+ylab("Phylogenetic signal (dA)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
    geom_hline(yintercept=0, linetype="dashed", color="#273746")+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    geom_boxplot(alpha=0.9,size=0.75, color="#273746")
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_",eco_matrix,"_size_dA.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)
  
  
  
  ggplot(results_all,aes(x=param,y=dB,fill=size))+xlab("Parameters")+ylab("Phylogenetic signal (dB)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
    geom_hline(yintercept=0, linetype="dashed", color="#273746")+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    geom_boxplot(alpha=0.9,size=0.75, color="#273746")
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_",eco_matrix,"_size_dB.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)
  
  
  
  ggplot(results_all,aes(x=param,y=Diff_MSE,fill=size))+xlab("Parameters")+ylab("MSE(Star) - MSE(Model)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
    geom_hline(yintercept=0, linetype="dashed", color="#273746")+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    geom_boxplot(alpha=0.9,size=0.75, color="#273746")
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_",eco_matrix,"_size_difference_MSE.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)
  
  
  ggplot(results_all,aes(x=param,y=Ratio_MSE,fill=size))+xlab("Parameters")+ylab("(MSE(Star) - MSE(Model))/MSE(Star)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
    geom_hline(yintercept=0, linetype="dashed", color="#273746")+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    geom_boxplot(alpha=0.9,size=0.75, color="#273746")
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_",eco_matrix,"_size_difference_MSE_ratio.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)
  
  
  
  ## No Size 
  
  ggplot(results_all,aes(x=param,y=dA,fill=param))+xlab("Parameters")+ylab("Phylogenetic signal (dA)") +labs(title=" ")+scale_fill_manual(values=c("#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#27ae60","#27ae60","#27ae60","#27ae60","#27ae60","#27ae60","#f4d03f"))+transparent_theme_no + 
    geom_hline(yintercept=0, linetype="dashed", color="#273746")+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    geom_boxplot(alpha=0.9,size=0.75, color="#273746")
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_",eco_matrix,"_dA.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)
  
  
  
  ggplot(results_all,aes(x=param,y=dB,fill=param))+xlab("Parameters")+ylab("Phylogenetic signal (dB)") +labs(title=" ")+scale_fill_manual(values=c("#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#27ae60","#27ae60","#27ae60","#27ae60","#27ae60","#27ae60","#f4d03f"))+transparent_theme_no + 
    geom_hline(yintercept=0, linetype="dashed", color="#273746")+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    geom_boxplot(alpha=0.9,size=0.75, color="#273746")
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_",eco_matrix,"_dB.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)
  
  
  
  ggplot(results_all,aes(x=param,y=Diff_MSE,fill=param))+xlab("Parameters")+ylab("MSE(Star) - MSE(Model)") +labs(title=" ")+scale_fill_manual(values=c("#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#27ae60","#27ae60","#27ae60","#27ae60","#27ae60","#27ae60","#f4d03f"))+transparent_theme_no + 
    geom_hline(yintercept=0, linetype="dashed", color="#273746")+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    geom_boxplot(alpha=0.9,size=0.75, color="#273746")
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_",eco_matrix,"_difference_MSE.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)
  
  
  ggplot(results_all,aes(x=param,y=Ratio_MSE,fill=param))+xlab("Parameters")+ylab("(MSE(Star) - MSE(Model))/MSE(Star)") +labs(title=" ")+scale_fill_manual(values=c("#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#d35400","#27ae60","#27ae60","#27ae60","#27ae60","#27ae60","#27ae60","#f4d03f"))+transparent_theme_no + 
    geom_hline(yintercept=0, linetype="dashed", color="#273746")+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    geom_boxplot(alpha=0.9,size=0.75, color="#273746")
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_",eco_matrix,"_difference_MSE_ratio.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)
  
  
  results_all_failed$Inference <- "Inference failed"
  results_all_failed$Inference[which(results_all_failed$Ratio_MSE<=0)] <- "a) Not significant signal"
  results_all_failed$Inference[which(results_all_failed$Ratio_MSE>0)] <- "b) Significant signal"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dA>0.05))] <- "c) Significant signal (d>0.05)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dB>0.05))] <- "c) Significant signal (d>0.05)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dA>0.15))] <- "d) Significant signal (d>0.15)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dB>0.15))] <- "d) Significant signal (d>0.15)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dA>1))] <- "e) Significant anti-signal (d>1)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dB>1))] <- "e) Significant anti-signal (d>1)"
  
  colors <- c()
  if ("a) Not significant signal" %in% results_all_failed$Inference) {colors <- c(colors, "#f7dc6f")}
  if ("b) Significant signal" %in% results_all_failed$Inference) {colors <- c(colors, "#7dcea0")}
  if ("c) Significant signal (d>0.05)" %in% results_all_failed$Inference) {colors <- c(colors, "#27ae60")}
  if ("d) Significant signal (d>0.15)" %in% results_all_failed$Inference) {colors <- c(colors, "#196f3d")}
  if ("e) Significant anti-signal (d>1)" %in% results_all_failed$Inference) {colors <- c(colors, "#943126")}
  
  
  library(dplyr)
  results_all_data <- results_all_failed %>% 
    group_by(param,Inference,size) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  brks <- c(0, 0.25, 0.5, 0.75, 1)
  ggplot(results_all_data,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=colors)+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    facet_grid(~size)
  
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_",eco_matrix,"_inference_size.pdf"), device=cairo_pdf, plot = last_plot(), width = 14, height = 7)
  
  
  ### More stringent
  
  results_all_failed$Inference <- "Inference failed"
  results_all_failed$Inference[which(results_all_failed$Ratio_MSE<=0.01)] <- "a) Not significant signal"
  results_all_failed$Inference[which(results_all_failed$Ratio_MSE>0.01)] <- "b) Significant signal"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0.01), which(results_all_failed$dA>0.05))] <- "c) Significant signal (d>0.05)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0.01), which(results_all_failed$dB>0.05))] <- "c) Significant signal (d>0.05)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0.01), which(results_all_failed$dA>0.15))] <- "d) Significant signal (d>0.15)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0.01), which(results_all_failed$dB>0.15))] <- "d) Significant signal (d>0.15)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0.01), which(results_all_failed$dA>1))] <- "e) Significant anti-signal (d>1)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0.01), which(results_all_failed$dB>1))] <- "e) Significant anti-signal (d>1)"
  
  
  colors <- c()
  if ("a) Not significant signal" %in% results_all_failed$Inference) {colors <- c(colors, "#f7dc6f")}
  if ("b) Significant signal" %in% results_all_failed$Inference) {colors <- c(colors, "#7dcea0")}
  if ("c) Significant signal (d>0.05)" %in% results_all_failed$Inference) {colors <- c(colors, "#27ae60")}
  if ("d) Significant signal (d>0.15)" %in% results_all_failed$Inference) {colors <- c(colors, "#196f3d")}
  if ("e) Significant anti-signal (d>1)" %in% results_all_failed$Inference) {colors <- c(colors, "#943126")}
  
  
  
  library(dplyr)
  results_all_data <- results_all_failed %>% 
    group_by(param,Inference,size) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  brks <- c(0, 0.25, 0.5, 0.75, 1)
  ggplot(results_all_data,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=colors)+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    facet_grid(~size)
  
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_",eco_matrix,"_inference_size_stringent.pdf"), device=cairo_pdf, plot = last_plot(), width = 14, height = 7)
  
  results_all_failed$Inference <- "Inference failed"
  results_all_failed$Inference[which(results_all_failed$Ratio_MSE<=0)] <- "a) Not significant signal"
  results_all_failed$Inference[which(results_all_failed$Ratio_MSE>0)] <- "b) Significant signal"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dA>0.05))] <- "c) Significant signal (d>0.05)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dB>0.05))] <- "c) Significant signal (d>0.05)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dA>0.15))] <- "d) Significant signal (d>0.15)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dB>0.15))] <- "d) Significant signal (d>0.15)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dA>1))] <- "e) Significant anti-signal (d>1)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dB>1))] <- "e) Significant anti-signal (d>1)"
  
  results_all_failed$size_tot <- "<150"
  results_all_failed$size_tot[which(results_all_failed$nb_A+results_all_failed$nb_B>=150)] <- "150-250"
  results_all_failed$size_tot[which(results_all_failed$nb_A+results_all_failed$nb_B>250)] <- ">250"
  
  colors <- c()
  if ("a) Not significant signal" %in% results_all_failed$Inference) {colors <- c(colors, "#f7dc6f")}
  if ("b) Significant signal" %in% results_all_failed$Inference) {colors <- c(colors, "#7dcea0")}
  if ("c) Significant signal (d>0.05)" %in% results_all_failed$Inference) {colors <- c(colors, "#27ae60")}
  if ("d) Significant signal (d>0.15)" %in% results_all_failed$Inference) {colors <- c(colors, "#196f3d")}
  if ("e) Significant anti-signal (d>1)" %in% results_all_failed$Inference) {colors <- c(colors, "#943126")}
  
  
  results_signal_all <- results_all_failed %>% 
    group_by(param,Inference,size_tot) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
  results_signal_all$size_tot <-  factor(results_signal_all$size_tot,c("<150"  ,  "150-250", ">250"  ))
  
  
  ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=colors)+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    facet_grid(~size_tot)
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_",eco_matrix,"_inference_size_tot.pdf"), device=cairo_pdf, plot = last_plot(), width = 10, height = 7)
  
  # no size
  results_signal_all <- results_all_failed %>% 
    group_by(param,Inference) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=colors)+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_",eco_matrix,"_inference_new.pdf"), device=cairo_pdf, plot = last_plot(), width = 5.6, height = 5.5)
  
  
  # stringent 5%
  results_all_failed$Inference <- "Inference failed"
  results_all_failed$Inference[which(results_all_failed$Ratio_MSE<=0.05)] <- "a) Not significant signal"
  results_all_failed$Inference[which(results_all_failed$Ratio_MSE>0.05)] <- "b) Significant signal"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0.05), which(results_all_failed$dA>0.05))] <- "c) Significant signal (d>0.05)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0.05), which(results_all_failed$dB>0.05))] <- "c) Significant signal (d>0.05)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0.05), which(results_all_failed$dA>0.15))] <- "d) Significant signal (d>0.15)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0.05), which(results_all_failed$dB>0.15))] <- "d) Significant signal (d>0.15)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0.05), which(results_all_failed$dA>1))] <- "e) Significant anti-signal (d>1)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0.05), which(results_all_failed$dB>1))] <- "e) Significant anti-signal (d>1)"
  
  
  results_all_failed$size_tot <- "<150"
  results_all_failed$size_tot[which(results_all_failed$nb_A+results_all_failed$nb_B>=150)] <- "150-250"
  results_all_failed$size_tot[which(results_all_failed$nb_A+results_all_failed$nb_B>250)] <- ">250"
  
  colors <- c()
  if ("a) Not significant signal" %in% results_all_failed$Inference) {colors <- c(colors, "#f7dc6f")}
  if ("b) Significant signal" %in% results_all_failed$Inference) {colors <- c(colors, "#7dcea0")}
  if ("c) Significant signal (d>0.05)" %in% results_all_failed$Inference) {colors <- c(colors, "#27ae60")}
  if ("d) Significant signal (d>0.15)" %in% results_all_failed$Inference) {colors <- c(colors, "#196f3d")}
  if ("e) Significant anti-signal (d>1)" %in% results_all_failed$Inference) {colors <- c(colors, "#943126")}
  
  results_signal_all <- results_all_failed %>% 
    group_by(param,Inference,size_tot) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
  results_signal_all$size_tot <-  factor(results_signal_all$size_tot,c("<150"  ,  "150-250", ">250"  ))
  
  
  ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=colors)+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    facet_grid(~size_tot)
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_",eco_matrix,"_inference_size_tot_stringeant_5.pdf"), device=cairo_pdf, plot = last_plot(), width = 10, height = 7)
  
  
  
  ##### Add signal traits 
  
  results_all_failed$Inference <- "Inference failed"
  results_all_failed$Inference[which(results_all_failed$Ratio_MSE<=0)] <- "a) Not significant signal"
  results_all_failed$Inference[which(results_all_failed$Ratio_MSE>0)] <- "b) Significant signal"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dA>0.05))] <- "c) Significant signal (d>0.05)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dB>0.05))] <- "c) Significant signal (d>0.05)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dA>0.15))] <- "d) Significant signal (d>0.15)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dB>0.15))] <- "d) Significant signal (d>0.15)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dA>1))] <- "e) Significant anti-signal (d>1)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dB>1))] <- "e) Significant anti-signal (d>1)"
  
  
  results_traits <- read.table(paste0("results_phylogenetic_signal_simul_1_in_traits.csv"), sep=";", header=T)
  
  results_traits$param <- as.character(results_traits$param)
  results_traits$mantel_cor_A <- as.numeric(as.character(results_traits$mantel_cor_A))
  results_traits$pvalue_high_A <- as.numeric(as.character(results_traits$pvalue_high_A))
  results_traits$mantel_cor_B <- as.numeric(as.character(results_traits$mantel_cor_B))
  results_traits$pvalue_high_B <- as.numeric(as.character(results_traits$pvalue_high_B))
  results_traits$nb_A <- as.numeric(as.character(results_traits$nb_A))
  results_traits$nb_B <- as.numeric(as.character(results_traits$nb_B))
  results_traits$name_simul_type <- paste0(results_traits$name, "_", results_traits$param, "_seed_", results_traits$seed )
  
  results_traits$Inference_traits_B <- "Not significant signal"
  results_traits$Inference_traits_B[which(results_traits$pvalue_high_B<0.05)] <- "Significant signal"
  results_traits$Inference_traits_B[which(results_traits$pvalue_low_B<0.05)] <- "Significant anti-signal"
  
  results_traits$Inference_traits_A <- "Not significant signal"
  results_traits$Inference_traits_A[which(results_traits$pvalue_high_A<0.05)] <- "Significant signal"
  results_traits$Inference_traits_A[which(results_traits$pvalue_low_A<0.05)] <- "Significant anti-signal"
  
  
  # make correspondance
  
  results_traits$name_simul_type <- as.character(results_traits$name_simul_type)
  results_all_failed$name <- as.character(results_all_failed$name)
  
  for (i in 1:nrow(results_all_failed)){
    index <- which(results_traits$name_simul_type==results_all_failed$name[i])
    results_all_failed$Inference_traits_A[i] <- results_traits$Inference_traits_A[index]
    results_all_failed$Inference_traits_B[i] <- results_traits$Inference_traits_B[index]
    results_all_failed$mantel_cor_traits_A[i] <- results_traits$mantel_cor_A[index]
    results_all_failed$mantel_cor_traits_B[i] <- results_traits$mantel_cor_B[index]
  }
  
  # expected signal
  results_all_failed$Expected_signal_A <- "Expected"
  results_all_failed$Expected_signal_A[which(results_all_failed$Inference_traits_A=="Significant anti-signal")] <- "Expected (negative)"
  results_all_failed$Expected_signal_A[which(results_all_failed$Inference_traits_A=="Not significant signal")] <- "Non-expected (traits)"
  results_all_failed$Expected_signal_A[grep("neutral", results_all_failed$name)] <- "Non-expected (neutral)"
  
  results_all_failed$Expected_signal_B <- "Expected"
  results_all_failed$Expected_signal_B[which(results_all_failed$Inference_traits_B=="Significant anti-signal")] <- "Expected (negative)"
  results_all_failed$Expected_signal_B[which(results_all_failed$Inference_traits_B=="Not significant signal")] <- "Non-expected (traits)"
  results_all_failed$Expected_signal_B[grep("neutral", results_all_failed$name)] <- "Non-expected (neutral)"
  
  # make stats in the labels of the plots 
  
  # only 21 (for A) and 12 (for B) have negative phylogenetic signal in the traits
  # remove them 
  
  results_all_failed <- results_all_failed[which(results_all_failed$Expected_signal_A!="Expected (negative)"),]
  results_all_failed <- results_all_failed[which(results_all_failed$Expected_signal_B!="Expected (negative)"),]
  
  # at least one with expected signal:
  results_all_failed$Expected_signal_A[which(results_all_failed$Expected_signal_B=="Expected")] <- "Expected"
  
  
  # no size
  results_signal_all <- results_all_failed %>% 
    group_by(Expected_signal_A,Inference) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  ggplot(results_signal_all,aes(x=Expected_signal_A, y = count))+  geom_bar( position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=colors)+ theme(rect = element_rect(fill = "transparent"))
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_",eco_matrix,"_inference_new_expected_traits.pdf"), device=cairo_pdf, plot = last_plot(), width = 4.5, height = 5.5,  bg = "transparent")
  
  
  
  # Correlation d (signal interactions) and Mantel corr (traits)
  
  results_all_failed$Inference_interactions <- "Not significant signal"
  results_all_failed$Inference_interactions[which(results_all_failed$Ratio_MSE>0)] <- "Significant signal"
  
  ggplot(results_all_failed,aes(x=mantel_cor_traits_A, y = dA))+  geom_point(alpha=0.7, aes(shape=Expected_signal_A,color=Inference_interactions))+
    xlab("Mantel correlation in traits")+ylab("d value of PBLM") +labs(title=" ")+scale_fill_manual(values=colors) + 
    theme_classic()+ geom_hline(yintercept = 0)+ ylim(c(0,1))+
    scale_color_manual(values=c("#f7dc6f", "#196f3d"))+ theme(rect = element_rect(fill = "transparent"))
  
  ggsave(filename=paste0("plot_correlation_d_PBLM_interactions_and_Mantel_traits_",eco_matrix,"_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 6.5, height = 4.5,   bg = "transparent")
  
  
  ggplot(results_all_failed,aes(x=mantel_cor_traits_B, y = dB))+  geom_point(alpha=0.7, aes(shape=Expected_signal_B,color=Inference_interactions))+
    xlab("Mantel correlation in traits")+ylab("d value of PBLM") +labs(title=" ")+scale_fill_manual(values=colors) + 
    theme_classic()+ geom_hline(yintercept = 0)+ ylim(c(0,1))+
    scale_color_manual(values=c("#f7dc6f", "#196f3d"))+ theme(rect = element_rect(fill = "transparent"))
  
  ggsave(filename=paste0("plot_correlation_d_PBLM_interactions_and_Mantel_traits_",eco_matrix,"_clade_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 6.5, height = 4.5,   bg = "transparent")
  
  
}



#####  Step 3-B': Plot outputs phylogenetic signal - simul 1 - PBLM with MSE with traits Pagel lambda  ######

# on cluster

rm(list=ls())

library(ggplot2)

setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1")


transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))


eco_matrix="PBLM"
eco_matrix="PBLM_binary"

for (eco_matrix in c( "PBLM", "PBLM_binary")){
  
  info = list.files(path = "results_PBLM/", pattern = paste0('_',eco_matrix,'.csv'))
  info <- info[grep(x=info,pattern = "seed")]
  
  results_all <- c()
  for (i in info){
    
    if (length(grep(x=i, pattern = "_neutral_")==1)>0) {param="neutral"}
    if (length(grep(x=i, pattern = "_mutualism_1_")==1)>0) {param="mutualism_1"}
    if (length(grep(x=i, pattern = "_mutualism_2_")==1)>0) {param="mutualism_2"}
    if (length(grep(x=i, pattern = "_mutualism_3_")==1)>0) {param="mutualism_3"}
    if (length(grep(x=i, pattern = "_mutualism_4_")==1)>0) {param="mutualism_4"}
    if (length(grep(x=i, pattern = "_mutualism_5_")==1)>0) {param="mutualism_5"}
    if (length(grep(x=i, pattern = "_mutualism_6_")==1)>0) {param="mutualism_6"}
    if (length(grep(x=i, pattern = "_antagonism_1_")==1)>0) {param="antagonism_1"}
    if (length(grep(x=i, pattern = "_antagonism_2_")==1)>0) {param="antagonism_2"}
    if (length(grep(x=i, pattern = "_antagonism_3_")==1)>0) {param="antagonism_3"}
    if (length(grep(x=i, pattern = "_antagonism_4_")==1)>0) {param="antagonism_4"}
    if (length(grep(x=i, pattern = "_antagonism_5_")==1)>0) {param="antagonism_5"}
    if (length(grep(x=i, pattern = "_antagonism_6_")==1)>0) {param="antagonism_6"}
    if (length(grep(x=i, pattern = "_antagonism_7_")==1)>0) {param="antagonism_7"}
    if (length(grep(x=i, pattern = "_antagonism_8_")==1)>0) {param="antagonism_8"}
    if (length(grep(x=i, pattern = "_antagonism_9_")==1)>0) {param="antagonism_9"}
    
    if (length(grep(x=i, pattern = "_A_")==1)>0) {size="3000"}
    if (length(grep(x=i, pattern = "_B_")==1)>0) {size="4000"}
    if (length(grep(x=i, pattern = "_C_")==1)>0) {size="5000"}
    if (length(grep(x=i, pattern = "_D_")==1)>0) {size="2000"}
    if (length(grep(x=i, pattern = "_E_")==1)>0) {size="1000"}
    if (length(grep(x=i, pattern = "_F_")==1)>0) {size="500"}
    
    results_seed <- read.table(paste0("results_PBLM/",i), header=T, sep=";")
    
    name_simul <- as.character(unlist(gsub("_PBLM.csv", "", gsub("_PBLM_binary.csv", "", gsub("results_phylogenetic_signal_", "", i)))))
    
    # pb with D and E
    if (ncol(results_seed)==1) {results_seed = t(results_seed)}
    
    results_all <- rbind(results_all, c(name_simul, param, size, results_seed))
  }
  
  results_all <- data.frame(results_all)
  
  colnames(results_all) <- c("name", "param","size", "nb_A", "nb_B", "dA", "dB", "MSETotal", "MSEFull", "MSEStar", "MSEBase")
  
  results_all$param <- as.character(results_all$param)
  results_all$size <- as.character(results_all$size)
  results_all$nb_A <- as.numeric(as.character(results_all$nb_A))
  results_all$nb_B <- as.numeric(as.character(results_all$nb_B))
  results_all$dA <- as.numeric(as.character(results_all$dA))
  results_all$dB <- as.numeric(as.character(results_all$dB))
  results_all$MSETotal <- as.numeric(as.character(results_all$MSETotal))
  results_all$MSEFull <- as.numeric(as.character(results_all$MSEFull))
  results_all$MSEStar <- as.numeric(as.character(results_all$MSEStar))
  results_all$MSEBase <- as.numeric(as.character(results_all$MSEBase))
  
  results_all$size <- as.factor(results_all$size)
  results_all$size <-  factor(results_all$size,levels(results_all$size)[c(5,1,2,3,4,6)])
  
  
  results_all$Diff_MSE <- results_all$MSEStar - results_all$MSEFull
  
  results_all$Ratio_MSE <- (results_all$MSEStar - results_all$MSEFull)/results_all$MSEStar
  
  results_all$name <- unlist(results_all$name)
  
  write.table(results_all,paste0("results_phylogenetic_signal_simul_1_",eco_matrix,".csv"), quote=F, sep=";",row.names=F)
  
  
  results_all_failed <- results_all
  results_all <- results_all[which(results_all$dA!="inference failed"),]
  
  
  results_all_failed$Inference <- "Inference failed"
  results_all_failed$Inference[which(results_all_failed$Ratio_MSE<=0)] <- "a) Not significant signal"
  results_all_failed$Inference[which(results_all_failed$Ratio_MSE>0)] <- "b) Significant signal"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dA>0.05))] <- "c) Significant signal (d>0.05)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dB>0.05))] <- "c) Significant signal (d>0.05)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dA>0.15))] <- "d) Significant signal (d>0.15)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dB>0.15))] <- "d) Significant signal (d>0.15)"
  
  colors <- c()
  if ("a) Not significant signal" %in% results_all_failed$Inference) {colors <- c(colors, "#f7dc6f")}
  if ("b) Significant signal" %in% results_all_failed$Inference) {colors <- c(colors, "#7dcea0")}
  if ("c) Significant signal (d>0.05)" %in% results_all_failed$Inference) {colors <- c(colors, "#27ae60")}
  if ("d) Significant signal (d>0.15)" %in% results_all_failed$Inference) {colors <- c(colors, "#196f3d")}
  if ("e) Significant anti-signal (d>1)" %in% results_all_failed$Inference) {colors <- c(colors, "#943126")}
  
  
  library(dplyr)
  brks <- c(0, 0.25, 0.5, 0.75, 1)
  
  results_all_failed$size_tot <- "<150"
  results_all_failed$size_tot[which(results_all_failed$nb_A+results_all_failed$nb_B>=150)] <- "150-250"
  results_all_failed$size_tot[which(results_all_failed$nb_A+results_all_failed$nb_B>250)] <- ">250"
  
  colors <- c()
  if ("a) Not significant signal" %in% results_all_failed$Inference) {colors <- c(colors, "#f7dc6f")}
  if ("b) Significant signal" %in% results_all_failed$Inference) {colors <- c(colors, "#7dcea0")}
  if ("c) Significant signal (d>0.05)" %in% results_all_failed$Inference) {colors <- c(colors, "#27ae60")}
  if ("d) Significant signal (d>0.15)" %in% results_all_failed$Inference) {colors <- c(colors, "#196f3d")}
  if ("e) Significant anti-signal (d>1)" %in% results_all_failed$Inference) {colors <- c(colors, "#943126")}
  
  
  ##### Add signal traits 
  results_traits <- read.table(paste0("results_phylogenetic_signal_simul_1_in_traits_Pagel_lambda.csv"), sep=";", header=T)
  
  results_traits$param <- as.character(results_traits$param)
  results_traits$lambda_A <- as.numeric(as.character(results_traits$lambda_A))
  results_traits$pvalue_high_A <- as.numeric(as.character(results_traits$p_LRT_A))
  results_traits$lambda_B <- as.numeric(as.character(results_traits$lambda_B))
  results_traits$pvalue_high_B <- as.numeric(as.character(results_traits$p_LRT_B))
  results_traits$nb_A <- as.numeric(as.character(results_traits$nb_A))
  results_traits$nb_B <- as.numeric(as.character(results_traits$nb_B))
  results_traits$name_simul_type <- paste0(results_traits$name, "_", results_traits$param, "_seed_", results_traits$seed )
  
  results_traits$Inference_traits_B <- "Not significant signal"
  results_traits$Inference_traits_B[which(results_traits$pvalue_high_B<0.05)] <- "Significant signal"
  
  results_traits$Inference_traits_A <- "Not significant signal"
  results_traits$Inference_traits_A[which(results_traits$pvalue_high_A<0.05)] <- "Significant signal"
  
  
  # make correspondance
  
  results_traits$name_simul_type <- as.character(results_traits$name_simul_type)
  results_all_failed$name <- as.character(results_all_failed$name)
  
  for (i in 1:nrow(results_all_failed)){
    index <- which(results_traits$name_simul_type==results_all_failed$name[i])
    results_all_failed$Inference_traits_A[i] <- results_traits$Inference_traits_A[index]
    results_all_failed$Inference_traits_B[i] <- results_traits$Inference_traits_B[index]
    results_all_failed$lambda_A[i] <- results_traits$lambda_A[index]
    results_all_failed$lambda_B[i] <- results_traits$lambda_B[index]
  }
  
  # expected signal
  results_all_failed$Expected_signal_A <- "Expected"
  results_all_failed$Expected_signal_A[which(results_all_failed$Inference_traits_A=="Not significant signal")] <- "Non-expected (traits)"
  results_all_failed$Expected_signal_A[grep("neutral", results_all_failed$name)] <- "Non-expected (neutral)"
  
  results_all_failed$Expected_signal_B <- "Expected"
  results_all_failed$Expected_signal_B[which(results_all_failed$Inference_traits_B=="Not significant signal")] <- "Non-expected (traits)"
  results_all_failed$Expected_signal_B[grep("neutral", results_all_failed$name)] <- "Non-expected (neutral)"
  
  # at least one with expected signal:
  results_all_failed$Expected_signal_A[which(results_all_failed$Expected_signal_B=="Expected")] <- "Expected"
  
  
  # no size
  results_signal_all <- results_all_failed %>% 
    group_by(Expected_signal_A,Inference) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  ggplot(results_signal_all,aes(x=Expected_signal_A, y = count))+  geom_bar( position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=colors)+ theme(rect = element_rect(fill = "transparent"))
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_",eco_matrix,"_inference_new_expected_traits_Pagel_lambda.pdf"), device=cairo_pdf, plot = last_plot(), width = 4.5, height = 5.5,  bg = "transparent")
  
  
  
  # Correlation d (signal interactions) and lambda (traits)
  
  results_all_failed$Inference_interactions <- "Not significant signal"
  results_all_failed$Inference_interactions[which(results_all_failed$Ratio_MSE>0)] <- "Significant signal"
  
  ggplot(results_all_failed,aes(x=lambda_A, y = dA))+  geom_point(alpha=0.7, aes(shape=Expected_signal_A,color=Inference_interactions))+
    xlab("Pagel's lambda in traits")+ylab("d value of PBLM") +labs(title=" ")+scale_fill_manual(values=colors) + 
    theme_classic()+ geom_hline(yintercept = 0)+ ylim(c(0,1))+
    scale_color_manual(values=c("#f7dc6f", "#196f3d"))+ theme(rect = element_rect(fill = "transparent"))
  
  ggsave(filename=paste0("plot_correlation_d_PBLM_interactions_and_Pagel_lambda_traits_",eco_matrix,"_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 6.5, height = 4.5,   bg = "transparent")
  
  
  ggplot(results_all_failed,aes(x=lambda_B, y = dB))+  geom_point(alpha=0.7, aes(shape=Expected_signal_B,color=Inference_interactions))+
    xlab("Pagel's lambda in traits")+ylab("d value of PBLM") +labs(title=" ")+scale_fill_manual(values=colors) + 
    theme_classic()+ geom_hline(yintercept = 0)+ ylim(c(0,1))+
    scale_color_manual(values=c("#f7dc6f", "#196f3d"))+ theme(rect = element_rect(fill = "transparent"))
  
  ggsave(filename=paste0("plot_correlation_d_PBLM_interactions_and_Pagel_lambda_traits_",eco_matrix,"_clade_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 6.5, height = 4.5,   bg = "transparent")
  
}



#####  Step 3-B: Plot outputs phylogenetic signal - simul 1 - PBLM with bootstrap ######

# on cluster

rm(list=ls())

library(ggplot2)

setwd("/Users/bperez/ownCloud/Recherche/These/ISYEB/signal_phylo/simul_1/")
setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1")


transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))


eco_matrix="PBLM"
eco_matrix="PBLM_binary"

for (eco_matrix in c("PBLM","PBLM_binary")){
  
  info = list.files(path = "./", pattern = "results_phylogenetic_signal_PBLM_bootstrap_")
  info <- info[grep(x=info,pattern = "seed")]
  info <- info[grep(x=info,pattern = paste0(eco_matrix,".csv"))]
  
  results_all <- c()
  for (i in info){
    
    if (length(grep(x=i, pattern = "_neutral_")==1)>0) {param="neutral"}
    if (length(grep(x=i, pattern = "_mutualism_1_")==1)>0) {param="mutualism_1"}
    if (length(grep(x=i, pattern = "_mutualism_2_")==1)>0) {param="mutualism_2"}
    if (length(grep(x=i, pattern = "_mutualism_3_")==1)>0) {param="mutualism_3"}
    if (length(grep(x=i, pattern = "_mutualism_4_")==1)>0) {param="mutualism_4"}
    if (length(grep(x=i, pattern = "_mutualism_5_")==1)>0) {param="mutualism_5"}
    if (length(grep(x=i, pattern = "_mutualism_6_")==1)>0) {param="mutualism_6"}
    if (length(grep(x=i, pattern = "_antagonism_1_")==1)>0) {param="antagonism_1"}
    if (length(grep(x=i, pattern = "_antagonism_2_")==1)>0) {param="antagonism_2"}
    if (length(grep(x=i, pattern = "_antagonism_3_")==1)>0) {param="antagonism_3"}
    if (length(grep(x=i, pattern = "_antagonism_4_")==1)>0) {param="antagonism_4"}
    if (length(grep(x=i, pattern = "_antagonism_5_")==1)>0) {param="antagonism_5"}
    if (length(grep(x=i, pattern = "_antagonism_6_")==1)>0) {param="antagonism_6"}
    if (length(grep(x=i, pattern = "_antagonism_7_")==1)>0) {param="antagonism_7"}
    if (length(grep(x=i, pattern = "_antagonism_8_")==1)>0) {param="antagonism_8"}
    if (length(grep(x=i, pattern = "_antagonism_9_")==1)>0) {param="antagonism_9"}
    
    if (length(grep(x=i, pattern = "_A_")==1)>0) {size="3000"}
    if (length(grep(x=i, pattern = "_B_")==1)>0) {size="4000"}
    if (length(grep(x=i, pattern = "_C_")==1)>0) {size="5000"}
    if (length(grep(x=i, pattern = "_D_")==1)>0) {size="2000"}
    if (length(grep(x=i, pattern = "_E_")==1)>0) {size="1000"}
    if (length(grep(x=i, pattern = "_F_")==1)>0) {size="500"}
    
    results_seed <- read.table(paste0(i), header=T, sep=";")
    
    
    # pb with D and E
    if (ncol(results_seed)==1) {results_seed = t(results_seed)}
    
    #if (results_seed[4]!="inference failed"){
    results_all <- rbind(results_all, c(param, size, results_seed))
    #}
  }
  
  results_all <- data.frame(results_all)
  colnames(results_all) <- c("param","size", "nb_A", "nb_B", "dA", "dB", "MSETotal", "MSEFull", "MSEStar", "MSEBase", "dA_min", "dA_max", "dB_min", "dB_max")
  
  results_all$param <- as.character(results_all$param)
  results_all$size <- as.character(results_all$size)
  results_all$nb_A <- as.numeric(as.character(results_all$nb_A))
  results_all$nb_B <- as.numeric(as.character(results_all$nb_B))
  results_all$dA <- as.numeric(as.character(results_all$dA))
  results_all$dB <- as.numeric(as.character(results_all$dB))
  results_all$MSETotal <- as.numeric(as.character(results_all$MSETotal))
  results_all$MSEFull <- as.numeric(as.character(results_all$MSEFull))
  results_all$MSEStar <- as.numeric(as.character(results_all$MSEStar))
  results_all$MSEBase <- as.numeric(as.character(results_all$MSEBase))
  results_all$dA_min <- as.numeric(as.character(results_all$dA_min))
  results_all$dB_min <- as.numeric(as.character(results_all$dB_min))
  results_all$dA_max <- as.numeric(as.character(results_all$dA_max))
  results_all$dB_max <- as.numeric(as.character(results_all$dB_max))
  
  results_all$size <- as.factor(results_all$size)
  
  
  results_all_failed <- results_all
  results_all <- results_all[which(results_all$dA!="inference failed"),]
  
  # Plot 
  library(dplyr)
  
  results_all_failed$Inference <- "Inference failed"
  results_all_failed$Inference[intersect(which(results_all_failed$dA_min<=0.05),which(results_all_failed$dB_min<=0.05))] <- "a) Not significant signal"
  results_all_failed$Inference[intersect(which(results_all_failed$dA_min>0.05),which(results_all_failed$dB_min<=0.05))] <- "b) Significant signal in guild A only"
  results_all_failed$Inference[intersect(which(results_all_failed$dB_min>0.05),which(results_all_failed$dA_min<=0.05))] <- "c) Significant signal in guild B only"
  results_all_failed$Inference[intersect(which(results_all_failed$dA_min>0.05),which(results_all_failed$dB_min>0.05))] <- "d) Significant signal in both guilds"
  
  
  colors <- c()
  if ("a) Not significant signal" %in% results_all_failed$Inference) {colors <- c(colors, "#f7dc6f")}
  if ("b) Significant signal in guild A only" %in% results_all_failed$Inference) {colors <- c(colors, "#7dcea0")}
  if ("c) Significant signal in guild B only" %in% results_all_failed$Inference) {colors <- c(colors, "#27ae60")}
  if ("d) Significant signal in both guilds" %in% results_all_failed$Inference) {colors <- c(colors, "#196f3d")}
  
  
  results_all_failed$size_tot <- "simul_F"
  
  results_signal_all <- results_all_failed %>% 
    group_by(param,Inference) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  brks <- c(0, 0.25, 0.5, 0.75, 1)
  ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=colors)+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_",eco_matrix,"_bootstrap.pdf"), device=cairo_pdf, plot = last_plot(), width = 10, height = 7)
  
  print(table(results_all_failed$param))
  
  
}


#####  Step 3-C: Plot connectance newtork - simul 1  ######


rm(list=ls())

library(ggplot2)

setwd("/Users/bperez/ownCloud/Recherche/These/ISYEB/signal_phylo/simul_1/")


transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))


info = list.files(pattern = paste0('.csv'))
info <- info[grep(x=info,pattern = "network")]
info <- info[grep(x=info,pattern = "seed")]


results_all <- c()
for (i in info){
  
  if (length(grep(x=i, pattern = "_neutral_")==1)>0) {param="neutral"}
  if (length(grep(x=i, pattern = "_mutualism_1_")==1)>0) {param="mutualism_1"}
  if (length(grep(x=i, pattern = "_mutualism_2_")==1)>0) {param="mutualism_2"}
  if (length(grep(x=i, pattern = "_mutualism_3_")==1)>0) {param="mutualism_3"}
  if (length(grep(x=i, pattern = "_mutualism_4_")==1)>0) {param="mutualism_4"}
  if (length(grep(x=i, pattern = "_mutualism_5_")==1)>0) {param="mutualism_5"}
  if (length(grep(x=i, pattern = "_mutualism_6_")==1)>0) {param="mutualism_6"}
  if (length(grep(x=i, pattern = "_antagonism_1_")==1)>0) {param="antagonism_1"}
  if (length(grep(x=i, pattern = "_antagonism_2_")==1)>0) {param="antagonism_2"}
  if (length(grep(x=i, pattern = "_antagonism_3_")==1)>0) {param="antagonism_3"}
  if (length(grep(x=i, pattern = "_antagonism_4_")==1)>0) {param="antagonism_4"}
  if (length(grep(x=i, pattern = "_antagonism_5_")==1)>0) {param="antagonism_5"}
  if (length(grep(x=i, pattern = "_antagonism_6_")==1)>0) {param="antagonism_6"}
  if (length(grep(x=i, pattern = "_antagonism_7_")==1)>0) {param="antagonism_7"}
  if (length(grep(x=i, pattern = "_antagonism_8_")==1)>0) {param="antagonism_8"}
  if (length(grep(x=i, pattern = "_antagonism_9_")==1)>0) {param="antagonism_9"}
  
  if (length(grep(x=i, pattern = "_A_")==1)>0) {size="3000"}
  if (length(grep(x=i, pattern = "_B_")==1)>0) {size="4000"}
  if (length(grep(x=i, pattern = "_C_")==1)>0) {size="5000"}
  if (length(grep(x=i, pattern = "_D_")==1)>0) {size="2000"}
  if (length(grep(x=i, pattern = "_E_")==1)>0) {size="1000"}
  if (length(grep(x=i, pattern = "_F_")==1)>0) {size="500"}
  
  results_seed <- read.table(i, header=T, sep=";")
  
  results_all <- rbind(results_all, c(param, size, ncol(results_seed), nrow(results_seed), length(which(results_seed>0))/(nrow(results_seed)*ncol(results_seed))  ))
}


results_all <- data.frame(results_all)

colnames(results_all) <- c("param","size", "nb_A", "nb_B", "connectance")


results_all$param <- as.character(results_all$param)

results_all$size <- as.character(results_all$size)
results_all$nb_A <- as.numeric(as.character(results_all$nb_A))
results_all$nb_B <- as.numeric(as.character(results_all$nb_B))
results_all$connectance <- as.numeric(as.character(results_all$connectance))

results_all$size <- as.factor(results_all$size)
results_all$size <-  factor(results_all$size,levels(results_all$size)[c(5,1,2,3,4,6)])


# Plot 
ggplot(results_all,aes(x=param,y=connectance,fill=size))+xlab("Parameters")+ylab("Connectance") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
  geom_hline(yintercept=0, linetype="dashed", color="#273746")+
  scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                              bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
  geom_boxplot(alpha=0.9,size=0.75, color="#273746")

ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_connectance.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)

results_all$nb_tot <- results_all$nb_A + results_all$nb_B

results_all$size_tot <- "<150"
results_all$size_tot[results_all$nb_tot>=150] <- "150-250"
results_all$size_tot[results_all$nb_tot>250] <- ">250"


ggplot(results_all,aes(x=param,y=connectance,fill=size_tot))+xlab("Parameters")+ylab("Connectance") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
  geom_hline(yintercept=0, linetype="dashed", color="#273746")+
  scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                              bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
  geom_boxplot(alpha=0.9,size=0.75, color="#273746")

ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_connectance_size.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)




########  Plot connectance assymetry (A - 10% most abundant)

# on cluster

rm(list=ls())

library(ggplot2)

setwd("/Users/bperez/ownCloud/Recherche/These/ISYEB/signal_phylo/simul_1/")


transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))


info = list.files(pattern = paste0('.csv'))
info <- info[grep(x=info,pattern = "network")]
info <- info[grep(x=info,pattern = "seed")]


results_all <- c()
for (i in info){
  
  if (length(grep(x=i, pattern = "_neutral_")==1)>0) {param="neutral"}
  if (length(grep(x=i, pattern = "_mutualism_1_")==1)>0) {param="mutualism_1"}
  if (length(grep(x=i, pattern = "_mutualism_2_")==1)>0) {param="mutualism_2"}
  if (length(grep(x=i, pattern = "_mutualism_3_")==1)>0) {param="mutualism_3"}
  if (length(grep(x=i, pattern = "_mutualism_4_")==1)>0) {param="mutualism_4"}
  if (length(grep(x=i, pattern = "_mutualism_5_")==1)>0) {param="mutualism_5"}
  if (length(grep(x=i, pattern = "_mutualism_6_")==1)>0) {param="mutualism_6"}
  if (length(grep(x=i, pattern = "_antagonism_1_")==1)>0) {param="antagonism_1"}
  if (length(grep(x=i, pattern = "_antagonism_2_")==1)>0) {param="antagonism_2"}
  if (length(grep(x=i, pattern = "_antagonism_3_")==1)>0) {param="antagonism_3"}
  if (length(grep(x=i, pattern = "_antagonism_4_")==1)>0) {param="antagonism_4"}
  if (length(grep(x=i, pattern = "_antagonism_5_")==1)>0) {param="antagonism_5"}
  if (length(grep(x=i, pattern = "_antagonism_6_")==1)>0) {param="antagonism_6"}
  if (length(grep(x=i, pattern = "_antagonism_7_")==1)>0) {param="antagonism_7"}
  if (length(grep(x=i, pattern = "_antagonism_8_")==1)>0) {param="antagonism_8"}
  if (length(grep(x=i, pattern = "_antagonism_9_")==1)>0) {param="antagonism_9"}
  
  if (length(grep(x=i, pattern = "_A_")==1)>0) {size="3000"}
  if (length(grep(x=i, pattern = "_B_")==1)>0) {size="4000"}
  if (length(grep(x=i, pattern = "_C_")==1)>0) {size="5000"}
  if (length(grep(x=i, pattern = "_D_")==1)>0) {size="2000"}
  if (length(grep(x=i, pattern = "_E_")==1)>0) {size="1000"}
  if (length(grep(x=i, pattern = "_F_")==1)>0) {size="500"}
  
  network <- read.table(i, header=T, sep=";")
  
  selected_A <- order(colSums(network),decreasing = T)[1:min(max(10,floor(ncol(network)/10)),ncol(network))]
  
  network <- network[,selected_A]
  network <- network[which(rowSums(network)>0),] # remove species from B without any interaction anymore
  
  results_all <- rbind(results_all, c(param, size, ncol(network), nrow(network), length(which(network>0))/(nrow(network)*ncol(network))  ))
}


results_all <- data.frame(results_all)

colnames(results_all) <- c("param","size", "nb_A", "nb_B", "connectance")


results_all$param <- as.character(results_all$param)

results_all$size <- as.character(results_all$size)
results_all$nb_A <- as.numeric(as.character(results_all$nb_A))
results_all$nb_B <- as.numeric(as.character(results_all$nb_B))
results_all$connectance <- as.numeric(as.character(results_all$connectance))

results_all$size <- as.factor(results_all$size)
results_all$size <-  factor(results_all$size,levels(results_all$size)[c(5,1,2,3,4,6)])

# Plot 
ggplot(results_all,aes(x=param,y=connectance,fill=size))+xlab("Parameters")+ylab("Connectance") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
  geom_hline(yintercept=0, linetype="dashed", color="#273746")+
  scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                              bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
  geom_boxplot(alpha=0.9,size=0.75, color="#273746")

ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_assymetry_A_connectance.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)

results_all$nb_tot <- results_all$nb_A + results_all$nb_B

results_all$size_tot <- "<80"
results_all$size_tot[results_all$nb_tot>=80] <- "80-140"
results_all$size_tot[results_all$nb_tot>140] <- ">140"


ggplot(results_all,aes(x=param,y=connectance,fill=size_tot))+xlab("Parameters")+ylab("Connectance") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
  geom_hline(yintercept=0, linetype="dashed", color="#273746")+
  scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                              bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
  geom_boxplot(alpha=0.9,size=0.75, color="#273746")

ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_assymetry_A_connectance_size.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)





#####  Step 3-D: Final plot phylogenetic signal - simul 1 - Mantel tests only - type I error / power  ######


rm(list=ls())

library(ggplot2)

setwd("/Users/bperez/ownCloud/Recherche/These/ISYEB/signal_phylo/simul_1/")

transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))


eco_matrix="Jaccard_weighted"
eco_matrix="Jaccard_binary"
eco_matrix="GUniFrac"
eco_matrix="UniFrac_unweighted"

method="Pearson"

results_performance <- c()

for (eco_matrix in c("Jaccard_weighted","Jaccard_binary","GUniFrac","UniFrac_unweighted")){
  
  for (method in c("Pearson", "Spearman", "Kendall")){
    
    print(eco_matrix)
    print(method)
    
    results_signal <- read.table(paste0("results_phylogenetic_signal_simul_1_",eco_matrix,"_", method, ".csv"), header=T, sep=";")
    
    results_signal$mantel_cor_A <- as.numeric(results_signal$mantel_cor_A)
    results_signal$mantel_cor_B <- as.numeric(results_signal$mantel_cor_B)
    results_signal$pvalue_high_A <- as.numeric(results_signal$pvalue_high_A)
    results_signal$pvalue_high_B <- as.numeric(results_signal$pvalue_high_B)
    
    # plot size A and B 
    results_signal$size <- as.character(results_signal$size)
    results_signal$size[results_signal$size=="A"] <- 3000
    results_signal$size[results_signal$size=="B"] <- 4000
    results_signal$size[results_signal$size=="C"] <- 5000
    results_signal$size[results_signal$size=="D"] <- 2000
    results_signal$size[results_signal$size=="E"] <- 1000
    results_signal$size[grep(pattern = "_F_", results_signal$name)] <- 500
    
    
    results_signal$size <- as.factor(results_signal$size)
    results_signal$size <-  factor(results_signal$size,levels(results_signal$size)[c(5,1,2,3,4,6)])
    
    results_signal$name_simul <- substr(results_signal$name, 1, 9)
    
    results_signal$name_simul_type <- paste0(results_signal$name_simul, "_", results_signal$param)
    
    results_signal$seed <- sapply(strsplit(split="seed_", as.character(results_signal$name)),"[[", 2)
    
    ##### signal traits 
    
    results_traits <- read.table(paste0("results_phylogenetic_signal_simul_1_in_traits.csv"), sep=";", header=T)
    
    results_traits$param <- as.character(results_traits$param)
    results_traits$mantel_cor_A <- as.numeric(as.character(results_traits$mantel_cor_A))
    results_traits$pvalue_high_A <- as.numeric(as.character(results_traits$pvalue_high_A))
    results_traits$mantel_cor_B <- as.numeric(as.character(results_traits$mantel_cor_B))
    results_traits$pvalue_high_B <- as.numeric(as.character(results_traits$pvalue_high_B))
    results_traits$nb_A <- as.numeric(as.character(results_traits$nb_A))
    results_traits$nb_B <- as.numeric(as.character(results_traits$nb_B))
    results_traits$name_simul_type <- paste0(results_traits$name, "_", results_traits$param, "_seed_", results_traits$seed )
    
    results_traits$Inference_traits_B <- "Not significant signal"
    results_traits$Inference_traits_B[which(results_traits$pvalue_high_B<0.05)] <- "Significant signal"
    results_traits$Inference_traits_B[which(results_traits$pvalue_low_B<0.05)] <- "Significant anti-signal"
    
    results_traits$Inference_traits_A <- "Not significant signal"
    results_traits$Inference_traits_A[which(results_traits$pvalue_high_A<0.05)] <- "Significant signal"
    results_traits$Inference_traits_A[which(results_traits$pvalue_low_A<0.05)] <- "Significant anti-signal"
    
    
    
    
    # make correspondance
    
    results_traits$name_simul_type <- as.character(results_traits$name_simul_type)
    results_signal$name <- as.character(results_signal$name)
    
    for (i in 1:nrow(results_signal)){
      index <- which(results_traits$name_simul_type==results_signal$name[i])
      results_signal$Inference_traits_A[i] <- results_traits$Inference_traits_A[index]
      results_signal$Inference_traits_B[i] <- results_traits$Inference_traits_B[index]
      results_signal$mantel_cor_traits_A[i] <- results_traits$mantel_cor_A[index]
      results_signal$mantel_cor_traits_B[i] <- results_traits$mantel_cor_B[index]
    }
    
    # expected signal
    results_signal$Expected_signal_A <- "Expected"
    results_signal$Expected_signal_A[which(results_signal$Inference_traits_A=="Significant anti-signal")] <- "Expected (negative)"
    results_signal$Expected_signal_A[which(results_signal$Inference_traits_A=="Not significant signal")] <- "Non-expected (traits)"
    results_signal$Expected_signal_A[grep("neutral", results_signal$name)] <- "Non-expected (neutral)"
    
    results_signal$Expected_signal_B <- "Expected"
    results_signal$Expected_signal_B[which(results_signal$Inference_traits_B=="Significant anti-signal")] <- "Expected (negative)"
    results_signal$Expected_signal_B[which(results_signal$Inference_traits_B=="Not significant signal")] <- "Non-expected (traits)"
    results_signal$Expected_signal_B[grep("neutral", results_signal$name)] <- "Non-expected (neutral)"
    
    
    results_signal <- results_signal[which(results_signal$Expected_signal_A!="Expected (negative)"),]
    results_signal <- results_signal[which(results_signal$Expected_signal_B!="Expected (negative)"),]
    
    
    ######  plots
    
    library(dplyr)
    brks <- c(0, 0.25, 0.5, 0.75, 1)
    
    
    ##  Resume for clade A
    
    results_signal$Inference <- "d) Not significant signal"
    results_signal$Inference[which(results_signal$pvalue_high_A<0.05)] <- "e) Significant signal"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.05))] <- "f) Significant signal (R>0.05)"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.15))] <- "g) Significant signal (R>0.15)"
    results_signal$Inference[which(results_signal$pvalue_low_A<0.05)] <- "c) Significant anti-signal"
    results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.05)))] <- "b) Significant anti-signal (R<-0.05)"
    results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.15)))] <- "a) Significant anti-signal (R<-0.15)"
    
    colors <- c()
    if ("a) Significant anti-signal (R<-0.15)" %in% results_signal$Inference) {colors <- c(colors, "#943126")}
    if ("b) Significant anti-signal (R<-0.05)" %in% results_signal$Inference) {colors <- c(colors, "#e74c3c")}
    if ("c) Significant anti-signal" %in% results_signal$Inference) {colors <- c(colors, "#f1948a")}
    if ("d) Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
    if ("e) Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
    if ("f) Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
    if ("g) Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}
    
    
    # no size
    results_signal_all <- results_signal %>% 
      group_by(Expected_signal_A,Inference) %>% 
      summarise(count=n()) %>% 
      mutate(perc=count/sum(count))
    ggplot(results_signal_all,aes(x=Expected_signal_A, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
      xlab(" ")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
      scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
      scale_color_manual(values=colors)+ theme(rect = element_rect(fill = "transparent"))
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Mantel_tests_",eco_matrix,"_", method,"_inference_with_traits_final_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 5, height = 5.5,  bg = "transparent")
    
    
    
    # Correlation R (signal interactions) and R (traits)
    
    results_signal$Inference_interactions_A <- "Not significant signal"
    results_signal$Inference_interactions_A[which(results_signal$pvalue_high_A<0.05)] <- "Significant signal"
    results_signal$Inference_interactions_A[which(results_signal$pvalue_low_A<0.05)] <- "Significant anti-signal"
    
    
    ggplot(results_signal,aes(x=mantel_cor_traits_A, y = mantel_cor_A))+  geom_point(alpha=0.7, aes(shape=Expected_signal_A,color=Inference_interactions_A))+
      xlab("Mantel correlation in traits")+ylab("Mantel correlation in species interactions") +labs(title=" ")+scale_fill_manual(values=colors) + 
      theme_classic()+ geom_hline(yintercept = 0)+
      scale_color_manual(values=c("#f7dc6f", "#943126", "#196f3d"))+ theme(rect = element_rect(fill = "transparent"))
    
    ggsave(filename=paste0("plot_correlation_Mantel_interactions_and_Mantel_traits_",eco_matrix,"_", method,"_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 6.5, height = 4.5,  bg = "transparent")
    
    
    
    
    
    ##  Resume for clade B
    
    results_signal$Inference <- "d) Not significant signal"
    results_signal$Inference[which(results_signal$pvalue_high_B<0.05)] <- "e) Significant signal"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.05))] <- "f) Significant signal (R>0.05)"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.15))] <- "g) Significant signal (R>0.15)"
    results_signal$Inference[which(results_signal$pvalue_low_B<0.05)] <- "c) Significant anti-signal"
    results_signal$Inference[intersect(which(results_signal$pvalue_low_B<0.05), which(results_signal$mantel_cor_B<(-0.05)))] <- "b) Significant anti-signal (R<-0.05)"
    results_signal$Inference[intersect(which(results_signal$pvalue_low_B<0.05), which(results_signal$mantel_cor_B<(-0.15)))] <- "a) Significant anti-signal (R<-0.15)"
    
    colors <- c()
    if ("a) Significant anti-signal (R<-0.15)" %in% results_signal$Inference) {colors <- c(colors, "#943126")}
    if ("b) Significant anti-signal (R<-0.05)" %in% results_signal$Inference) {colors <- c(colors, "#e74c3c")}
    if ("c) Significant anti-signal" %in% results_signal$Inference) {colors <- c(colors, "#f1948a")}
    if ("d) Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
    if ("e) Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
    if ("f) Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
    if ("g) Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}
    
    
    # no size
    results_signal_all <- results_signal %>% 
      group_by(Expected_signal_B,Inference) %>% 
      summarise(count=n()) %>% 
      mutate(perc=count/sum(count))
    ggplot(results_signal_all,aes(x=Expected_signal_B, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
      xlab(" ")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
      scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
      scale_color_manual(values=colors) + theme(rect = element_rect(fill = "transparent"))
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Mantel_tests_",eco_matrix,"_", method,"_inference_with_traits_final_clade_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 5, height = 5.5,  bg = "transparent")
    
    
    
    # Correlation R (signal interactions) and R (traits)
    
    results_signal$Inference_interactions_B <- "Not significant signal"
    results_signal$Inference_interactions_B[which(results_signal$pvalue_high_B<0.05)] <- "Significant signal"
    results_signal$Inference_interactions_B[which(results_signal$pvalue_low_B<0.05)] <- "Significant anti-signal"
    
    
    ggplot(results_signal,aes(x=mantel_cor_traits_B, y = mantel_cor_B))+  geom_point(alpha=0.7, aes(shape=Expected_signal_B,color=Inference_interactions_B))+
      xlab("Mantel correlation in traits")+ylab("Mantel correlation in species interactions") +labs(title=" ")+scale_fill_manual(values=colors) + 
      theme_classic()+ geom_hline(yintercept = 0)+
      scale_color_manual(values=c("#f7dc6f", "#943126", "#196f3d"))+ theme(rect = element_rect(fill = "transparent"))
    
    ggsave(filename=paste0("plot_correlation_Mantel_interactions_and_Mantel_traits_",eco_matrix,"_", method,"_clade_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 6.5, height = 4.5,  bg = "transparent")
    
    
    
    # type-I error
    
    print("type 1 error")
    results_signal$size_tot <- "<150"
    results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>=150)] <- "150-250"
    results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>250)] <- ">250"
    
    
    print("A")
    pA1 <- print(length(intersect(which(results_signal$pvalue_high_A<0.05),intersect(which(results_signal$Expected_signal_A %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot=="<150"))))/length(intersect(which(results_signal$Expected_signal_A %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot=="<150"))))
    pA2 <- print(length(intersect(which(results_signal$pvalue_high_A<0.05),intersect(which(results_signal$Expected_signal_A %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot=="150-250"))))/length(intersect(which(results_signal$Expected_signal_A %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot=="150-250"))))
    pA3 <- print(length(intersect(which(results_signal$pvalue_high_A<0.05),intersect(which(results_signal$Expected_signal_A %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot==">250"))))/length(intersect(which(results_signal$Expected_signal_A %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot==">250"))))
    
    print("B")
    pB1 <- print(length(intersect(which(results_signal$pvalue_high_B<0.05),intersect(which(results_signal$Expected_signal_B %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot=="<150"))))/length(intersect(which(results_signal$Expected_signal_B %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot=="<150"))))
    pB2 <- print(length(intersect(which(results_signal$pvalue_high_B<0.05),intersect(which(results_signal$Expected_signal_B %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot=="150-250"))))/length(intersect(which(results_signal$Expected_signal_B %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot=="150-250"))))
    pB3 <- print(length(intersect(which(results_signal$pvalue_high_B<0.05),intersect(which(results_signal$Expected_signal_B %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot==">250"))))/length(intersect(which(results_signal$Expected_signal_B %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot==">250"))))
    
    results_performance <- rbind(results_performance, c(eco_matrix, method, pA1, pB1))
    results_performance <- rbind(results_performance, c(eco_matrix, method, pA2, pB2))
    results_performance <- rbind(results_performance, c(eco_matrix, method, pA3, pB3))
    
    
    
    # Power 
    
    if (method=="Pearson"){
      
      print("Power")
      
      results_all <- read.table(paste0("results_correlation_degree_phylo_simul_1.csv"), header=T,sep=";")
      
      results_all$param <- as.character(results_all$param)
      results_all$size <- as.character(results_all$size)
      results_all$mantel_cor_A <- as.numeric(as.character(results_all$mantel_cor_A))
      results_all$pvalue_high_A <- as.numeric(as.character(results_all$pvalue_high_A))
      results_all$mantel_cor_B <- as.numeric(as.character(results_all$mantel_cor_B))
      results_all$pvalue_high_B <- as.numeric(as.character(results_all$pvalue_high_B))
      results_all$nb_A <- as.numeric(as.character(results_all$nb_A))
      results_all$nb_B <- as.numeric(as.character(results_all$nb_B))
      
      results_all <- results_all[which(results_all$name %in%results_signal$name),]
      print(all(results_all$name==results_signal$name))
      
      # Partial Mantel test
      
      results_partial <- read.table(paste0("results_partial_Mantel_test_degree_simul_1_",eco_matrix,".csv"), header=T, sep=";")
      results_partial$mantel_cor_A <- as.numeric(results_partial$mantel_cor_A)
      results_partial$mantel_cor_B <- as.numeric(results_partial$mantel_cor_B)
      results_partial$pvalue_high_A <- as.numeric(results_partial$pvalue_high_A)
      results_partial$pvalue_high_B <- as.numeric(results_partial$pvalue_high_B)
      
      results_partial <- results_partial[which(results_partial$name %in%results_signal$name),]
      print(all(results_partial$name==results_signal$name))
      
      
      print("Power regular Mantel test:")
      print(length(which(results_signal$pvalue_high_A[which(results_signal$Expected_signal_A=="Expected")]<0.05))/length(which(results_signal$Expected_signal_A=="Expected")))
      print(length(which(results_signal$pvalue_high_B[which(results_signal$Expected_signal_B=="Expected")]<0.05))/length(which(results_signal$Expected_signal_B=="Expected")))
      
      print("Power partial Mantel test:")
      print(length(which(results_partial$pvalue_high_A[which(results_signal$Expected_signal_A=="Expected")]<0.05))/length(which(results_signal$Expected_signal_A=="Expected")))
      print(length(which(results_partial$pvalue_high_B[which(results_signal$Expected_signal_B=="Expected")]<0.05))/length(which(results_signal$Expected_signal_B=="Expected")))
      
      print("Power separate Mantel tests:")
      print(length(intersect(which(results_signal$pvalue_high_A[which(results_signal$Expected_signal_A=="Expected")]<0.05), which(results_all$pvalue_high_A[which(results_signal$Expected_signal_A=="Expected")]>0.05)))/length(which(results_signal$Expected_signal_A=="Expected")))
      print(length(intersect(which(results_signal$pvalue_high_B[which(results_signal$Expected_signal_B=="Expected")]<0.05), which(results_all$pvalue_high_B[which(results_signal$Expected_signal_B=="Expected")]>0.05)))/length(which(results_signal$Expected_signal_B=="Expected")))
      
      
    }
    
  }
}


write.table(results_performance,"results_typeI_error_simul_1.csv", sep=";", row.names=F, quote=F, dec=",")




#####  Step 3-D': Final plot phylogenetic signal - simul 1 - Pagel lambda + Mantel tests  - type I error / power  ######


rm(list=ls())

library(ggplot2)

setwd("/Users/bperez/ownCloud/Recherche/These/ISYEB/signal_phylo/simul_1/")

transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))


eco_matrix="Jaccard_weighted"
eco_matrix="Jaccard_binary"
eco_matrix="GUniFrac"
eco_matrix="UniFrac_unweighted"

method="Pearson"

results_performance <- c()

for (eco_matrix in c("Jaccard_weighted","Jaccard_binary","GUniFrac","UniFrac_unweighted")){
  
  for (method in c("Pearson", "Spearman", "Kendall")){
    
    print(eco_matrix)
    print(method)
    
    results_signal <- read.table(paste0("results_phylogenetic_signal_simul_1_",eco_matrix,"_", method, ".csv"), header=T, sep=";")
    
    results_signal$mantel_cor_A <- as.numeric(results_signal$mantel_cor_A)
    results_signal$mantel_cor_B <- as.numeric(results_signal$mantel_cor_B)
    results_signal$pvalue_high_A <- as.numeric(results_signal$pvalue_high_A)
    results_signal$pvalue_high_B <- as.numeric(results_signal$pvalue_high_B)
    
    # plot size A and B 
    results_signal$size <- as.character(results_signal$size)
    results_signal$size[results_signal$size=="A"] <- 3000
    results_signal$size[results_signal$size=="B"] <- 4000
    results_signal$size[results_signal$size=="C"] <- 5000
    results_signal$size[results_signal$size=="D"] <- 2000
    results_signal$size[results_signal$size=="E"] <- 1000
    results_signal$size[grep(pattern = "_F_", results_signal$name)] <- 500
    
    
    results_signal$size <- as.factor(results_signal$size)
    results_signal$size <-  factor(results_signal$size,levels(results_signal$size)[c(5,1,2,3,4,6)])
    
    results_signal$name_simul <- substr(results_signal$name, 1, 9)
    
    results_signal$name_simul_type <- paste0(results_signal$name_simul, "_", results_signal$param)
    
    results_signal$seed <- sapply(strsplit(split="seed_", as.character(results_signal$name)),"[[", 2)
    
    ##### signal traits 
    
    results_traits <- read.table(paste0("results_phylogenetic_signal_simul_1_in_traits_Pagel_lambda.csv"), sep=";", header=T)
    
    results_traits$param <- as.character(results_traits$param)
    results_traits$lambda_A <- as.numeric(as.character(results_traits$lambda_A))
    results_traits$pvalue_high_A <- as.numeric(as.character(results_traits$p_LRT_A))
    results_traits$lambda_B <- as.numeric(as.character(results_traits$lambda_B))
    results_traits$pvalue_high_B <- as.numeric(as.character(results_traits$p_LRT_B))
    results_traits$nb_A <- as.numeric(as.character(results_traits$nb_A))
    results_traits$nb_B <- as.numeric(as.character(results_traits$nb_B))
    results_traits$name_simul_type <- paste0(results_traits$name, "_", results_traits$param, "_seed_", results_traits$seed )
    
    results_traits$Inference_traits_B <- "Not significant signal"
    results_traits$Inference_traits_B[which(results_traits$pvalue_high_B<0.05)] <- "Significant signal"
    
    results_traits$Inference_traits_A <- "Not significant signal"
    results_traits$Inference_traits_A[which(results_traits$pvalue_high_A<0.05)] <- "Significant signal"
    
    
    
    # make correspondance
    
    results_traits$name_simul_type <- as.character(results_traits$name_simul_type)
    results_signal$name <- as.character(results_signal$name)
    
    for (i in 1:nrow(results_signal)){
      index <- which(results_traits$name_simul_type==results_signal$name[i])
      results_signal$Inference_traits_A[i] <- results_traits$Inference_traits_A[index]
      results_signal$Inference_traits_B[i] <- results_traits$Inference_traits_B[index]
      results_signal$lambda_A[i] <- results_traits$lambda_A[index]
      results_signal$lambda_B[i] <- results_traits$lambda_B[index]
    }
    
    # expected signal
    results_signal$Expected_signal_A <- "Expected"
    results_signal$Expected_signal_A[which(results_signal$Inference_traits_A=="Not significant signal")] <- "Non-expected (traits)"
    results_signal$Expected_signal_A[grep("neutral", results_signal$name)] <- "Non-expected (neutral)"
    
    results_signal$Expected_signal_B <- "Expected"
    results_signal$Expected_signal_B[which(results_signal$Inference_traits_B=="Not significant signal")] <- "Non-expected (traits)"
    results_signal$Expected_signal_B[grep("neutral", results_signal$name)] <- "Non-expected (neutral)"
    
    # make stats in the labels of the plots 
    
    
    ######  plots
    
    library(dplyr)
    brks <- c(0, 0.25, 0.5, 0.75, 1)
    
    
    ##  Resume for clade A
    
    results_signal$Inference <- "d) Not significant signal"
    results_signal$Inference[which(results_signal$pvalue_high_A<0.05)] <- "e) Significant signal"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.05))] <- "f) Significant signal (R>0.05)"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.15))] <- "g) Significant signal (R>0.15)"
    results_signal$Inference[which(results_signal$pvalue_low_A<0.05)] <- "c) Significant anti-signal"
    results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.05)))] <- "b) Significant anti-signal (R<-0.05)"
    results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.15)))] <- "a) Significant anti-signal (R<-0.15)"
    
    colors <- c()
    if ("a) Significant anti-signal (R<-0.15)" %in% results_signal$Inference) {colors <- c(colors, "#943126")}
    if ("b) Significant anti-signal (R<-0.05)" %in% results_signal$Inference) {colors <- c(colors, "#e74c3c")}
    if ("c) Significant anti-signal" %in% results_signal$Inference) {colors <- c(colors, "#f1948a")}
    if ("d) Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
    if ("e) Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
    if ("f) Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
    if ("g) Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}
    
    
    # no size
    results_signal_all <- results_signal %>% 
      group_by(Expected_signal_A,Inference) %>% 
      summarise(count=n()) %>% 
      mutate(perc=count/sum(count))
    ggplot(results_signal_all,aes(x=Expected_signal_A, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
      xlab(" ")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
      scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
      scale_color_manual(values=colors)+ theme(rect = element_rect(fill = "transparent"))
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Mantel_tests_",eco_matrix,"_", method,"_inference_with_traits_final_Pagel_lambda_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 5, height = 5.5,  bg = "transparent")
    
    
    # Correlation R (signal interactions) and lambda (traits)
    
    results_signal$Inference_interactions_A <- "Not significant signal"
    results_signal$Inference_interactions_A[which(results_signal$pvalue_high_A<0.05)] <- "Significant signal"
    results_signal$Inference_interactions_A[which(results_signal$pvalue_low_A<0.05)] <- "Significant anti-signal"
    
    
    ggplot(results_signal,aes(x=lambda_A, y = mantel_cor_A))+  geom_point(alpha=0.7, aes(shape=Expected_signal_A,color=Inference_interactions_A))+
      xlab("Pagel's lambda in traits")+ylab("Mantel correlation in species interactions") +labs(title=" ")+scale_fill_manual(values=colors) + 
      theme_classic()+ geom_hline(yintercept = 0)+
      scale_color_manual(values=c("#f7dc6f", "#943126", "#196f3d"))+ theme(rect = element_rect(fill = "transparent"))
    
    ggsave(filename=paste0("plot_correlation_Mantel_interactions_and_Pagel_lambda_traits_",eco_matrix,"_", method,"_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 6.5, height = 4.5,   bg = "transparent")
    
    
    
    
    ##  Resume for clade B
    
    results_signal$Inference <- "d) Not significant signal"
    results_signal$Inference[which(results_signal$pvalue_high_B<0.05)] <- "e) Significant signal"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.05))] <- "f) Significant signal (R>0.05)"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.15))] <- "g) Significant signal (R>0.15)"
    results_signal$Inference[which(results_signal$pvalue_low_B<0.05)] <- "c) Significant anti-signal"
    results_signal$Inference[intersect(which(results_signal$pvalue_low_B<0.05), which(results_signal$mantel_cor_B<(-0.05)))] <- "b) Significant anti-signal (R<-0.05)"
    results_signal$Inference[intersect(which(results_signal$pvalue_low_B<0.05), which(results_signal$mantel_cor_B<(-0.15)))] <- "a) Significant anti-signal (R<-0.15)"
    
    colors <- c()
    if ("a) Significant anti-signal (R<-0.15)" %in% results_signal$Inference) {colors <- c(colors, "#943126")}
    if ("b) Significant anti-signal (R<-0.05)" %in% results_signal$Inference) {colors <- c(colors, "#e74c3c")}
    if ("c) Significant anti-signal" %in% results_signal$Inference) {colors <- c(colors, "#f1948a")}
    if ("d) Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
    if ("e) Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
    if ("f) Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
    if ("g) Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}
    
    
    # no size
    results_signal_all <- results_signal %>% 
      group_by(Expected_signal_B,Inference) %>% 
      summarise(count=n()) %>% 
      mutate(perc=count/sum(count))
    ggplot(results_signal_all,aes(x=Expected_signal_B, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
      xlab(" ")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
      scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
      scale_color_manual(values=colors) + theme(rect = element_rect(fill = "transparent"))
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Mantel_tests_",eco_matrix,"_", method,"_inference_with_traits_final_Pagel_lambda_clade_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 5, height = 5.5,  bg = "transparent")
    
    
    # Correlation R (signal interactions) and lambda (traits)
    
    results_signal$Inference_interactions_B <- "Not significant signal"
    results_signal$Inference_interactions_B[which(results_signal$pvalue_high_B<0.05)] <- "Significant signal"
    results_signal$Inference_interactions_B[which(results_signal$pvalue_low_B<0.05)] <- "Significant anti-signal"
    
    
    ggplot(results_signal,aes(x=lambda_B, y = mantel_cor_B))+  geom_point(alpha=0.7, aes(shape=Expected_signal_B,color=Inference_interactions_B))+
      xlab("Pagel's lambda in traits")+ylab("Mantel correlation in species interactions") +labs(title=" ")+scale_fill_manual(values=colors) + 
      theme_classic()+ geom_hline(yintercept = 0)+
      scale_color_manual(values=c("#f7dc6f", "#943126", "#196f3d"))+ theme(rect = element_rect(fill = "transparent"))
    
    ggsave(filename=paste0("plot_correlation_Mantel_interactions_and_Pagel_lambda_traits_",eco_matrix,"_", method,"_clade_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 6.5, height = 4.5,   bg = "transparent")
    
    
    
    
    # type-I error
    
    print("type 1 error")
    results_signal$size_tot <- "<150"
    results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>=150)] <- "150-250"
    results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>250)] <- ">250"
    
    
    print("A")
    pA1 <- print(length(intersect(which(results_signal$pvalue_high_A<0.05),intersect(which(results_signal$Expected_signal_A %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot=="<150"))))/length(intersect(which(results_signal$Expected_signal_A %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot=="<150"))))
    pA2 <- print(length(intersect(which(results_signal$pvalue_high_A<0.05),intersect(which(results_signal$Expected_signal_A %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot=="150-250"))))/length(intersect(which(results_signal$Expected_signal_A %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot=="150-250"))))
    pA3 <- print(length(intersect(which(results_signal$pvalue_high_A<0.05),intersect(which(results_signal$Expected_signal_A %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot==">250"))))/length(intersect(which(results_signal$Expected_signal_A %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot==">250"))))
    
    print("B")
    pB1 <- print(length(intersect(which(results_signal$pvalue_high_B<0.05),intersect(which(results_signal$Expected_signal_B %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot=="<150"))))/length(intersect(which(results_signal$Expected_signal_B %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot=="<150"))))
    pB2 <- print(length(intersect(which(results_signal$pvalue_high_B<0.05),intersect(which(results_signal$Expected_signal_B %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot=="150-250"))))/length(intersect(which(results_signal$Expected_signal_B %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot=="150-250"))))
    pB3 <- print(length(intersect(which(results_signal$pvalue_high_B<0.05),intersect(which(results_signal$Expected_signal_B %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot==">250"))))/length(intersect(which(results_signal$Expected_signal_B %in% c("Non-expected (neutral)", "Non-expected (traits)")), which(results_signal$size_tot==">250"))))
    
    results_performance <- rbind(results_performance, c(eco_matrix, method, pA1, pB1))
    results_performance <- rbind(results_performance, c(eco_matrix, method, pA2, pB2))
    results_performance <- rbind(results_performance, c(eco_matrix, method, pA3, pB3))
    
    
    
    # Power 
    
    if (method=="Pearson"){
      
      print("Power")
      
      results_all <- read.table(paste0("results_correlation_degree_phylo_simul_1.csv"), header=T,sep=";")
      
      results_all$param <- as.character(results_all$param)
      results_all$size <- as.character(results_all$size)
      results_all$mantel_cor_A <- as.numeric(as.character(results_all$mantel_cor_A))
      results_all$pvalue_high_A <- as.numeric(as.character(results_all$pvalue_high_A))
      results_all$mantel_cor_B <- as.numeric(as.character(results_all$mantel_cor_B))
      results_all$pvalue_high_B <- as.numeric(as.character(results_all$pvalue_high_B))
      results_all$nb_A <- as.numeric(as.character(results_all$nb_A))
      results_all$nb_B <- as.numeric(as.character(results_all$nb_B))
      
      results_all <- results_all[which(results_all$name %in%results_signal$name),]
      print(all(results_all$name==results_signal$name))
      
      # Partial Mantel test
      
      results_partial <- read.table(paste0("results_partial_Mantel_test_degree_simul_1_",eco_matrix,".csv"), header=T, sep=";")
      results_partial$mantel_cor_A <- as.numeric(results_partial$mantel_cor_A)
      results_partial$mantel_cor_B <- as.numeric(results_partial$mantel_cor_B)
      results_partial$pvalue_high_A <- as.numeric(results_partial$pvalue_high_A)
      results_partial$pvalue_high_B <- as.numeric(results_partial$pvalue_high_B)
      
      results_partial <- results_partial[which(results_partial$name %in%results_signal$name),]
      print(all(results_partial$name==results_signal$name))
      
      
      print("Power regular Mantel test:")
      print(length(which(results_signal$pvalue_high_A[which(results_signal$Expected_signal_A=="Expected")]<0.05))/length(which(results_signal$Expected_signal_A=="Expected")))
      print(length(which(results_signal$pvalue_high_B[which(results_signal$Expected_signal_B=="Expected")]<0.05))/length(which(results_signal$Expected_signal_B=="Expected")))
      
      print("Power partial Mantel test:")
      print(length(which(results_partial$pvalue_high_A[which(results_signal$Expected_signal_A=="Expected")]<0.05))/length(which(results_signal$Expected_signal_A=="Expected")))
      print(length(which(results_partial$pvalue_high_B[which(results_signal$Expected_signal_B=="Expected")]<0.05))/length(which(results_signal$Expected_signal_B=="Expected")))
      
      print("Power separate Mantel tests:")
      print(length(intersect(which(results_signal$pvalue_high_A[which(results_signal$Expected_signal_A=="Expected")]<0.05), which(results_all$pvalue_high_A[which(results_signal$Expected_signal_A=="Expected")]>0.05)))/length(which(results_signal$Expected_signal_A=="Expected")))
      print(length(intersect(which(results_signal$pvalue_high_B[which(results_signal$Expected_signal_B=="Expected")]<0.05), which(results_all$pvalue_high_B[which(results_signal$Expected_signal_B=="Expected")]>0.05)))/length(which(results_signal$Expected_signal_B=="Expected")))
      
      
    }
    
  }
}


write.table(results_performance,"results_typeI_error_simul_1_Pagel_lambda.csv", sep=";", row.names=F, quote=F, dec=",")


results_traits_pagel <- read.table(paste0("results_phylogenetic_signal_simul_1_in_traits_Pagel_lambda.csv"), sep=";", header=T)
results_traits_mantel <- read.table(paste0("results_phylogenetic_signal_simul_1_in_traits.csv"), sep=";", header=T)


results_traits_pagel$param <- as.character(results_traits_pagel$param)
results_traits_pagel$pvalue_high_A <- as.numeric(as.character(results_traits_pagel$p_LRT_A))
results_traits_pagel$pvalue_high_B <- as.numeric(as.character(results_traits_pagel$p_LRT_B))
results_traits_pagel$name_simul_type <- paste0(results_traits_pagel$name, "_", results_traits_pagel$param, "_seed_", results_traits_pagel$seed )

results_traits_pagel$Inference_traits_B <- "Not significant signal (Pagel)"
results_traits_pagel$Inference_traits_B[which(results_traits_pagel$pvalue_high_B<0.05)] <- "Significant signal (Pagel)"

results_traits_pagel$Inference_traits_A <- "Not significant signal (Pagel)"
results_traits_pagel$Inference_traits_A[which(results_traits_pagel$pvalue_high_A<0.05)] <- "Significant signal (Pagel)"

results_traits_pagel$pvalue_high_mantel_A <- as.numeric(as.character(results_traits_mantel$pvalue_high_A))
results_traits_pagel$pvalue_low_mantel_A <- as.numeric(as.character(results_traits_mantel$pvalue_low_A))
results_traits_pagel$pvalue_high_mantel_B <- as.numeric(as.character(results_traits_mantel$pvalue_high_B))
results_traits_pagel$pvalue_low_mantel_B <- as.numeric(as.character(results_traits_mantel$pvalue_low_B))


results_traits_pagel$Inference_traits_mantel_B <- "Not significant signal (Mantel)"
results_traits_pagel$Inference_traits_mantel_B[which(results_traits_pagel$pvalue_high_mantel_B<0.05)] <- "Significant signal (Mantel)"

results_traits_pagel$Inference_traits_mantel_A <- "Not significant signal (Mantel)"
results_traits_pagel$Inference_traits_mantel_A[which(results_traits_pagel$pvalue_high_mantel_A<0.05)] <- "Significant signal (Mantel)"


table(results_traits_pagel$Inference_traits_mantel_A, results_traits_pagel$Inference_traits_A)
table(results_traits_pagel$Inference_traits_mantel_B, results_traits_pagel$Inference_traits_B)
(2056+242)/2400
(2114+197)/2400

results_traits_mut <- results_traits_pagel[grep("mutualism",results_traits_pagel$param),]

table(results_traits_mut$Inference_traits_mantel_A, results_traits_mut$Inference_traits_A)
384/(28+384) # 93%
(384+69)/720
(384+28)/720
(239+28)/720
(239+69)/720

table(results_traits_mut$Inference_traits_mantel_B, results_traits_mut$Inference_traits_B)
450/(17+450) # 93%
(450+62)/720
(450+17)/720
(191+17)/720
(191+62)/720

#####  Step 4-A: Compute phylogenetic signal with assymetry (10%) - simul 1 ######


# on cluster
rm(list=ls())

library(vegan)
library(phytools)
library(ggplot2)
library(GUniFrac,lib.loc="/users/biodiv/bperez/packages/")
library(parallel)
require(R.utils)

setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1")

transparent_theme <- theme(
  panel.grid = element_blank(),
  axis.line = element_line("black"),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))


source("../script/function_phylo_signal_network.R")

dyn.load("../script/permute.so")

name="simul_1"

list_networks <- as.vector(outer(c("simul_1_A","simul_1_B","simul_1_C","simul_1_D","simul_1_E", "simul_1_F"), c(paste0("neutral_seed_",1:100), paste0("mutualism_", as.vector(outer(1:6, 1:20, paste, sep="_seed_"))), paste0("antagonism_",as.vector(outer(1:9, 1:20, paste, sep="_seed_")))), paste, sep="_"))

compute_phylo_signal <- function(name, method, correlation, nperm){
  
  if (!file.exists(paste0("results_PBLM_assymetry/results_phylogenetic_signal_assymetry_A_",name,"_",method,".csv"))){
    
    network <- read.table(paste0("network_",name,".csv"), header=T,sep=";")
    
    colnames(network) <- seq(1:ncol(network))
    rownames(network) <- seq(1:nrow(network))
    
    tree_A <- read.tree(paste0("network_tree_guild_A_",name,".tre"))
    tree_B <- read.tree(paste0("network_tree_guild_B_",name,".tre"))
    
    # Select 10% of the most abundant species in clade A (A in columns and B in rows) - or minimum 10
    
    selected_A <- order(colSums(network),decreasing = T)[1:min(max(10,floor(ncol(network)/10)),ncol(network))]
    
    network <- network[,selected_A]
    network <- network[which(rowSums(network)>0),] # remove species from B without any interaction anymore
    
    tree_A <- drop.tip(tree_A, tip=tree_A$tip.label[which(!tree_A$tip.label %in% colnames(network))])
    tree_B <- drop.tip(tree_B, tip=tree_B$tip.label[which(!tree_B$tip.label %in% rownames(network))])
    
    
    # check that there are several species
    if ((nrow(network)>1)&(ncol(network)!=1)){
      
      # Not all species interactins with each others
      if (length(which(network==0))>0){
        
        output <- phylosignal_network(network, tree_A, tree_B, method = method, nperm = nperm, correlation = correlation)
        
        if (method %in% c("PBLM", "PBLM_binary")){
          names(output) <- c("nb_A", "nb_B", "dA", "dB", "MSETotal", "MSEFull", "MSEStar", "MSEBase")
          write.table(output,paste0("results_PBLM_assymetry/results_phylogenetic_signal_assymetry_A_",name,"_",method,".csv"), quote=F, sep=";",row.names=F)}
        
        return(output)
        
      }else{
        return(c(nb_A, nb_B, 0,0,0,0,0,0))}
    }else{
      return(c(NA,NA, NA,NA,NA,NA,NA,NA))}
  }else{
    return(c(NA,NA,NA,NA,NA,NA,NA,NA))}
}


correlation="Pearson"

method="PBLM"
method="PBLM_binary"

method="Jaccard_weighted"
method="Jaccard_binary"
method="GUniFrac"
method="UniFrac_unweighted"

if (method %in% c("Jaccard_weighted", "Jaccard_binary", "GUniFrac", "UniFrac_unweighted")){
  
  for (correlation in c("Pearson", "Spearman", "Kendall")){
    
    nperm=10000
    if (correlation=="Kendall") {nperm=100}
    
    results_signal <- matrix(unlist(mclapply(list_networks,compute_phylo_signal,mc.cores=20, method=method, mc.preschedule=F, correlation=correlation, nperm=nperm)),nrow=length(list_networks),byrow=T)
    results_signal <- cbind(list_networks, results_signal)
    results_signal <- data.frame(results_signal, stringsAsFactors =F)
    colnames(results_signal) <- c("name","nb_A", "nb_B", "mantel_cor_A","pvalue_high_A","pvalue_low_A", "mantel_cor_B", "pvalue_high_B", "pvalue_low_B")
    
    results_signal$size <- "A"
    results_signal$size[grep(pattern = "_B_", results_signal$name)] <- "B"
    results_signal$size[grep(pattern = "_C_", results_signal$name)] <- "C"
    results_signal$size[grep(pattern = "_D_", results_signal$name)] <- "D"
    results_signal$size[grep(pattern = "_E_", results_signal$name)] <- "E"
    results_signal$size[grep(pattern = "_F_", results_signal$name)] <- "F"
    results_signal$param <- "neutral"
    results_signal$param[grep(pattern = "_antagonism_1_", results_signal$name)] <- "antagonism_1"
    results_signal$param[grep(pattern = "_antagonism_2_", results_signal$name)] <- "antagonism_2"
    results_signal$param[grep(pattern = "_antagonism_3_", results_signal$name)] <- "antagonism_3"
    results_signal$param[grep(pattern = "_antagonism_4_", results_signal$name)] <- "antagonism_4"
    results_signal$param[grep(pattern = "_antagonism_5_", results_signal$name)] <- "antagonism_5"
    results_signal$param[grep(pattern = "_antagonism_6_", results_signal$name)] <- "antagonism_6"
    results_signal$param[grep(pattern = "_antagonism_7_", results_signal$name)] <- "antagonism_7"
    results_signal$param[grep(pattern = "_antagonism_8_", results_signal$name)] <- "antagonism_8"
    results_signal$param[grep(pattern = "_antagonism_9_", results_signal$name)] <- "antagonism_9"
    results_signal$param[grep(pattern = "_mutualism_1_", results_signal$name)] <- "mutualism_1"
    results_signal$param[grep(pattern = "_mutualism_2_", results_signal$name)] <- "mutualism_2"
    results_signal$param[grep(pattern = "_mutualism_3_", results_signal$name)] <- "mutualism_3"
    results_signal$param[grep(pattern = "_mutualism_4_", results_signal$name)] <- "mutualism_4"
    results_signal$param[grep(pattern = "_mutualism_5_", results_signal$name)] <- "mutualism_5"
    results_signal$param[grep(pattern = "_mutualism_6_", results_signal$name)] <- "mutualism_6"
    
    write.table(results_signal,paste0("results_phylogenetic_signal_simul_1_assymetry_A_",method,"_", correlation,".csv"), quote=F, sep=";",row.names=F)
    
  }
}


if (method %in% c("PBLM_binary", "PBLM")){
  
  nperm=10000 # useless
  correlation="Kendall" # useless
  
  results_signal <- matrix(unlist(mclapply(list_networks,compute_phylo_signal,mc.cores=20, mc.preschedule=F, method=method, correlation=correlation, nperm=nperm)),nrow=length(list_networks),byrow=T)
  results_signal <- cbind(list_networks, results_signal)
  results_signal <- data.frame(results_signal, stringsAsFactors =F)
  colnames(results_signal) <- c("name","nb_A", "nb_B", "dA", "dB", "MSETotal", "MSEFull", "MSEStar", "MSEBase")
  
  results_signal$size <- "A"
  results_signal$size[grep(pattern = "_B_", results_signal$name)] <- "B"
  results_signal$size[grep(pattern = "_C_", results_signal$name)] <- "C"
  results_signal$size[grep(pattern = "_D_", results_signal$name)] <- "D"
  results_signal$size[grep(pattern = "_E_", results_signal$name)] <- "E"
  results_signal$size[grep(pattern = "_F_", results_signal$name)] <- "F"
  results_signal$param <- "neutral"
  results_signal$param[grep(pattern = "_antagonism_1_", results_signal$name)] <- "antagonism_1"
  results_signal$param[grep(pattern = "_antagonism_2_", results_signal$name)] <- "antagonism_2"
  results_signal$param[grep(pattern = "_antagonism_3_", results_signal$name)] <- "antagonism_3"
  results_signal$param[grep(pattern = "_antagonism_4_", results_signal$name)] <- "antagonism_4"
  results_signal$param[grep(pattern = "_antagonism_5_", results_signal$name)] <- "antagonism_5"
  results_signal$param[grep(pattern = "_antagonism_6_", results_signal$name)] <- "antagonism_6"
  results_signal$param[grep(pattern = "_antagonism_7_", results_signal$name)] <- "antagonism_7"
  results_signal$param[grep(pattern = "_antagonism_8_", results_signal$name)] <- "antagonism_8"
  results_signal$param[grep(pattern = "_antagonism_9_", results_signal$name)] <- "antagonism_9"
  results_signal$param[grep(pattern = "_mutualism_1_", results_signal$name)] <- "mutualism_1"
  results_signal$param[grep(pattern = "_mutualism_2_", results_signal$name)] <- "mutualism_2"
  results_signal$param[grep(pattern = "_mutualism_3_", results_signal$name)] <- "mutualism_3"
  results_signal$param[grep(pattern = "_mutualism_4_", results_signal$name)] <- "mutualism_4"
  results_signal$param[grep(pattern = "_mutualism_5_", results_signal$name)] <- "mutualism_5"
  results_signal$param[grep(pattern = "_mutualism_6_", results_signal$name)] <- "mutualism_6"
  
  write.table(results_signal,paste0("results_phylogenetic_signal_simul_1_assymetry_A_",method,".csv"), quote=F, sep=";",row.names=F)
  
}


#####  Step 4-B: Plot outputs phylogenetic signal with assymetry (10%) - simul 1 - Mantel tests ######


rm(list=ls())

library(ggplot2)

setwd("/Users/bperez/ownCloud/Recherche/These/ISYEB/signal_phylo/simul_1/")

transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))


eco_matrix="Jaccard_weighted"
eco_matrix="Jaccard_binary"
eco_matrix="GUniFrac"
eco_matrix="UniFrac_unweighted"

method="Pearson"
method="Kendall"

for (eco_matrix in c("Jaccard_weighted","Jaccard_binary","GUniFrac","UniFrac_unweighted")){
  
  for (method in c("Pearson", "Spearman", "Kendall")){
    
    results_signal <- read.table(paste0("results_phylogenetic_signal_simul_1_assymetry_A_",eco_matrix,"_", method, ".csv"), header=T, sep=";")
    
    results_signal$mantel_cor_A <- as.numeric(results_signal$mantel_cor_A)
    results_signal$mantel_cor_B <- as.numeric(results_signal$mantel_cor_B)
    results_signal$pvalue_high_A <- as.numeric(results_signal$pvalue_high_A)
    results_signal$pvalue_high_B <- as.numeric(results_signal$pvalue_high_B)
    
    # plot size A and B 
    results_signal$size <- as.character(results_signal$size)
    results_signal$size[results_signal$size=="A"] <- 3000
    results_signal$size[results_signal$size=="B"] <- 4000
    results_signal$size[results_signal$size=="C"] <- 5000
    results_signal$size[results_signal$size=="D"] <- 2000
    results_signal$size[results_signal$size=="E"] <- 1000
    results_signal$size[grep(pattern = "_F_", results_signal$name)] <- 500
    
    
    results_signal$size <- as.factor(results_signal$size)
    results_signal$size <-  factor(results_signal$size,levels(results_signal$size)[c(5,1,2,3,4,6)])
    
    
    ggplot(results_signal,aes(x=param,y=nb_A,fill=size))+xlab("Parameters")+ylab("Number of species (A)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
      ylim(0, max(results_signal$nb_A))+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      geom_boxplot(alpha=0.9,size=0.75, color="#273746")
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_assymetry_A_size_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 10, height = 5)
    
    ggplot(results_signal,aes(x=param,y=nb_B,fill=size))+xlab("Parameters")+ylab("Number of species (B)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
      ylim(0, max(results_signal$nb_B))+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      geom_boxplot(alpha=0.9,size=0.75, color="#273746")
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_assymetry_A_size_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 10, height = 5)
    
    
    
    ##  Resume for clade A
    
    results_signal$Inference <- "Not significant signal"
    results_signal$Inference[which(results_signal$pvalue_high_A<0.05)] <- "Significant signal"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.05))] <- "Significant signal (R>0.05)"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.15))] <- "Significant signal (R>0.15)"
    
    colors <- c()
    if ("Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
    if ("Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
    if ("Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
    if ("Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}
    
    library(dplyr)
    results_signal_all <- results_signal %>% 
      group_by(param,Inference,size) %>% 
      summarise(count=n()) %>% 
      mutate(perc=count/sum(count))
    
    brks <- c(0, 0.25, 0.5, 0.75, 1)
    ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar( position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
      xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
      scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
      scale_color_manual(values=colors)+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      facet_grid(~size)
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_assymetry_A_Mantel_tests_",eco_matrix,"_", method,"_inference_size_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 14, height = 7)
    
    
    results_signal$size_tot <- "<80"
    results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>=75)] <- "80-140"
    results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>140)] <- ">140"
    
    results_signal_all <- results_signal %>% 
      group_by(param,Inference,size_tot) %>% 
      summarise(count=n()) %>% 
      mutate(perc=count/sum(count))
    
    results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
    results_signal_all$size_tot <-  factor(results_signal_all$size_tot,levels(results_signal_all$size_tot)[c(1,3,2)])
    
    
    ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
      xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
      scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
      scale_color_manual(values=colors)+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      facet_grid(~size_tot)
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_assymetry_A_Mantel_tests_",eco_matrix,"_", method,"_inference_size_tot_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 10, height = 7)
    
    
    ### Add anticorrelation
    
    results_signal$Inference <- "d) Not significant signal"
    results_signal$Inference[which(results_signal$pvalue_high_A<0.05)] <- "e) Significant signal"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.05))] <- "f) Significant signal (R>0.05)"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.15))] <- "g) Significant signal (R>0.15)"
    results_signal$Inference[which(results_signal$pvalue_low_A<0.05)] <- "c) Significant anti-signal"
    results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.05)))] <- "b) Significant anti-signal (R<-0.05)"
    results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.15)))] <- "a) Significant anti-signal (R<-0.15)"
    
    colors <- c()
    if ("a) Significant anti-signal (R<-0.15)" %in% results_signal$Inference) {colors <- c(colors, "#943126")}
    if ("b) Significant anti-signal (R<-0.05)" %in% results_signal$Inference) {colors <- c(colors, "#e74c3c")}
    if ("c) Significant anti-signal" %in% results_signal$Inference) {colors <- c(colors, "#f1948a")}
    if ("d) Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
    if ("e) Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
    if ("f) Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
    if ("g) Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}
    
    results_signal_all <- results_signal %>% 
      group_by(param,Inference,size_tot) %>% 
      summarise(count=n()) %>% 
      mutate(perc=count/sum(count))
    
    results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
    results_signal_all$size_tot <-  factor(results_signal_all$size_tot,levels(results_signal_all$size_tot)[c(1,3,2)])
    
    ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
      xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
      scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
      scale_color_manual(values=colors)+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      facet_grid(~size_tot)
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_assymetry_A_Mantel_tests_",eco_matrix,"_", method,"_inference_size_tot_anticorrelation_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 11, height = 7)
    
    
    
    ##  Resume for clade B
    
    results_signal$Inference <- "Not significant signal"
    results_signal$Inference[which(results_signal$pvalue_high_B<0.05)] <- "Significant signal"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.05))] <- "Significant signal (R>0.05)"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.15))] <- "Significant signal (R>0.15)"
    
    colors <- c()
    if ("Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
    if ("Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
    if ("Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
    if ("Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}
    
    library(dplyr)
    results_signal_all <- results_signal %>% 
      group_by(param,Inference,size) %>% 
      summarise(count=n()) %>% 
      mutate(perc=count/sum(count))
    
    brks <- c(0, 0.25, 0.5, 0.75, 1)
    ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
      xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
      scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
      scale_color_manual(values=colors)+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      facet_grid(~size)
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_assymetry_A_Mantel_tests_",eco_matrix,"_", method,"_inference_size_clade_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 14, height = 7)
    
    
    results_signal$size_tot <- "<80"
    results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>=75)] <- "80-140"
    results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>140)] <- ">140"
    
    results_signal_all <- results_signal %>% 
      group_by(param,Inference,size_tot) %>% 
      summarise(count=n()) %>% 
      mutate(perc=count/sum(count))
    
    results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
    results_signal_all$size_tot <-  factor(results_signal_all$size_tot,levels(results_signal_all$size_tot)[c(1,3,2)])
    
    
    ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
      xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
      scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
      scale_color_manual(values=colors)+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      facet_grid(~size_tot)
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_assymetry_A_Mantel_tests_",eco_matrix,"_", method,"_inference_size_tot_clade_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 10, height = 7)
    
    
    ### Add anticorrelation
    
    results_signal$Inference <- "d) Not significant signal"
    results_signal$Inference[which(results_signal$pvalue_high_B<0.05)] <- "e) Significant signal"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.05))] <- "f) Significant signal (R>0.05)"
    results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.15))] <- "g) Significant signal (R>0.15)"
    results_signal$Inference[which(results_signal$pvalue_low_B<0.05)] <- "c) Significant anti-signal"
    results_signal$Inference[intersect(which(results_signal$pvalue_low_B<0.05), which(results_signal$mantel_cor_B<(-0.05)))] <- "b) Significant anti-signal (R<-0.05)"
    results_signal$Inference[intersect(which(results_signal$pvalue_low_B<0.05), which(results_signal$mantel_cor_B<(-0.15)))] <- "a) Significant anti-signal (R<-0.15)"
    
    #colors=c("#943126","#e74c3c","#f1948a","#f7dc6f","#7dcea0","#27ae60","#196f3d")
    colors <- c()
    if ("a) Significant anti-signal (R<-0.15)" %in% results_signal$Inference) {colors <- c(colors, "#943126")}
    if ("b) Significant anti-signal (R<-0.05)" %in% results_signal$Inference) {colors <- c(colors, "#e74c3c")}
    if ("c) Significant anti-signal" %in% results_signal$Inference) {colors <- c(colors, "#f1948a")}
    if ("d) Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
    if ("e) Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
    if ("f) Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
    if ("g) Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}
    
    
    results_signal_all <- results_signal %>% 
      group_by(param,Inference,size_tot) %>% 
      summarise(count=n()) %>% 
      mutate(perc=count/sum(count))
    
    results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
    results_signal_all$size_tot <-  factor(results_signal_all$size_tot,levels(results_signal_all$size_tot)[c(1,3,2)])
    
    ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
      xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
      scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
      scale_color_manual(values=colors)+
      scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                  bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                  bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                  bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                  bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
      facet_grid(~size_tot)
    
    ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_assymetry_A_Mantel_tests_",eco_matrix,"_", method,"_inference_size_tot_anticorrelation_clade_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 11, height = 7)
    
    print(table(results_signal$param,results_signal$size_tot))
    
  }
}



#####  Step 4-C: Plot outputs phylogenetic signal with assymetry (10%) - simul 1 - PBLM ######

# on cluster

rm(list=ls())

library(ggplot2)

setwd("/Users/bperez/ownCloud/Recherche/These/ISYEB/signal_phylo/simul_1/")
setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1")


transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))


eco_matrix="PBLM"
eco_matrix="PBLM_binary"

for (eco_matrix in c("PBLM","PBLM_binary")){
  
  info = list.files(path = "results_PBLM_assymetry/", pattern = paste0('_',eco_matrix,'.csv'))
  info <- info[grep(x=info,pattern = "seed")]
  
  results_all <- c()
  for (i in info){
    
    if (length(grep(x=i, pattern = "_neutral_")==1)>0) {param="neutral"}
    if (length(grep(x=i, pattern = "_mutualism_1_")==1)>0) {param="mutualism_1"}
    if (length(grep(x=i, pattern = "_mutualism_2_")==1)>0) {param="mutualism_2"}
    if (length(grep(x=i, pattern = "_mutualism_3_")==1)>0) {param="mutualism_3"}
    if (length(grep(x=i, pattern = "_mutualism_4_")==1)>0) {param="mutualism_4"}
    if (length(grep(x=i, pattern = "_mutualism_5_")==1)>0) {param="mutualism_5"}
    if (length(grep(x=i, pattern = "_mutualism_6_")==1)>0) {param="mutualism_6"}
    if (length(grep(x=i, pattern = "_antagonism_1_")==1)>0) {param="antagonism_1"}
    if (length(grep(x=i, pattern = "_antagonism_2_")==1)>0) {param="antagonism_2"}
    if (length(grep(x=i, pattern = "_antagonism_3_")==1)>0) {param="antagonism_3"}
    if (length(grep(x=i, pattern = "_antagonism_4_")==1)>0) {param="antagonism_4"}
    if (length(grep(x=i, pattern = "_antagonism_5_")==1)>0) {param="antagonism_5"}
    if (length(grep(x=i, pattern = "_antagonism_6_")==1)>0) {param="antagonism_6"}
    if (length(grep(x=i, pattern = "_antagonism_7_")==1)>0) {param="antagonism_7"}
    if (length(grep(x=i, pattern = "_antagonism_8_")==1)>0) {param="antagonism_8"}
    if (length(grep(x=i, pattern = "_antagonism_9_")==1)>0) {param="antagonism_9"}
    
    if (length(grep(x=i, pattern = "_A_")==1)>0) {size="3000"}
    if (length(grep(x=i, pattern = "_B_")==1)>0) {size="4000"}
    if (length(grep(x=i, pattern = "_C_")==1)>0) {size="5000"}
    if (length(grep(x=i, pattern = "_D_")==1)>0) {size="2000"}
    if (length(grep(x=i, pattern = "_E_")==1)>0) {size="1000"}
    if (length(grep(x=i, pattern = "_F_")==1)>0) {size="500"}
    
    results_seed <- read.table(paste0("results_PBLM_assymetry/",i), header=T, sep=";")
    
    # pb with D and E
    if (ncol(results_seed)==1) {results_seed = t(results_seed)}
    
    #if (results_seed[4]!="inference failed"){
    results_all <- rbind(results_all, c(param, size, results_seed))
    #}
  }
  
  results_all <- data.frame(results_all)
  colnames(results_all) <- c("param","size", "nb_A", "nb_B", "dA", "dB", "MSETotal", "MSEFull", "MSEStar", "MSEBase")
  
  results_all$param <- as.character(results_all$param)
  results_all$size <- as.character(results_all$size)
  results_all$nb_A <- as.numeric(as.character(results_all$nb_A))
  results_all$nb_B <- as.numeric(as.character(results_all$nb_B))
  results_all$dA <- as.numeric(as.character(results_all$dA))
  results_all$dB <- as.numeric(as.character(results_all$dB))
  results_all$MSETotal <- as.numeric(as.character(results_all$MSETotal))
  results_all$MSEFull <- as.numeric(as.character(results_all$MSEFull))
  results_all$MSEStar <- as.numeric(as.character(results_all$MSEStar))
  results_all$MSEBase <- as.numeric(as.character(results_all$MSEBase))
  
  results_all$size <- as.factor(results_all$size)
  results_all$size <-  factor(results_all$size,levels(results_all$size)[c(5,1,2,3,4,6)])
  
  
  results_all$Diff_MSE <- results_all$MSEStar - results_all$MSEFull
  
  results_all$Ratio_MSE <- (results_all$MSEStar - results_all$MSEFull)/results_all$MSEStar
  
  write.table(results_all,paste0("results_phylogenetic_signal_simul_1_assymetry_A_",eco_matrix,".csv"), quote=F, sep=";",row.names=F)
  
  
  results_all_failed <- results_all
  results_all <- results_all[which(results_all$dA!="inference failed"),]
  
  # Plot 
  
  results_all_failed$Inference <- "Inference failed"
  results_all_failed$Inference[which(results_all_failed$Ratio_MSE<=0)] <- "Not significant signal"
  results_all_failed$Inference[which(results_all_failed$Ratio_MSE>0)] <- "Significant signal"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dA>0.05))] <- "Significant signal (d>0.05)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dB>0.05))] <- "Significant signal (d>0.05)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dA>0.15))] <- "Significant signal (d>0.15)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dB>0.15))] <- "Significant signal (d>0.15)"
  
  colors <- c()
  if ("Not significant signal" %in% results_all_failed$Inference) {colors <- c(colors, "#f7dc6f")}
  if ("Significant signal" %in% results_all_failed$Inference) {colors <- c(colors, "#7dcea0")}
  if ("Significant signal (d>0.05)" %in% results_all_failed$Inference) {colors <- c(colors, "#27ae60")}
  if ("Significant signal (d>0.15)" %in% results_all_failed$Inference) {colors <- c(colors, "#196f3d")}
  
  library(dplyr)
  results_all_data <- results_all_failed %>% 
    group_by(param,Inference,size) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  brks <- c(0, 0.25, 0.5, 0.75, 1)
  ggplot(results_all_data,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=colors)+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    facet_grid(~size)
  
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_assymetry_A_",eco_matrix,"_inference_size.pdf"), device=cairo_pdf, plot = last_plot(), width = 14, height = 7)
  
  
  ### More stringent
  
  results_all_failed$Inference <- "Inference failed"
  results_all_failed$Inference[which(results_all_failed$Ratio_MSE<=0.01)] <- "Not significant signal"
  results_all_failed$Inference[which(results_all_failed$Ratio_MSE>0.01)] <- "Significant signal"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0.01), which(results_all_failed$dA>0.05))] <- "Significant signal (d>0.05)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0.01), which(results_all_failed$dB>0.05))] <- "Significant signal (d>0.05)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0.01), which(results_all_failed$dA>0.15))] <- "Significant signal (d>0.15)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0.01), which(results_all_failed$dB>0.15))] <- "Significant signal (d>0.15)"
  
  colors <- c()
  if ("Not significant signal" %in% results_all_failed$Inference) {colors <- c(colors, "#f7dc6f")}
  if ("Significant signal" %in% results_all_failed$Inference) {colors <- c(colors, "#7dcea0")}
  if ("Significant signal (d>0.05)" %in% results_all_failed$Inference) {colors <- c(colors, "#27ae60")}
  if ("Significant signal (d>0.15)" %in% results_all_failed$Inference) {colors <- c(colors, "#196f3d")}
  
  
  library(dplyr)
  results_all_data <- results_all_failed %>% 
    group_by(param,Inference,size) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  brks <- c(0, 0.25, 0.5, 0.75, 1)
  ggplot(results_all_data,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=colors)+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    facet_grid(~size)
  
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_assymetry_A_",eco_matrix,"_inference_size_stringent.pdf"), device=cairo_pdf, plot = last_plot(), width = 14, height = 7)
  
  results_all_failed$Inference <- "Inference failed"
  results_all_failed$Inference[which(results_all_failed$Ratio_MSE<=0)] <- "Not significant signal"
  results_all_failed$Inference[which(results_all_failed$Ratio_MSE>0)] <- "Significant signal"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dA>0.05))] <- "Significant signal (d>0.05)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dB>0.05))] <- "Significant signal (d>0.05)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dA>0.15))] <- "Significant signal (d>0.15)"
  results_all_failed$Inference[intersect(which(results_all_failed$Ratio_MSE>0), which(results_all_failed$dB>0.15))] <- "Significant signal (d>0.15)"
  
  results_all_failed$size_tot <- "<80"
  results_all_failed$size_tot[which(results_all_failed$nb_A+results_all_failed$nb_B>=80)] <- "80-140"
  results_all_failed$size_tot[which(results_all_failed$nb_A+results_all_failed$nb_B>140)] <- ">140"
  
  colors <- c()
  if ("Not significant signal" %in% results_all_failed$Inference) {colors <- c(colors, "#f7dc6f")}
  if ("Significant signal" %in% results_all_failed$Inference) {colors <- c(colors, "#7dcea0")}
  if ("Significant signal (d>0.05)" %in% results_all_failed$Inference) {colors <- c(colors, "#27ae60")}
  if ("Significant signal (d>0.15)" %in% results_all_failed$Inference) {colors <- c(colors, "#196f3d")}
  
  results_signal_all <- results_all_failed %>% 
    group_by(param,Inference,size_tot) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
  results_signal_all$size_tot <-  factor(results_signal_all$size_tot,c("<80"  ,  "80-140", ">140"  ))
  
  
  ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=colors)+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    facet_grid(~size_tot)
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_assymetry_A_",eco_matrix,"_inference_size_tot_stringent.pdf"), device=cairo_pdf, plot = last_plot(), width = 10, height = 7)
  
  
}

# scp bperez@jord.biologie.ens.fr:/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1/*assymetry*.pdf  /Users/bperez/ownCloud/Recherche/These/ISYEB/signal_phylo/simul_1/




#### Step 5-A: Add Noise in UniFrac with Normal noise #####

# on cluster   test_noise.R

rm(list=ls())

library(ggplot2)
library(vegan)
library(phytools)
library(ggplot2)

library(GUniFrac,lib.loc="/users/biodiv/bperez/packages/")
library(parallel)
require(R.utils)

setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1")

dyn.load("../script/permute.so")
source("../script/function_phylo_signal_network.R")

transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))

info = list.files(pattern = paste0('.csv'))
info <- info[grep(x=info,pattern = "network")]
info <- info[grep(x=info,pattern = "seed")]


results_all <- c()

name="simul_1"
method="GUniFrac"
correlation="Pearson"
nperm=10000

#sigma in 0.25, 0.5, 0.75, 1
sigma=0.25

for (sigma in c(0.25, 0.5, 0.75, 1)){
  for (i in info){
    
    if (length(grep(x=i, pattern = "_neutral_")==1)>0) {param="neutral"}
    if (length(grep(x=i, pattern = "_mutualism_1_")==1)>0) {param="mutualism_1"}
    if (length(grep(x=i, pattern = "_mutualism_2_")==1)>0) {param="mutualism_2"}
    if (length(grep(x=i, pattern = "_mutualism_3_")==1)>0) {param="mutualism_3"}
    if (length(grep(x=i, pattern = "_mutualism_4_")==1)>0) {param="mutualism_4"}
    if (length(grep(x=i, pattern = "_mutualism_5_")==1)>0) {param="mutualism_5"}
    if (length(grep(x=i, pattern = "_mutualism_6_")==1)>0) {param="mutualism_6"}
    if (length(grep(x=i, pattern = "_antagonism_1_")==1)>0) {param="antagonism_1"}
    if (length(grep(x=i, pattern = "_antagonism_2_")==1)>0) {param="antagonism_2"}
    if (length(grep(x=i, pattern = "_antagonism_3_")==1)>0) {param="antagonism_3"}
    if (length(grep(x=i, pattern = "_antagonism_4_")==1)>0) {param="antagonism_4"}
    if (length(grep(x=i, pattern = "_antagonism_5_")==1)>0) {param="antagonism_5"}
    if (length(grep(x=i, pattern = "_antagonism_6_")==1)>0) {param="antagonism_6"}
    if (length(grep(x=i, pattern = "_antagonism_7_")==1)>0) {param="antagonism_7"}
    if (length(grep(x=i, pattern = "_antagonism_8_")==1)>0) {param="antagonism_8"}
    if (length(grep(x=i, pattern = "_antagonism_9_")==1)>0) {param="antagonism_9"}
    
    if (length(grep(x=i, pattern = "_A_")==1)>0) {size_tot="3000"}
    if (length(grep(x=i, pattern = "_B_")==1)>0) {size_tot="4000"}
    if (length(grep(x=i, pattern = "_C_")==1)>0) {size_tot="5000"}
    if (length(grep(x=i, pattern = "_D_")==1)>0) {size_tot="2000"}
    if (length(grep(x=i, pattern = "_E_")==1)>0) {size_tot="1000"}
    if (length(grep(x=i, pattern = "_F_")==1)>0) {size_tot="500"}
    
    if (length(grep(x=i, pattern = "_A_")==1)>0) {size="A"}
    if (length(grep(x=i, pattern = "_B_")==1)>0) {size="B"}
    if (length(grep(x=i, pattern = "_C_")==1)>0) {size="C"}
    if (length(grep(x=i, pattern = "_D_")==1)>0) {size="D"}
    if (length(grep(x=i, pattern = "_E_")==1)>0) {size="E"}
    if (length(grep(x=i, pattern = "_F_")==1)>0) {size="F"}
    
    network <- read.table(i, header=T,sep=";")
    
    colnames(network) <- seq(1:ncol(network))
    rownames(network) <- seq(1:nrow(network))
    
    #"network_simul_1_E_neutral_seed_9.csv"
    
    
    tree_A <- read.tree(paste0("network_tree_guild_A_",gsub(x=gsub(x=i,pattern="network_",replacement = ""),pattern=".csv",replacement = ""),".tre"))
    tree_B <- read.tree(paste0("network_tree_guild_B_",gsub(x=gsub(x=i,pattern="network_",replacement = ""),pattern=".csv",replacement = ""),".tre"))
    
    # Change B - Signal on A
    set.seed(1)
    cophe_B <- cophenetic(tree_B)
    cophe_B_noise <- as.matrix(as.dist(cophe_B + matrix(rnorm(n=nrow(cophe_B)*ncol(cophe_B), mean = 0, sd= sigma*sd(cophe_B) ),nrow = nrow(cophe_B))))
    tree_B_noise <- nj(cophe_B_noise)
    output_A <- phylosignal_network(network, tree_A, tree_B_noise, method = method, nperm = nperm, correlation = correlation)
    
    set.seed(1)
    cophe_A <- cophenetic(tree_A)
    cophe_A_noise <- as.matrix(as.dist(cophe_A + matrix(rnorm(n=nrow(cophe_A)*ncol(cophe_A), mean = 0, sd= sigma*sd(cophe_A) ),nrow = nrow(cophe_A))))
    tree_A_noise <- nj(cophe_A_noise)
    output_B <- phylosignal_network(network, tree_A_noise, tree_B, method = method, nperm = nperm, correlation = correlation)
    
    results_all <- rbind(results_all, c(param, size, size_tot, sigma, output_A[3:5], output_B[6:8] ))
  }
  
  
  results_all <- data.frame(results_all)
  colnames(results_all) <- c("param","size_name","size",  "sigma", "mantel_cor_A", "pvalue_high_A", "pvalue_low_A", "mantel_cor_B", "pvalue_high_B", "pvalue_low_B")
  
  results_all$param <- as.character(results_all$param)
  results_all$size <- as.character(results_all$size)
  results_all$mantel_cor_A <- as.numeric(as.character(results_all$mantel_cor_A))
  results_all$pvalue_high_A <- as.numeric(as.character(results_all$pvalue_high_A))
  results_all$mantel_cor_B <- as.numeric(as.character(results_all$mantel_cor_B))
  results_all$pvalue_high_B <- as.numeric(as.character(results_all$pvalue_high_B))
  
  # Plot 
  results_signal <- results_all
  
  # Size plot
  ggplot(results_signal,aes(x=param,y=mantel_cor_A,fill=size))+xlab("Parameters")+ylab("Mantel correlation (A)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
    geom_hline(yintercept=0, linetype="dashed", color="#273746")+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    geom_boxplot(alpha=0.9,size=0.75, color="#273746")
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Noise_Mantel_tests_",method,"_", correlation,"_sigma_", sigma ,"_size_correlation_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)
  
  ggplot(results_signal,aes(x=param,y=mantel_cor_B,fill=size))+xlab("Parameters")+ylab("Mantel correlation (B)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
    geom_hline(yintercept=0, linetype="dashed", color="#273746")+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    geom_boxplot(alpha=0.9,size=0.75, color="#273746")
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Noise_Mantel_tests_",method,"_", correlation,"_sigma_", sigma ,"_size_correlation_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)
  
  
  ggplot(results_signal,aes(x=param,y=pvalue_high_A,fill=size))+xlab("Parameters")+ylab("pvalue - Mantel correlation (A)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
    geom_hline(yintercept=0, linetype="dashed", color="#273746")+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    geom_boxplot(alpha=0.9,size=0.75, color="#273746")
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Noise_Mantel_tests_",method,"_", correlation,"_sigma_", sigma ,"_size_pvalue_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)
  
  ggplot(results_signal,aes(x=param,y=pvalue_high_B,fill=size))+xlab("Parameters")+ylab("pvalue - Mantel correlation (B)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
    geom_hline(yintercept=0, linetype="dashed", color="#273746")+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    geom_boxplot(alpha=0.9,size=0.75, color="#273746")
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Noise_Mantel_tests_",method,"_", correlation,"_sigma_", sigma ,"_size_pvalue_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)
  
  
  
  ##  Resume for clade A
  
  results_signal$Inference <- "Not significant signal"
  results_signal$Inference[which(results_signal$pvalue_high_A<0.05)] <- "Significant signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.05))] <- "Significant signal (R>0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.15))] <- "Significant signal (R>0.15)"
  
  library(dplyr)
  results_signal_all <- results_signal %>% 
    group_by(param,Inference,size) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  brks <- c(0, 0.25, 0.5, 0.75, 1)
  ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=c("#f7dc6f","#7dcea0","#27ae60","#196f3d"))+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=c("#f7dc6f","#7dcea0","#27ae60","#196f3d"))+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    facet_grid(~size)
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Noise_Mantel_tests_",method,"_", correlation,"_sigma_", sigma ,"_inference_size_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 14, height = 7)
  
  
  ##  Resume for clade B
  
  results_signal$Inference <- "Not significant signal"
  results_signal$Inference[which(results_signal$pvalue_high_B<0.05)] <- "Significant signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.05))] <- "Significant signal (R>0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.15))] <- "Significant signal (R>0.15)"
  
  library(dplyr)
  results_signal_all <- results_signal %>% 
    group_by(param,Inference,size) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  brks <- c(0, 0.25, 0.5, 0.75, 1)
  ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=c("#f7dc6f","#7dcea0","#27ae60","#196f3d"))+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=c("#f7dc6f","#7dcea0","#27ae60","#196f3d"))+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    facet_grid(~size)
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Noise_Mantel_tests_",method,"_", correlation,"_sigma_", sigma ,"_inference_size_clade_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 14, height = 7)
  
}


#### Step 5-B: Add Noise in UniFrac with simulations of sequences #####

# on cluster

rm(list=ls())

library(ggplot2)
library(vegan)
library(phytools)
library(ggplot2)
library(GUniFrac,lib.loc="/users/biodiv/bperez/packages/")
library(parallel)
require(R.utils)

setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1")

dyn.load("../script/permute.so")
source("../script/function_phylo_signal_network.R")

dir.create(file.path("./simulations/"), showWarnings = FALSE)


simulate_alignment <- function(name, host_tree, mu=1.5, N=300,proportion_variant=0.3){
  
  host_tree <- ladderize(host_tree)
  host_tree <- multi2di(host_tree)
  host_tree$edge.length <- host_tree$edge.length/sum(host_tree$edge.length) #host tree scaled with total branch length=1
  host_tree$edge.length[host_tree$edge.length==0] <- 0.00001
  host_tree <- force.ultrametric(host_tree,method = "extend")
  
  host_tree$edge.length <- host_tree$edge.length/sum(host_tree$edge.length) #host tree scaled with total branch length=1
  n <- Ntip(host_tree)
  
  set.seed(3)
  
  simulated_mu <- mu
  maxlen <- max(node.depth.edgelength(host_tree))
  tree <- host_tree
  
  maxlen <- max(node.depth.edgelength(tree))
  a <- 1
  b <- 4
  c <- 1
  d <- 1
  e <- 4
  f <- 1
  propinv <- c(0.25,0.25,0.25,0.25)
  Q <-  (as.matrix(rbind(c(0,a,b,c)*t(propinv),c(a,0,d,e)*t(propinv),c(b,d,0,f)*t(propinv),c(c,e,f,0)*t(propinv))))
  diag(Q) <- - apply(Q,1,sum)
  Q <- -Q/Q[4,4]
  eigQ <- eigen(Q)
  eig_val <- eigQ$values
  eig_vect <- eigQ$vectors
  ivp <- solve(eig_vect)
  nodes<-order(node.depth.edgelength(tree),decreasing=F)[-1]
  original_sequences <-  matrix(0,nrow=n,ncol=N)
  for (nu in 1:N){
    if (runif(1,0,1)<proportion_variant){
      L <-  matrix(0,nrow=(2*n-1),ncol=4)
      L[n+1,sample(1:4,size=1)] <- 1
      for(i in nodes) {
        v <- tree$edge[which(tree$edge[,2]==i),1]
        t <- as.numeric(tree$edge.length[which(tree$edge[,2]==i)])
        proba_nu <-  L[v,]%*%eigQ$vectors%*%diag(exp(t*mu*eigQ$values))%*%ivp
        L[i,sample(1:4,size=1,prob=proba_nu)] <- 1}
      n_seq <- rep(c("a","c","g","t"),n)[as.logical(t(L[1:n,]))]
      original_sequences[,nu]<- t(t(n_seq))
    } else {original_sequences[1:n,nu] <- sample(c("a","c","g","t"),size=1)}}
  row.names(original_sequences) <- tree$tip.label
  
  return(original_sequences)
}


transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))

info = list.files(pattern = paste0('.csv'))
info <- info[grep(x=info,pattern = "network")]
info <- info[grep(x=info,pattern = "seed")]


name="simul_1"
method="GUniFrac"
correlation="Pearson"
nperm=10000

N=150

results_noise <- function(N){
  
  print(N)
  
  result_noise_i <- function(i,N){
    
    if (length(grep(x=i, pattern = "_neutral_")==1)>0) {param="neutral"}
    if (length(grep(x=i, pattern = "_mutualism_1_")==1)>0) {param="mutualism_1"}
    if (length(grep(x=i, pattern = "_mutualism_2_")==1)>0) {param="mutualism_2"}
    if (length(grep(x=i, pattern = "_mutualism_3_")==1)>0) {param="mutualism_3"}
    if (length(grep(x=i, pattern = "_mutualism_4_")==1)>0) {param="mutualism_4"}
    if (length(grep(x=i, pattern = "_mutualism_5_")==1)>0) {param="mutualism_5"}
    if (length(grep(x=i, pattern = "_mutualism_6_")==1)>0) {param="mutualism_6"}
    if (length(grep(x=i, pattern = "_antagonism_1_")==1)>0) {param="antagonism_1"}
    if (length(grep(x=i, pattern = "_antagonism_2_")==1)>0) {param="antagonism_2"}
    if (length(grep(x=i, pattern = "_antagonism_3_")==1)>0) {param="antagonism_3"}
    if (length(grep(x=i, pattern = "_antagonism_4_")==1)>0) {param="antagonism_4"}
    if (length(grep(x=i, pattern = "_antagonism_5_")==1)>0) {param="antagonism_5"}
    if (length(grep(x=i, pattern = "_antagonism_6_")==1)>0) {param="antagonism_6"}
    if (length(grep(x=i, pattern = "_antagonism_7_")==1)>0) {param="antagonism_7"}
    if (length(grep(x=i, pattern = "_antagonism_8_")==1)>0) {param="antagonism_8"}
    if (length(grep(x=i, pattern = "_antagonism_9_")==1)>0) {param="antagonism_9"}
    
    if (length(grep(x=i, pattern = "_A_")==1)>0) {size_tot="3000"}
    if (length(grep(x=i, pattern = "_B_")==1)>0) {size_tot="4000"}
    if (length(grep(x=i, pattern = "_C_")==1)>0) {size_tot="5000"}
    if (length(grep(x=i, pattern = "_D_")==1)>0) {size_tot="2000"}
    if (length(grep(x=i, pattern = "_E_")==1)>0) {size_tot="1000"}
    if (length(grep(x=i, pattern = "_F_")==1)>0) {size_tot="500"}
    
    if (length(grep(x=i, pattern = "_A_")==1)>0) {size="A"}
    if (length(grep(x=i, pattern = "_B_")==1)>0) {size="B"}
    if (length(grep(x=i, pattern = "_C_")==1)>0) {size="C"}
    if (length(grep(x=i, pattern = "_D_")==1)>0) {size="D"}
    if (length(grep(x=i, pattern = "_E_")==1)>0) {size="E"}
    if (length(grep(x=i, pattern = "_F_")==1)>0) {size="F"}
    
    network <- read.table(i, header=T,sep=";")
    
    colnames(network) <- seq(1:ncol(network))
    rownames(network) <- seq(1:nrow(network))
    
    tree_A <- read.tree(paste0("network_tree_guild_A_",gsub(x=gsub(x=i,pattern="network_",replacement = ""),pattern=".csv",replacement = ""),".tre"))
    tree_B <- read.tree(paste0("network_tree_guild_B_",gsub(x=gsub(x=i,pattern="network_",replacement = ""),pattern=".csv",replacement = ""),".tre"))
    
    nb_A <- Ntip(tree_A)
    nb_B <- Ntip(tree_B)
    
    name=gsub(x=gsub(x=i,pattern="network_",replacement = ""),pattern=".csv",replacement = "")
    
    # Change B - Signal on A
    set.seed(3)
    alignment <- simulate_alignment(name, host_tree=tree_B, N=N)
    tree_B_noise <- midpoint.root(nj(dist.dna(as.DNAbin(alignment), model = "K80")))
    write.dna(alignment,paste("simulations/alignment_B_",name,"_",N,"bp.fas",sep=""), format = "fasta",nbcol = -1,colsep="",colw=N)
    write.tree(tree_B_noise, paste("simulations/tree_noise_B_",name,"_",N,"bp.tre",sep=""))
    output_A <- phylosignal_network(network, tree_A, tree_B_noise, method = method, nperm = nperm, correlation = correlation)
    
    set.seed(3)
    alignment <- simulate_alignment(name, host_tree=tree_A, N=N)
    tree_A_noise <- midpoint.root(nj(dist.dna(as.DNAbin(alignment), model = "K80")))
    write.dna(alignment,paste("simulations/alignment_A_",name,"_",N,"bp.fas",sep=""), format = "fasta",nbcol = -1,colsep="",colw=N)
    write.tree(tree_A_noise, paste("simulations/tree_noise_A_",name,"_",N,"bp.tre",sep=""))
    output_B <- phylosignal_network(network, tree_A_noise, tree_B, method = method, nperm = nperm, correlation = correlation)
    
    #results_all <- rbind(results_all, c(param, size, size_tot, nb_A,nb_B, N, output_A[3:5], output_B[6:8] ))
    return(c(param, size, size_tot, nb_A,nb_B, N, output_A[3:5], output_B[6:8] ))
  }
  
  results_all <- matrix(unlist(mclapply(info,result_noise_i,mc.cores=5,mc.preschedule=F, N=N)),nrow=length(info),byrow=T)
  results_all <- data.frame(results_all)
  colnames(results_all) <- c("param","size_name","size", "nb_A", "nb_B",  "N", "mantel_cor_A", "pvalue_high_A", "pvalue_low_A", "mantel_cor_B", "pvalue_high_B", "pvalue_low_B")
  print("worked")
  
  results_all$param <- as.character(results_all$param)
  results_all$size <- as.character(results_all$size)
  results_all$mantel_cor_A <- as.numeric(as.character(results_all$mantel_cor_A))
  results_all$pvalue_high_A <- as.numeric(as.character(results_all$pvalue_high_A))
  results_all$mantel_cor_B <- as.numeric(as.character(results_all$mantel_cor_B))
  results_all$pvalue_high_B <- as.numeric(as.character(results_all$pvalue_high_B))
  
  write.table(results_all, paste0("results_phylogenetic_signal_simul_1_Noise_Mantel_tests_",method,"_", correlation,"_N_", N ,".csv"), sep=";", quote=F,row.names = F)
  print(N)
}

mclapply(c(75, 150, 300, 600, 1200),results_noise,mc.cores = 5)


# Plot results

for (N in c(75, 150, 300, 600, 1200)){
  
  results_all <- read.table(paste0("results_phylogenetic_signal_simul_1_Noise_Mantel_tests_",method,"_", correlation,"_N_", N ,".csv"), sep=";", header=T)
  
  results_all$param <- as.character(results_all$param)
  results_all$size <- as.character(results_all$size)
  results_all$mantel_cor_A <- as.numeric(as.character(results_all$mantel_cor_A))
  results_all$pvalue_high_A <- as.numeric(as.character(results_all$pvalue_high_A))
  results_all$mantel_cor_B <- as.numeric(as.character(results_all$mantel_cor_B))
  results_all$pvalue_high_B <- as.numeric(as.character(results_all$pvalue_high_B))
  results_signal <- results_all
  
  ##  Resume for clade A
  
  results_signal$Inference <- "Not significant signal"
  results_signal$Inference[which(results_signal$pvalue_high_A<0.05)] <- "Significant signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.05))] <- "Significant signal (R>0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.15))] <- "Significant signal (R>0.15)"
  
  library(dplyr)
  results_signal_all <- results_signal %>% 
    group_by(param,Inference,size) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  brks <- c(0, 0.25, 0.5, 0.75, 1)
  ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=c("#f7dc6f","#7dcea0","#27ae60","#196f3d"))+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=c("#f7dc6f","#7dcea0","#27ae60","#196f3d"))+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    facet_grid(~size)
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Noise_Mantel_tests_",method,"_", correlation,"_N_", N ,"_inference_size_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 14, height = 7)
  
  
  results_signal$size_tot <- "<150"
  results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>=150)] <- "150-250"
  results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>250)] <- ">250"
  
  results_signal_all <- results_signal %>% 
    group_by(param,Inference,size_tot) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
  results_signal_all$size_tot <-  factor(results_signal_all$size_tot,levels(results_signal_all$size_tot)[c(1,3,2)])
  
  
  ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=c("#f7dc6f","#7dcea0","#27ae60","#196f3d"))+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=c("#f7dc6f","#7dcea0","#27ae60","#196f3d"))+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    facet_grid(~size_tot)
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Noise_Mantel_tests_",method,"_", correlation,"_N_", N ,"_inference_size_tot_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 10, height = 7)
  
  
  ### Add anticorrelation
  
  results_signal$Inference <- "d) Not significant signal"
  results_signal$Inference[which(results_signal$pvalue_high_A<0.05)] <- "e) Significant signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.05))] <- "f) Significant signal (R>0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.15))] <- "g) Significant signal (R>0.15)"
  results_signal$Inference[which(results_signal$pvalue_low_A<0.05)] <- "c) Significant anti-signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.05)))] <- "b) Significant anti-signal (R<-0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.15)))] <- "a) Significant anti-signal (R<-0.15)"
  
  colors=c("#943126","#e74c3c","#f1948a","#f7dc6f","#7dcea0","#27ae60","#196f3d")
  
  results_signal_all <- results_signal %>% 
    group_by(param,Inference,size_tot) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
  results_signal_all$size_tot <-  factor(results_signal_all$size_tot,levels(results_signal_all$size_tot)[c(1,3,2)])
  
  ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=colors)+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    facet_grid(~size_tot)
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Noise_Mantel_tests_",method,"_", correlation,"_N_", N ,"_inference_size_tot_anticorrelation_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 11, height = 7)
  
  
  
  ##  Resume for clade B
  
  results_signal$Inference <- "Not significant signal"
  results_signal$Inference[which(results_signal$pvalue_high_B<0.05)] <- "Significant signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.05))] <- "Significant signal (R>0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.15))] <- "Significant signal (R>0.15)"
  
  library(dplyr)
  results_signal_all <- results_signal %>% 
    group_by(param,Inference,size) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  brks <- c(0, 0.25, 0.5, 0.75, 1)
  ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=c("#f7dc6f","#7dcea0","#27ae60","#196f3d"))+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=c("#f7dc6f","#7dcea0","#27ae60","#196f3d"))+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    facet_grid(~size)
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Noise_Mantel_tests_",method,"_", correlation,"_N_", N ,"_inference_size_clade_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 14, height = 7)
  
  
  results_signal$size_tot <- "<150"
  results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>=150)] <- "150-250"
  results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>250)] <- ">250"
  
  results_signal_all <- results_signal %>% 
    group_by(param,Inference,size_tot) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
  results_signal_all$size_tot <-  factor(results_signal_all$size_tot,levels(results_signal_all$size_tot)[c(1,3,2)])
  
  
  ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=c("#f7dc6f","#7dcea0","#27ae60","#196f3d"))+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=c("#f7dc6f","#7dcea0","#27ae60","#196f3d"))+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    facet_grid(~size_tot)
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Noise_Mantel_tests_",method,"_", correlation,"_N_", N ,"_inference_size_tot_clade_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 10, height = 7)
  
  
  ### Add anticorrelation
  
  results_signal$Inference <- "d) Not significant signal"
  results_signal$Inference[which(results_signal$pvalue_high_B<0.05)] <- "e) Significant signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.05))] <- "f) Significant signal (R>0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.15))] <- "g) Significant signal (R>0.15)"
  results_signal$Inference[which(results_signal$pvalue_low_B<0.05)] <- "c) Significant anti-signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_low_B<0.05), which(results_signal$mantel_cor_B<(-0.05)))] <- "b) Significant anti-signal (R<-0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_low_B<0.05), which(results_signal$mantel_cor_B<(-0.15)))] <- "a) Significant anti-signal (R<-0.15)"
  
  colors=c("#943126","#e74c3c","#f1948a","#f7dc6f","#7dcea0","#27ae60","#196f3d")
  
  results_signal_all <- results_signal %>% 
    group_by(param,Inference,size_tot) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
  results_signal_all$size_tot <-  factor(results_signal_all$size_tot,levels(results_signal_all$size_tot)[c(1,3,2)])
  
  ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=colors)+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    facet_grid(~size_tot)
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_simul_1_Noise_Mantel_tests_",method,"_", correlation,"_N_", N ,"_inference_size_tot_anticorrelation_clade_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 11, height = 7)
  
}


tree_B <- ladderize(read.tree("network_tree_guild_B_simul_1_A_mutualism_5_seed_7.tre"))
tree_B_noise_75 <- ladderize(read.tree(paste("simulations/tree_noise_B_simul_1_A_mutualism_5_seed_7_75bp.tre",sep="")))
tree_B_noise_150 <- ladderize(read.tree(paste("simulations/tree_noise_B_simul_1_A_mutualism_5_seed_7_150bp.tre",sep="")))
tree_B_noise_300 <- ladderize(read.tree(paste("simulations/tree_noise_B_simul_1_A_mutualism_5_seed_7_300bp.tre",sep="")))
tree_B_noise_600 <- ladderize(read.tree(paste("simulations/tree_noise_B_simul_1_A_mutualism_5_seed_7_600bp.tre",sep="")))
tree_B_noise_1200 <- ladderize(read.tree(paste("simulations/tree_noise_B_simul_1_A_mutualism_5_seed_7_1200bp.tre",sep="")))


pdf("phylo_matrix_test_network_simul_1_A_mutualism_5_seed_7_guild_B.pdf")
image(cophenetic(tree_B))
image(cophenetic(tree_B_noise_1200)[rownames(cophenetic(tree_B)),colnames(cophenetic(tree_B))])
image(cophenetic(tree_B_noise_600)[rownames(cophenetic(tree_B)),colnames(cophenetic(tree_B))])
image(cophenetic(tree_B_noise_300)[rownames(cophenetic(tree_B)),colnames(cophenetic(tree_B))])
image(cophenetic(tree_B_noise_150)[rownames(cophenetic(tree_B)),colnames(cophenetic(tree_B))])
image(cophenetic(tree_B_noise_75)[rownames(cophenetic(tree_B)),colnames(cophenetic(tree_B))])
dev.off()


###### Setp 5-C: Check correlation between original tree and simulated trees #####

rm(list=ls()) # on cluster

library(ggplot2)
library(vegan)
library(phytools)
library(ggplot2)
library(GUniFrac,lib.loc="/users/biodiv/bperez/packages/")
library(parallel)
require(R.utils)

setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1")

dyn.load("../script/permute.so")
source("../script/function_phylo_signal_network.R")

dir.create(file.path("./simulations/"), showWarnings = FALSE)



transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))

info = list.files(pattern = paste0('.csv'))
info <- info[grep(x=info,pattern = "network")]
info <- info[grep(x=info,pattern = "seed")]


name="simul_1"
method="GUniFrac"
correlation="Pearson"
nperm=10000

N=150

for (N in c(75, 150, 300, 600, 1200)){
  #results_noise <- function(N){
  
  results_all <- c()
  
  for (i in info){
    
    if (length(grep(x=i, pattern = "_neutral_")==1)>0) {param="neutral"}
    if (length(grep(x=i, pattern = "_mutualism_1_")==1)>0) {param="mutualism_1"}
    if (length(grep(x=i, pattern = "_mutualism_2_")==1)>0) {param="mutualism_2"}
    if (length(grep(x=i, pattern = "_mutualism_3_")==1)>0) {param="mutualism_3"}
    if (length(grep(x=i, pattern = "_mutualism_4_")==1)>0) {param="mutualism_4"}
    if (length(grep(x=i, pattern = "_mutualism_5_")==1)>0) {param="mutualism_5"}
    if (length(grep(x=i, pattern = "_mutualism_6_")==1)>0) {param="mutualism_6"}
    if (length(grep(x=i, pattern = "_antagonism_1_")==1)>0) {param="antagonism_1"}
    if (length(grep(x=i, pattern = "_antagonism_2_")==1)>0) {param="antagonism_2"}
    if (length(grep(x=i, pattern = "_antagonism_3_")==1)>0) {param="antagonism_3"}
    if (length(grep(x=i, pattern = "_antagonism_4_")==1)>0) {param="antagonism_4"}
    if (length(grep(x=i, pattern = "_antagonism_5_")==1)>0) {param="antagonism_5"}
    if (length(grep(x=i, pattern = "_antagonism_6_")==1)>0) {param="antagonism_6"}
    if (length(grep(x=i, pattern = "_antagonism_7_")==1)>0) {param="antagonism_7"}
    if (length(grep(x=i, pattern = "_antagonism_8_")==1)>0) {param="antagonism_8"}
    if (length(grep(x=i, pattern = "_antagonism_9_")==1)>0) {param="antagonism_9"}
    
    if (length(grep(x=i, pattern = "_A_")==1)>0) {size_tot="3000"}
    if (length(grep(x=i, pattern = "_B_")==1)>0) {size_tot="4000"}
    if (length(grep(x=i, pattern = "_C_")==1)>0) {size_tot="5000"}
    if (length(grep(x=i, pattern = "_D_")==1)>0) {size_tot="2000"}
    if (length(grep(x=i, pattern = "_E_")==1)>0) {size_tot="1000"}
    if (length(grep(x=i, pattern = "_F_")==1)>0) {size_tot="500"}
    
    if (length(grep(x=i, pattern = "_A_")==1)>0) {size="A"}
    if (length(grep(x=i, pattern = "_B_")==1)>0) {size="B"}
    if (length(grep(x=i, pattern = "_C_")==1)>0) {size="C"}
    if (length(grep(x=i, pattern = "_D_")==1)>0) {size="D"}
    if (length(grep(x=i, pattern = "_E_")==1)>0) {size="E"}
    if (length(grep(x=i, pattern = "_F_")==1)>0) {size="F"}
    
    network <- read.table(i, header=T,sep=";")
    
    colnames(network) <- seq(1:ncol(network))
    rownames(network) <- seq(1:nrow(network))
    
    #"network_simul_1_E_neutral_seed_9.csv"
    
    
    tree_A <- read.tree(paste0("network_tree_guild_A_",gsub(x=gsub(x=i,pattern="network_",replacement = ""),pattern=".csv",replacement = ""),".tre"))
    tree_B <- read.tree(paste0("network_tree_guild_B_",gsub(x=gsub(x=i,pattern="network_",replacement = ""),pattern=".csv",replacement = ""),".tre"))
    
    nb_A <- Ntip(tree_A)
    nb_B <- Ntip(tree_B)
    
    name=gsub(x=gsub(x=i,pattern="network_",replacement = ""),pattern=".csv",replacement = "")
    
    # Change B - Signal on A
    tree_B_noise <- read.tree(paste("simulations/tree_noise_B_",name,"_",N,"bp.tre",sep=""))
    cophe_B <- cophenetic.phylo(tree_B)
    cophe_B_noise <- cophenetic.phylo(tree_B_noise)
    cophe_B_noise <- cophe_B_noise[rownames(cophe_B),colnames(cophe_B)]
    mantel_B <- mantel(as.dist(cophe_B_noise) ~ as.dist(cophe_B),  nperm = 10000, correlation="Pearson")
    
    # Change A - Signal on B
    tree_A_noise <- read.tree(paste("simulations/tree_noise_A_",name,"_",N,"bp.tre",sep=""))
    cophe_A <- cophenetic.phylo(tree_A)
    cophe_A_noise <- cophenetic.phylo(tree_A_noise)
    cophe_A_noise <- cophe_A_noise[rownames(cophe_A),colnames(cophe_A)]
    mantel_A <- mantel(as.dist(cophe_A_noise) ~ as.dist(cophe_A),  nperm = 10000, correlation="Pearson")
    
    results_all <- rbind(results_all, c(param, size, size_tot, nb_A, nb_B, N, mantel_A[1:3], mantel_B[1:3] ))
  }
  
  
  results_all <- data.frame(results_all)
  colnames(results_all) <- c("param","size_name","size", "nb_A", "nb_B", "N", "mantel_cor_A", "pvalue_high_A", "pvalue_low_A", "mantel_cor_B", "pvalue_high_B", "pvalue_low_B")
  
  results_all$param <- as.character(results_all$param)
  results_all$size <- as.character(results_all$size)
  results_all$mantel_cor_A <- as.numeric(as.character(results_all$mantel_cor_A))
  results_all$pvalue_high_A <- as.numeric(as.character(results_all$pvalue_high_A))
  results_all$mantel_cor_B <- as.numeric(as.character(results_all$mantel_cor_B))
  results_all$pvalue_high_B <- as.numeric(as.character(results_all$pvalue_high_B))
  results_all$nb_A <- as.numeric(as.character(results_all$nb_A))
  results_all$nb_B <- as.numeric(as.character(results_all$nb_B))
  
  
  # Plot 
  results_signal <- results_all
  
  results_signal$size_tot <- "<150"
  results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>=150)] <- "150-250"
  results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>250)] <- ">250"
  
  
  # for tree A
  
  ggplot(results_signal,aes(x=param,y=mantel_cor_A,fill=size_tot))+xlab("Parameters")+ylab("Mantel correlation (A)") +labs(title=" ")+scale_fill_manual(values=c("#7dcea0","#f5b041","#e67e22"))+transparent_theme_no + 
    geom_hline(yintercept=0, linetype="dashed", color="#273746")+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    geom_boxplot(alpha=0.9,size=0.75, color="#273746")
  
  ggsave(filename=paste0("boxplot_correlation_trees_Noise_Mantel_tests_Pearson_N_",N,"_size_correlation_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)
  
  
  # for tree B
  
  ggplot(results_signal,aes(x=param,y=mantel_cor_B,fill=size_tot))+xlab("Parameters")+ylab("Mantel correlation (B)") +labs(title=" ")+scale_fill_manual(values=c("#7dcea0","#f5b041","#e67e22"))+transparent_theme_no + 
    geom_hline(yintercept=0, linetype="dashed", color="#273746")+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    geom_boxplot(alpha=0.9,size=0.75, color="#273746")
  
  ggsave(filename=paste0("boxplot_correlation_trees_Noise_Mantel_tests_Pearson_N_",N,"_size_correlation_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)
  
  
}



###### Step 6-A: check monophyly - all simulations ########

rm(list=ls())

setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/script/")


library(coda)
library(ape)
library(fields)
library(Matrix)
library(RPANDA, lib.loc="/users/biodiv/bperez/packages/")
library(bipartite)
library(apTreeshape)
library(ade4)
library(parallel)
library(vegan)
library(phytools)
library(ggplot2)
library(GUniFrac,lib.loc="/users/biodiv/bperez/packages/")
library(parallel)
require(R.utils)


library(treebase, lib.loc="/users/biodiv/bperez/packages/")
library(adephylo, lib.loc="/users/biodiv/bperez/packages/")


dyn.load("fitness.so")
dyn.load("fitness_death.so")
source("ComputeMetricsEtablissement.R")
source("Rcode_comment_new.R") # compile c functions
source("beta.R")

dyn.load("../script/permute.so")
source("../script/function_phylo_signal_network.R")

setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1/")


seed=as.integer(1)

library(ggplot2)
transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))


info = list.files(pattern = paste0('.Rdata'))
info <- info[grep(x=info,pattern = "network_simul_1_")]


# name="simul_1_D"
name="simul_1"

results_all <- c()

for (i in info){
  
  if (length(grep(x=i, pattern = "_neutral_")==1)>0) {param="neutral"}
  if (length(grep(x=i, pattern = "_mutualism_1_")==1)>0) {param="mutualism_1"}
  if (length(grep(x=i, pattern = "_mutualism_2_")==1)>0) {param="mutualism_2"}
  if (length(grep(x=i, pattern = "_mutualism_3_")==1)>0) {param="mutualism_3"}
  if (length(grep(x=i, pattern = "_mutualism_4_")==1)>0) {param="mutualism_4"}
  if (length(grep(x=i, pattern = "_mutualism_5_")==1)>0) {param="mutualism_5"}
  if (length(grep(x=i, pattern = "_mutualism_6_")==1)>0) {param="mutualism_6"}
  if (length(grep(x=i, pattern = "_antagonism_1_")==1)>0) {param="antagonism_1"}
  if (length(grep(x=i, pattern = "_antagonism_2_")==1)>0) {param="antagonism_2"}
  if (length(grep(x=i, pattern = "_antagonism_3_")==1)>0) {param="antagonism_3"}
  if (length(grep(x=i, pattern = "_antagonism_4_")==1)>0) {param="antagonism_4"}
  if (length(grep(x=i, pattern = "_antagonism_5_")==1)>0) {param="antagonism_5"}
  if (length(grep(x=i, pattern = "_antagonism_6_")==1)>0) {param="antagonism_6"}
  if (length(grep(x=i, pattern = "_antagonism_7_")==1)>0) {param="antagonism_7"}
  if (length(grep(x=i, pattern = "_antagonism_8_")==1)>0) {param="antagonism_8"}
  if (length(grep(x=i, pattern = "_antagonism_9_")==1)>0) {param="antagonism_9"}
  
  if (length(grep(x=i, pattern = "_simul_1_A_")==1)>0) {name="simul_1_A"}
  if (length(grep(x=i, pattern = "_simul_1_B_")==1)>0) {name="simul_1_B"}
  if (length(grep(x=i, pattern = "_simul_1_C_")==1)>0) {name="simul_1_C"}
  if (length(grep(x=i, pattern = "_simul_1_D_")==1)>0) {name="simul_1_D"}
  if (length(grep(x=i, pattern = "_simul_1_E_")==1)>0) {name="simul_1_E"}
  if (length(grep(x=i, pattern = "_simul_1_F_")==1)>0) {name="simul_1_F"}
  
  
  seed=strsplit(i, split="_")[[1]][length(strsplit(i, split="_")[[1]])]
  seed=strsplit(seed,split=".",fixed = T)[[1]][1]
  
  set.seed(seed)
  load(file = paste0("network_",name,"_", param,"_seed_",seed,".Rdata"))
  
  tree_A <- read.tree(paste0("network_tree_guild_A_",name,"_", param,"_seed_",seed,".tre"))
  tree_B <- read.tree(paste0("network_tree_guild_B_",name,"_", param,"_seed_",seed,".tre"))
  
  
  phylo <- define.species(gen,threshold=1, monophyly = FALSE, seed=5)
  
  
  # tree A
  cophe_A <- cophenetic(tree_A)
  cophe_A_2 <- cophenetic(phylo$Hphylo$tree)
  cophe_A_2 <- cophe_A_2[rownames(cophe_A),colnames(cophe_A)]
  mantel_A <- RPANDA::mantel_test(as.dist(cophe_A) ~ as.dist(cophe_A_2), correlation="Pearson", nperm = 10000)
  
  # tree B
  cophe_B <- cophenetic(tree_B)
  cophe_B_2 <- cophenetic(phylo$Pphylo$tree)
  cophe_B_2 <- cophe_B_2[rownames(cophe_B),colnames(cophe_B)]
  mantel_B <- RPANDA::mantel_test(as.dist(cophe_B) ~ as.dist(cophe_B_2), correlation="Pearson", nperm = 10000)
  
  
  
  results_all <- rbind(results_all, c(name, param, Ntip(tree_A), Ntip(tree_B), mantel_A[1:3], mantel_B[1:3] ))
}

name="simul_1"

results_all <- data.frame(results_all)
colnames(results_all) <- c("name","param", "nb_A", "nb_B", "mantel_cor_A", "pvalue_high_A", "pvalue_low_A", "mantel_cor_B", "pvalue_high_B", "pvalue_low_B")

results_all$param <- as.character(results_all$param)
results_all$mantel_cor_A <- as.numeric(as.character(results_all$mantel_cor_A))
results_all$pvalue_high_A <- as.numeric(as.character(results_all$pvalue_high_A))
results_all$mantel_cor_B <- as.numeric(as.character(results_all$mantel_cor_B))
results_all$pvalue_high_B <- as.numeric(as.character(results_all$pvalue_high_B))
results_all$nb_A <- as.numeric(as.character(results_all$nb_A))
results_all$nb_B <- as.numeric(as.character(results_all$nb_B))

save(list = ls(), file=paste0("test_trees_", name,"_monophyly.Rdata"))


results_signal <- results_all

results_signal$size_tot <- "<150"
results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>=150)] <- "150-250"
results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>250)] <- ">250"


results_signal$size_tot <- "<70"
results_signal$size_tot[which(results_signal$nb_A>=70)] <- "70-100"
results_signal$size_tot[which(results_signal$nb_A>100)] <- ">100"


# for tree A

ggplot(results_signal,aes(x=param,y=mantel_cor_A,fill=size_tot))+xlab("Parameters")+ylab("Mantel correlation (A)") +labs(title=" ")+scale_fill_manual(values=c("#7dcea0","#f5b041","#e67e22"))+transparent_theme_no + 
  geom_hline(yintercept=0, linetype="dashed", color="#273746")+
  scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                              bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
  geom_boxplot(alpha=0.9,size=0.75, color="#273746")+
  facet_grid(~size_tot)

ggsave(filename=paste0("boxplot_correlation_trees_simul_1_Mantel_tests_Pearson_size_correlation_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 11, height = 7)


# for tree B


results_signal$size_tot <- "<70"
results_signal$size_tot[which(results_signal$nb_B>=70)] <- "70-100"
results_signal$size_tot[which(results_signal$nb_B>100)] <- ">100"


ggplot(results_signal,aes(x=param,y=mantel_cor_B,fill=size_tot))+xlab("Parameters")+ylab("Mantel correlation (B)") +labs(title=" ")+scale_fill_manual(values=c("#7dcea0","#f5b041","#e67e22"))+transparent_theme_no + 
  geom_hline(yintercept=0, linetype="dashed", color="#273746")+
  scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                              bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
  geom_boxplot(alpha=0.9,size=0.75, color="#273746")+
  facet_grid(~size_tot)

ggsave(filename=paste0("boxplot_correlation_trees_simul_1_Mantel_tests_Pearson_size_correlation_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 11, height = 7)



#####  Make final plot

rm(list=ls())

setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1/")


seed=as.integer(1)

library(ggplot2)
library(dplyr)
transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))

name="simul_1"
load(file=paste0("test_trees_", name,"_monophyly.Rdata"))

results_signal <- results_all

results_signal$size_tot <- "<150"
results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>=150)] <- "150-250"
results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>250)] <- ">250"


results_signal$size_tot <- "<70"
results_signal$size_tot[which(results_signal$nb_A>=70)] <- "70-100"
results_signal$size_tot[which(results_signal$nb_A>100)] <- ">100"

results_signal$size_tot <- as.factor(results_signal$size_tot)
results_signal$size_tot <-  factor(results_signal$size_tot,levels(results_signal$size_tot)[c(2,3,1)])


# for tree A

ggplot(results_signal,aes(x=param,y=mantel_cor_A,fill=size_tot))+xlab("Parameters")+ylab("Mantel correlation (A)") +labs(title=" ")+scale_fill_manual(values=c("#7dcea0","#f5b041","#e67e22"))+transparent_theme_no + 
  geom_hline(yintercept=0, linetype="dashed", color="#273746")+
  scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                              bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
  geom_boxplot(alpha=0.9,size=0.75, color="#273746")+
  facet_grid(~size_tot)

ggsave(filename=paste0("boxplot_correlation_trees_simul_1_Mantel_tests_Pearson_size_correlation_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 11, height = 7)


# for tree B


results_signal$size_tot <- "<70"
results_signal$size_tot[which(results_signal$nb_B>=70)] <- "70-100"
results_signal$size_tot[which(results_signal$nb_B>100)] <- ">100"

results_signal$size_tot <- as.factor(results_signal$size_tot)
results_signal$size_tot <-  factor(results_signal$size_tot,levels(results_signal$size_tot)[c(2,3,1)])

ggplot(results_signal,aes(x=param,y=mantel_cor_B,fill=size_tot))+xlab("Parameters")+ylab("Mantel correlation (B)") +labs(title=" ")+scale_fill_manual(values=c("#7dcea0","#f5b041","#e67e22"))+transparent_theme_no + 
  geom_hline(yintercept=0, linetype="dashed", color="#273746")+
  scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                              bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
  geom_boxplot(alpha=0.9,size=0.75, color="#273746")+
  facet_grid(~size_tot)

ggsave(filename=paste0("boxplot_correlation_trees_simul_1_Mantel_tests_Pearson_size_correlation_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 11, height = 7)



###### Step 6-B: check the signal in the traits - all simulations - Mantel tests ########

rm(list=ls())

setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/script/")


library(coda)
library(ape)
library(fields)
library(Matrix)
library(RPANDA, lib.loc="/users/biodiv/bperez/packages/")
library(bipartite)
library(apTreeshape)
library(ade4)
library(parallel)
library(vegan)
library(phytools)
library(ggplot2)
library(GUniFrac,lib.loc="/users/biodiv/bperez/packages/")
library(parallel)
require(R.utils)
library(dplyr)

library(treebase, lib.loc="/users/biodiv/bperez/packages/")
library(adephylo, lib.loc="/users/biodiv/bperez/packages/")


dyn.load("fitness.so")
dyn.load("fitness_death.so")
source("ComputeMetricsEtablissement.R")
source("Rcode_comment_new.R") # compile c functions
source("beta.R")

dyn.load("../script/permute.so")
source("../script/function_phylo_signal_network.R")

setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1/")


seed=as.integer(1)

library(ggplot2)
transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))

results_signal_all_names <- c()

name=c("_D_")
for (name in c("_A_","_B_","_C_","_D_","_E_","_F_" )){
  
  print(name)
  
  info = list.files(pattern = paste0('.Rdata'))
  info <- info[grep(x=info,pattern = name)]
  info <- info[grep(x=info,pattern = "_seed_")]
  
  if (name=="_A_") name="simul_1_A"
  if (name=="_B_") name="simul_1_B"
  if (name=="_C_") name="simul_1_C"
  if (name=="_D_") name="simul_1_D"
  if (name=="_E_") name="simul_1_E"
  if (name=="_F_") name="simul_1_F"
  
  results_all <- c()
  
  for (i in info){
    
    if (length(grep(x=i, pattern = "_neutral_")==1)>0) {param="neutral"}
    if (length(grep(x=i, pattern = "_mutualism_1_")==1)>0) {param="mutualism_1"}
    if (length(grep(x=i, pattern = "_mutualism_2_")==1)>0) {param="mutualism_2"}
    if (length(grep(x=i, pattern = "_mutualism_3_")==1)>0) {param="mutualism_3"}
    if (length(grep(x=i, pattern = "_mutualism_4_")==1)>0) {param="mutualism_4"}
    if (length(grep(x=i, pattern = "_mutualism_5_")==1)>0) {param="mutualism_5"}
    if (length(grep(x=i, pattern = "_mutualism_6_")==1)>0) {param="mutualism_6"}
    if (length(grep(x=i, pattern = "_antagonism_1_")==1)>0) {param="antagonism_1"}
    if (length(grep(x=i, pattern = "_antagonism_2_")==1)>0) {param="antagonism_2"}
    if (length(grep(x=i, pattern = "_antagonism_3_")==1)>0) {param="antagonism_3"}
    if (length(grep(x=i, pattern = "_antagonism_4_")==1)>0) {param="antagonism_4"}
    if (length(grep(x=i, pattern = "_antagonism_5_")==1)>0) {param="antagonism_5"}
    if (length(grep(x=i, pattern = "_antagonism_6_")==1)>0) {param="antagonism_6"}
    if (length(grep(x=i, pattern = "_antagonism_7_")==1)>0) {param="antagonism_7"}
    if (length(grep(x=i, pattern = "_antagonism_8_")==1)>0) {param="antagonism_8"}
    if (length(grep(x=i, pattern = "_antagonism_9_")==1)>0) {param="antagonism_9"}
    
    seed=strsplit(i, split="_")[[1]][length(strsplit(i, split="_")[[1]])]
    seed=strsplit(seed,split=".",fixed = T)[[1]][1]
    
    set.seed(seed)
    load(file = paste0("network_",name,"_", param,"_seed_",seed,".Rdata"))
    
    tree_A <- read.tree(paste0("network_tree_guild_A_",name,"_", param,"_seed_",seed,".tre"))
    tree_B <- read.tree(paste0("network_tree_guild_B_",name,"_", param,"_seed_",seed,".tre"))
    
    cophe_A <- cophenetic(tree_A)
    cophe_B <- cophenetic(tree_B)
    
    if (name=="simul_1_D") {phy1 <- define.species(gen,threshold=1, monophyly = FALSE, seed=1)}
    
    # tree A
    traits_A <- phy1$Hphylo$mean.trait 
    colnames(traits_A) <- as.character(phy1$Hphylo$tree$tip.label)
    dist_traits_A <- dist(t(traits_A))
    dist_traits_A <- as.matrix(dist_traits_A)[rownames(cophe_A),colnames(cophe_A)]
    if (all(dist_traits_A==0)){mantel_A <- c("same_traits","NA","NA")
    }else{mantel_A <- RPANDA::mantel_test(as.dist(cophe_A) ~ as.dist(dist_traits_A), correlation="Pearson", nperm = 10000)}
    
    # tree B
    traits_B <- phy1$Pphylo$mean.trait 
    colnames(traits_B) <- as.character(phy1$Pphylo$tree$tip.label) 
    dist_traits_B <- dist(t(traits_B))
    dist_traits_B <- as.matrix(dist_traits_B)[rownames(cophe_B),colnames(cophe_B)]
    if (all(dist_traits_B==0)){mantel_B <- c("same_traits","NA","NA")
    }else{mantel_B <- RPANDA::mantel_test(as.dist(cophe_B) ~ as.dist(dist_traits_B), correlation="Pearson", nperm = 10000)}
    
    
    results_all <- rbind(results_all, c(name, param, seed, Ntip(tree_A), Ntip(tree_B), mantel_A[1:3], mantel_B[1:3] ))
  }
  
  
  results_all <- data.frame(results_all)
  colnames(results_all) <- c("name","param", "seed", "nb_A", "nb_B", "mantel_cor_A", "pvalue_high_A", "pvalue_low_A", "mantel_cor_B", "pvalue_high_B", "pvalue_low_B")
  
  # remove lines with no trait variation:
  results_all <- results_all[which(results_all$mantel_cor_A!="same_traits"),]
  results_all <- results_all[which(results_all$mantel_cor_B!="same_traits"),]
  
  # Plot 
  save(list = ls(), file=paste0("test_traits_", name,".Rdata"))
  
  results_all$param <- as.character(results_all$param)
  results_all$mantel_cor_A <- as.numeric(as.character(results_all$mantel_cor_A))
  results_all$pvalue_high_A <- as.numeric(as.character(results_all$pvalue_high_A))
  results_all$mantel_cor_B <- as.numeric(as.character(results_all$mantel_cor_B))
  results_all$pvalue_high_B <- as.numeric(as.character(results_all$pvalue_high_B))
  results_all$nb_A <- as.numeric(as.character(results_all$nb_A))
  results_all$nb_B <- as.numeric(as.character(results_all$nb_B))
  
  
  results_signal <- results_all
  
  print(table(results_signal$param))
  
  results_signal_all_names <- rbind(results_signal_all_names, results_signal)
  
  write.table(results_signal, paste0("results_phylogenetic_signal_simul_1_in_traits_",name,".csv"), sep=";", quote=F,row.names = F)
  
  
}

write.table(results_signal_all_names, paste0("results_phylogenetic_signal_simul_1_in_traits.csv"), sep=";", quote=F,row.names = F)


#####  Make final plot

rm(list=ls())

setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1/")


seed=as.integer(1)

library(ggplot2)
library(dplyr)
transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))


results_all <- read.table(paste0("results_phylogenetic_signal_simul_1_in_traits.csv"), sep=";", header=T)

results_all$param <- as.character(results_all$param)
results_all$mantel_cor_A <- as.numeric(as.character(results_all$mantel_cor_A))
results_all$pvalue_high_A <- as.numeric(as.character(results_all$pvalue_high_A))
results_all$mantel_cor_B <- as.numeric(as.character(results_all$mantel_cor_B))
results_all$pvalue_high_B <- as.numeric(as.character(results_all$pvalue_high_B))
results_all$nb_A <- as.numeric(as.character(results_all$nb_A))
results_all$nb_B <- as.numeric(as.character(results_all$nb_B))




# for tree A

results_signal <- results_all

results_signal$size_tot <- "<70"
results_signal$size_tot[which(results_signal$nb_A>=70)] <- "70-100"
results_signal$size_tot[which(results_signal$nb_A>100)] <- ">100"
results_signal$size_tot <- as.factor(results_signal$size_tot)

results_signal$Inference <- "d) Not significant signal"
results_signal$Inference[which(results_signal$pvalue_high_A<0.05)] <- "e) Significant signal"
results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.05))] <- "f) Significant signal (R>0.05)"
results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.15))] <- "g) Significant signal (R>0.15)"
results_signal$Inference[which(results_signal$pvalue_low_A<0.05)] <- "c) Significant anti-signal"
results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.05)))] <- "b) Significant anti-signal (R<-0.05)"
results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.15)))] <- "a) Significant anti-signal (R<-0.15)"

colors <- c()
if ("a) Significant anti-signal (R<-0.15)" %in% results_signal$Inference) {colors <- c(colors, "#943126")}
if ("b) Significant anti-signal (R<-0.05)" %in% results_signal$Inference) {colors <- c(colors, "#e74c3c")}
if ("c) Significant anti-signal" %in% results_signal$Inference) {colors <- c(colors, "#f1948a")}
if ("d) Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
if ("e) Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
if ("f) Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
if ("g) Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}

brks <- c(0, 0.25, 0.5, 0.75, 1)

results_signal_all <- results_signal %>% 
  group_by(param,Inference,size_tot) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))

results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
results_signal_all$size_tot <-  factor(results_signal_all$size_tot,levels(results_signal_all$size_tot)[c(2,3,1)])

ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
  xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
  scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
  scale_color_manual(values=colors)+
  scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                              bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
  facet_grid(~size_tot)

ggsave(filename=paste0("boxplot_correlation_traits_phylo_size_inference_Mantel_test_anticorrelation_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 11, height = 7)
# remove "device=cairo_pdf," on 06/01/2021


# for tree B

results_signal <- results_all

results_signal$size_tot <- "<70"
results_signal$size_tot[which(results_signal$nb_B>=70)] <- "70-100"
results_signal$size_tot[which(results_signal$nb_B>100)] <- ">100"
results_signal$size_tot <- as.factor(results_signal$size_tot)

results_signal$Inference <- "d) Not significant signal"
results_signal$Inference[which(results_signal$pvalue_high_B<0.05)] <- "e) Significant signal"
results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.05))] <- "f) Significant signal (R>0.05)"
results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.15))] <- "g) Significant signal (R>0.15)"
results_signal$Inference[which(results_signal$pvalue_low_B<0.05)] <- "c) Significant anti-signal"
results_signal$Inference[intersect(which(results_signal$pvalue_low_B<0.05), which(results_signal$mantel_cor_B<(-0.05)))] <- "b) Significant anti-signal (R<-0.05)"
results_signal$Inference[intersect(which(results_signal$pvalue_low_B<0.05), which(results_signal$mantel_cor_B<(-0.15)))] <- "a) Significant anti-signal (R<-0.15)"

colors <- c()
if ("a) Significant anti-signal (R<-0.15)" %in% results_signal$Inference) {colors <- c(colors, "#943126")}
if ("b) Significant anti-signal (R<-0.05)" %in% results_signal$Inference) {colors <- c(colors, "#e74c3c")}
if ("c) Significant anti-signal" %in% results_signal$Inference) {colors <- c(colors, "#f1948a")}
if ("d) Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
if ("e) Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
if ("f) Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
if ("g) Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}

brks <- c(0, 0.25, 0.5, 0.75, 1)

results_signal_all <- results_signal %>% 
  group_by(param,Inference,size_tot) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))

results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
results_signal_all$size_tot <-  factor(results_signal_all$size_tot,levels(results_signal_all$size_tot)[c(2,3,1)])

ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
  xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
  scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
  scale_color_manual(values=colors)+
  scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                              bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
  facet_grid(~size_tot)

ggsave(filename=paste0("boxplot_correlation_traits_phylo_size_inference_Mantel_test_anticorrelation_clade_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 11, height = 7)




###### Step 6-C: check the signal in the traits - all simulations - Pagel lambda ########


rm(list=ls())

setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/script/")


library(ape)
library(phytools)
library(ggplot2)
library(RPANDA, lib.loc="/users/biodiv/bperez/packages/")
library(GUniFrac,lib.loc="/users/biodiv/bperez/packages/")
library(parallel)
library(mvMORPH)
library(dplyr)

library(phylocurve)

source("../script/ComputeMetricsEtablissement.R")
source("../script/Rcode_comment_new.R") # compile c functions
source("../script/beta.R")

source("../script/function_phylo_signal_network.R")

setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1/")


seed=as.integer(1)

transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))

results_signal_all_names <- c()

name=c("_A_")
name=c("_D_")
for (name in c("_A_","_B_","_C_","_D_","_E_","_F_" )){
  
  print(name)
  
  info = list.files(pattern = paste0('.Rdata'))
  info <- info[grep(x=info,pattern = name)]
  info <- info[grep(x=info,pattern = "_seed_")]
  
  if (name=="_A_") name="simul_1_A"
  if (name=="_B_") name="simul_1_B"
  if (name=="_C_") name="simul_1_C"
  if (name=="_D_") name="simul_1_D"
  if (name=="_E_") name="simul_1_E"
  if (name=="_F_") name="simul_1_F"
  
  results_all <- c()
  
  for (i in info){
    
    if (length(grep(x=i, pattern = "_neutral_")==1)>0) {param="neutral"}
    if (length(grep(x=i, pattern = "_mutualism_1_")==1)>0) {param="mutualism_1"}
    if (length(grep(x=i, pattern = "_mutualism_2_")==1)>0) {param="mutualism_2"}
    if (length(grep(x=i, pattern = "_mutualism_3_")==1)>0) {param="mutualism_3"}
    if (length(grep(x=i, pattern = "_mutualism_4_")==1)>0) {param="mutualism_4"}
    if (length(grep(x=i, pattern = "_mutualism_5_")==1)>0) {param="mutualism_5"}
    if (length(grep(x=i, pattern = "_mutualism_6_")==1)>0) {param="mutualism_6"}
    if (length(grep(x=i, pattern = "_antagonism_1_")==1)>0) {param="antagonism_1"}
    if (length(grep(x=i, pattern = "_antagonism_2_")==1)>0) {param="antagonism_2"}
    if (length(grep(x=i, pattern = "_antagonism_3_")==1)>0) {param="antagonism_3"}
    if (length(grep(x=i, pattern = "_antagonism_4_")==1)>0) {param="antagonism_4"}
    if (length(grep(x=i, pattern = "_antagonism_5_")==1)>0) {param="antagonism_5"}
    if (length(grep(x=i, pattern = "_antagonism_6_")==1)>0) {param="antagonism_6"}
    if (length(grep(x=i, pattern = "_antagonism_7_")==1)>0) {param="antagonism_7"}
    if (length(grep(x=i, pattern = "_antagonism_8_")==1)>0) {param="antagonism_8"}
    if (length(grep(x=i, pattern = "_antagonism_9_")==1)>0) {param="antagonism_9"}
    
    seed=strsplit(i, split="_")[[1]][length(strsplit(i, split="_")[[1]])]
    seed=strsplit(seed,split=".",fixed = T)[[1]][1]
    
    set.seed(seed)
    load(file = paste0("network_",name,"_", param,"_seed_",seed,".Rdata"))
    
    tree_A <- read.tree(paste0("network_tree_guild_A_",name,"_", param,"_seed_",seed,".tre"))
    tree_B <- read.tree(paste0("network_tree_guild_B_",name,"_", param,"_seed_",seed,".tre"))
    
    cophe_A <- cophenetic(tree_A)
    cophe_B <- cophenetic(tree_B)
    
    
    if (name=="simul_1_D") {phy1 <- define.species(gen,threshold=1, monophyly = FALSE, seed=1)}
    
    
    # tree A
    traits_A <- phy1$Hphylo$mean.trait 
    colnames(traits_A) <- as.character(phy1$Hphylo$tree$tip.label)
    dist_traits_A <- dist(t(traits_A))
    dist_traits_A <- as.matrix(dist_traits_A)[rownames(cophe_A),colnames(cophe_A)]
    if (all(dist_traits_A==0)){test_A <- c("same_traits","NA","NA")
    }else{
      fit <- evo.model(tree = tree_A, Y = t(traits_A), model="lambda", method = "Pairwise ML")
      fit0 <- evo.model(tree = geiger::rescale(tree_A, model="lambda", lambda=0), Y = t(traits_A), model="lambda", method = "Pairwise ML")
      LR <- 2*(fit$logL-fit0$logL)
      p_LRT <- pchisq(LR,df=1,lower.tail=FALSE)
      test_A <- c(fit$model.par,LR,p_LRT)
    }
    
    # tree B
    traits_B <- phy1$Pphylo$mean.trait 
    colnames(traits_B) <- as.character(phy1$Pphylo$tree$tip.label) 
    dist_traits_B <- dist(t(traits_B))
    dist_traits_B <- as.matrix(dist_traits_B)[rownames(cophe_B),colnames(cophe_B)]
    if (all(dist_traits_B==0)){test_B <- c("same_traits","NA","NA")
    }else{
      fit <- evo.model(tree = tree_B, Y = t(traits_B), model="lambda", method = "Pairwise ML")
      fit0 <- evo.model(tree = geiger::rescale(tree_B, model="lambda", lambda=0), Y = t(traits_B), model="lambda", method = "Pairwise ML")
      LR <- 2*(fit$logL-fit0$logL)
      p_LRT <- pchisq(LR,df=1,lower.tail=FALSE)
      test_B <- c(fit$model.par,LR,p_LRT)
    }
    
    results_all <- rbind(results_all, c(name, param, seed, Ntip(tree_A), Ntip(tree_B), test_A, test_B ))
  }
  
  
  results_all <- data.frame(results_all)
  colnames(results_all) <- c("name","param", "seed", "nb_A", "nb_B", "lambda_A", "LR_A", "p_LRT_A", "lambda_B", "LR_B", "p_LRT_B")
  
  # remove lines with no trait variation:
  results_all <- results_all[which(results_all$lambda_A!="same_traits"),]
  results_all <- results_all[which(results_all$lambda_B!="same_traits"),]
  
  # Plot 
  save(list = ls(), file=paste0("test_traits_Pagel_lambda_", name,".Rdata"))
  
  
  results_all$param <- as.character(results_all$param)
  results_all$lambda_A <- as.numeric(as.character(results_all$lambda_A))
  results_all$p_LRT_A <- as.numeric(as.character(results_all$p_LRT_A))
  results_all$lambda_B <- as.numeric(as.character(results_all$lambda_B))
  results_all$p_LRT_B <- as.numeric(as.character(results_all$p_LRT_B))
  results_all$nb_A <- as.numeric(as.character(results_all$nb_A))
  results_all$nb_B <- as.numeric(as.character(results_all$nb_B))
  
  
  results_signal <- results_all
  
  print(table(results_signal$param))
  
  results_signal_all_names <- rbind(results_signal_all_names, results_signal)
  
  write.table(results_signal, paste0("results_phylogenetic_signal_simul_1_in_traits_Pagel_lambda_",name,".csv"), sep=";", quote=F,row.names = F)
  
  
}

write.table(results_signal_all_names, paste0("results_phylogenetic_signal_simul_1_in_traits_Pagel_lambda.csv"), sep=";", quote=F,row.names = F)


# scp bperez@jord.biologie.ens.fr:/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1/boxplot_correlation_traits_phylo_*.pdf  /Users/bperez/ownCloud/Recherche/These/ISYEB/signal_phylo/simul_1/


#####  Make final plot

rm(list=ls())

setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1/")


seed=as.integer(1)

library(ggplot2)
library(dplyr)
transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))


results_all <- read.table(paste0("results_phylogenetic_signal_simul_1_in_traits_Pagel_lambda.csv"), sep=";", header=T)

results_all$param <- as.character(results_all$param)
results_all$lambda_A <- as.numeric(as.character(results_all$lambda_A))
results_all$pvalue_high_A <- results_all$p_LRT_A <- as.numeric(as.character(results_all$p_LRT_A))
results_all$lambda_B <- as.numeric(as.character(results_all$lambda_B))
results_all$pvalue_high_B <- results_all$p_LRT_B <- as.numeric(as.character(results_all$p_LRT_B))
results_all$nb_A <- as.numeric(as.character(results_all$nb_A))
results_all$nb_B <- as.numeric(as.character(results_all$nb_B))




# for tree A

results_signal <- results_all

results_signal$size_tot <- "<70"
results_signal$size_tot[which(results_signal$nb_A>=70)] <- "70-100"
results_signal$size_tot[which(results_signal$nb_A>100)] <- ">100"
results_signal$size_tot <- as.factor(results_signal$size_tot)

results_signal$Inference <- "a) Not significant signal"
results_signal$Inference[which(results_signal$pvalue_high_A<0.05)] <- "b) Significant signal"
results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$lambda_A>0.5))] <- "c) Significant signal (lambda>0.5)"
results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$lambda_A>0.8))] <- "d) Significant signal (lambda>0.8)"

colors <- c()
if ("a) Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
if ("b) Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
if ("c) Significant signal (lambda>0.5)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
if ("d) Significant signal (lambda>0.8)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}

brks <- c(0, 0.25, 0.5, 0.75, 1)

results_signal_all <- results_signal %>% 
  group_by(param,Inference,size_tot) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))

results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
results_signal_all$size_tot <-  factor(results_signal_all$size_tot,levels(results_signal_all$size_tot)[c(2,3,1)])

ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
  xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
  scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
  scale_color_manual(values=colors)+
  scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                              bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
  facet_grid(~size_tot)

ggsave(filename=paste0("boxplot_correlation_traits_phylo_size_inference_Pagel_lambda_anticorrelation_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 11, height = 7)

# for tree B

results_signal <- results_all

results_signal$size_tot <- "<70"
results_signal$size_tot[which(results_signal$nb_B>=70)] <- "70-100"
results_signal$size_tot[which(results_signal$nb_B>100)] <- ">100"
results_signal$size_tot <- as.factor(results_signal$size_tot)

results_signal$Inference <- "a) Not significant signal"
results_signal$Inference[which(results_signal$pvalue_high_B<0.05)] <- "b) Significant signal"
results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$lambda_B>0.5))] <- "c) Significant signal (lambda>0.5)"
results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$lambda_B>0.8))] <- "d) Significant signal (lambda>0.8)"

colors <- c()
if ("a) Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
if ("b) Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
if ("c) Significant signal (lambda>0.5)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
if ("d) Significant signal (lambda>0.8)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}

brks <- c(0, 0.25, 0.5, 0.75, 1)

results_signal_all <- results_signal %>% 
  group_by(param,Inference,size_tot) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))

results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
results_signal_all$size_tot <-  factor(results_signal_all$size_tot,levels(results_signal_all$size_tot)[c(2,3,1)])

ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
  xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
  scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
  scale_color_manual(values=colors)+
  scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                              bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
  facet_grid(~size_tot)

ggsave(filename=paste0("boxplot_correlation_traits_phylo_size_inference_Pagel_lambda_anticorrelation_clade_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 11, height = 7)




###### Step 7-A: Compute clade-specific phylogenetic signal ########

rm(list=ls())

library(ape)

setwd("/Users/bperez/ownCloud/Recherche/These/ISYEB/signal_phylo/simul_1/")


dyn.load("../script/permute.so")
source("../script/function_phylo_signal_network.R")


name="simul_1_A_antagonism_1_seed_5"

network <- read.table(paste0("network_",name,".csv"), header=T,sep=";")
colnames(network) <- seq(1:ncol(network))
rownames(network) <- seq(1:nrow(network))

tree_A <- read.tree(paste0("network_tree_guild_A_",name,".tre"))
tree_B <- read.tree(paste0("network_tree_guild_B_",name,".tre"))


network <- network[tree_B$tip.label,tree_A$tip.label]


# for tree_A
host_tree <- tree_A
symbiont_tree <- tree_B
network # host in columns and symbionts in rows
minimum=10
method="GUniFrac"
correlation="Pearson"
nperm=100000


results_sub_clades <- phylosignal_sub_network(network, host_tree, symbiont_tree, method = "GUniFrac", nperm = 100000, correlation = "Pearson", minimum=10)


pdf("test_simul_1_A_antagonism_1_seed_5_sub_clades.pdf")
plot_phylosignal_sub_network(host_tree, results_sub_clades, legend=TRUE, show.tip.label=FALSE, where="bottomleft")
dev.off()



### mutualism 5 on neutral #####

rm(list=ls()) # on cluster

library(ggplot2)
library(vegan)
library(phytools)
library(ggplot2)

library(GUniFrac,lib.loc="/users/biodiv/bperez/packages/")
library(parallel)
require(R.utils)
library(RPANDA,lib.loc="/users/biodiv/bperez/packages/")

#setwd("/Users/bperez/ownCloud/Recherche/These/ISYEB/signal_phylo/simul_1/")
setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1")

dyn.load("../script/permute.so")
source("../script/function_phylo_signal_network.R")


size="A"
seed=1
seed=14
size="B"
seed=8

results <- c()

for (size in c("A","B","C","D","E","F")){
  
  seed_parrallel <- function(seed){
    
    # mutualism
    network <- read.table(paste0("network_simul_1_",size,"_mutualism_5_seed_",seed,".csv"), header=T,sep=";")
    colnames(network) <- seq(1:ncol(network))
    rownames(network) <- seq(1:nrow(network))
    tree_A <- read.tree(paste0("network_tree_guild_A_simul_1_",size,"_mutualism_5_seed_",seed,".tre"))
    tree_B <- read.tree(paste0("network_tree_guild_B_simul_1_",size,"_mutualism_5_seed_",seed,".tre"))
    network <- network[tree_B$tip.label,tree_A$tip.label]
    rownames(network) <- paste0("mut_", rownames(network))
    colnames(network) <- paste0("mut_", colnames(network))
    tree_A$tip.label <- paste0("mut_", tree_A$tip.label)
    tree_B$tip.label <- paste0("mut_", tree_B$tip.label)
    
    # neutral 
    network_large <- read.table(paste0("network_simul_1_",size,"_neutral_seed_",seed,".csv"), header=T,sep=";")
    colnames(network_large) <- seq(1:ncol(network_large))
    rownames(network_large) <- seq(1:nrow(network_large))
    tree_large_A <- read.tree(paste0("network_tree_guild_A_simul_1_",size,"_neutral_seed_",seed,".tre"))
    tree_large_B <- read.tree(paste0("network_tree_guild_B_simul_1_",size,"_neutral_seed_",seed,".tre"))
    network_large <- network_large[tree_large_B$tip.label,tree_large_A$tip.label]
    
    
    # Combine network
    network_combined <- rbind(cbind(network_large, matrix(0, nrow = nrow(network_large), ncol = ncol(network), dimnames = list(rownames(network_large), colnames(network)))), 
                              cbind(matrix(0, nrow = nrow(network), ncol = ncol(network_large), dimnames = list(rownames(network), colnames(network_large)) ), network))
    
    # Combine tree A
    set.seed(seed)
    node <- sample(size=1, (Ntip(tree_large_A)+1):(2*Ntip(tree_large_A)-1)) # pick node at random
    tree_A$root.edge <- tree_large_A$edge.length[which(tree_large_A$edge[,2]==node)]/2
    scale <- max(node.depth.edgelength(tree_large_A)) - node.depth.edgelength(tree_large_A)[node] #- tree_A$root.edge
    tree_A$edge.length <- tree_A$edge.length/max(node.depth.edgelength(tree_A))*scale
    tree_A_combined <- bind.tree(x=tree_large_A, y=tree_A, where=node, position = tree_large_A$edge.length[which(tree_large_A$edge[,2]==node)]/2 )
    
    # Combine tree B
    set.seed(seed)
    node <- sample(size=1, (Ntip(tree_large_B)+1):(2*Ntip(tree_large_B)-1)) # pick node at random that can not be the root
    while (node==(Ntip(tree_large_B)+1)){
      node <- sample(size=1, (Ntip(tree_large_B)+1):(2*Ntip(tree_large_B)-1)) 
    }
    tree_B$root.edge <- tree_large_B$edge.length[which(tree_large_B$edge[,2]==node)]/2
    scale <- max(node.depth.edgelength(tree_large_B)) - node.depth.edgelength(tree_large_B)[node] #- tree_A$root.edge
    tree_B$edge.length <- tree_B$edge.length/max(node.depth.edgelength(tree_B))*scale
    tree_B_combined <- bind.tree(x=tree_large_B, y=tree_B, where=node, position = tree_large_B$edge.length[which(tree_large_B$edge[,2]==node)]/2 )
    
    # test sub-clade A
    results_sub_clades <- phylosignal_sub_network(network_combined, tree_A_combined, tree_B_combined, method = "GUniFrac", nperm = 100000, correlation = "Pearson", minimum=10)
    pdf(paste0("sub_clades/test_Mantel_sub_clade_simul_1_guild_A_simul_1_",size,"_mutualism_5_seed_",seed,".pdf"))
    plot_phylosignal_sub_network(tree_A_combined, results_sub_clades, legend=TRUE, show.tip.label=FALSE, where="bottomleft")
    dev.off()
    MRCA_mut <- getMRCA(tree_A_combined, tip = tree_A_combined$tip.label[grep("mut_", tree_A_combined$tip.label)] )
    if (results_sub_clades$pvalue_high_corrected[which(results_sub_clades$node==MRCA_mut)]<=0.05) {work_A <- 1} else {work_A <- 0}
    
    # number of ascending nodes and independent nodes with signal 
    ascending_nodes <- c()
    indep_nodes <- c()
    for (i in (Ntip(tree_A_combined)+1):(Ntip(tree_A_combined)+Nnode(tree_A_combined))){
      extracted_tree <- extract.clade(tree_A_combined,  node=i)
      if (all(tree_A$tip.label %in% extracted_tree$tip.label)){
        if (length(tree_A$tip.label<length(extracted_tree$tip.label))){
          ascending_nodes <- c(ascending_nodes, i)
        }
      }
      if (length(which(tree_A$tip.label %in% extracted_tree$tip.label))==0){
        indep_nodes <- c(indep_nodes, i)
      }
    }
    
    work_A_ascending <- length(which(results_sub_clades$pvalue_high_corrected[which(results_sub_clades$node %in% ascending_nodes)]<=0.05))
    work_A_indep <- length(which(results_sub_clades$pvalue_high_corrected[which(results_sub_clades$node %in% indep_nodes)]<=0.05))
    frac_A_ascending <- work_A_ascending/length(ascending_nodes)
    frac_A_indep <- work_A_indep/length(indep_nodes)
    
    
    
    # test sub-clade B
    results_sub_clades <- phylosignal_sub_network(t(network_combined), tree_B_combined, tree_A_combined, method = "GUniFrac", nperm = 100000, correlation = "Pearson", minimum=10)
    pdf(paste0("sub_clades/test_Mantel_sub_clade_simul_1_guild_B_simul_1_",size,"_mutualism_5_seed_",seed,".pdf"))
    print(plot_phylosignal_sub_network(tree_B_combined, results_sub_clades, legend=TRUE, show.tip.label=FALSE, where="bottomleft"))
    dev.off()
    MRCA_mut <- getMRCA(tree_B_combined, tip = tree_B_combined$tip.label[grep("mut_", tree_B_combined$tip.label)] )
    if (results_sub_clades$pvalue_high_corrected[which(results_sub_clades$node==MRCA_mut)]<=0.05) {work_B <- 1} else {work_B <- 0}
    
    # number of ascending nodes and independent nodes with signal 
    ascending_nodes <- c()
    indep_nodes <- c()
    for (i in (Ntip(tree_B_combined)+1):(Ntip(tree_B_combined)+Nnode(tree_B_combined))){
      extracted_tree <- extract.clade(tree_B_combined,  node=i)
      if (all(tree_B$tip.label %in% extracted_tree$tip.label)){
        if (length(tree_B$tip.label<length(extracted_tree$tip.label))){
          ascending_nodes <- c(ascending_nodes, i)
        }
      }
      if (length(which(tree_B$tip.label %in% extracted_tree$tip.label))==0){
        indep_nodes <- c(indep_nodes, i)
      }
    }
    
    work_B_ascending <- length(which(results_sub_clades$pvalue_high_corrected[which(results_sub_clades$node %in% ascending_nodes)]<=0.05))
    work_B_indep <- length(which(results_sub_clades$pvalue_high_corrected[which(results_sub_clades$node %in% indep_nodes)]<=0.05))
    frac_B_ascending <- work_B_ascending/length(ascending_nodes)
    frac_B_indep <- work_B_indep/length(indep_nodes)
    
    write.table(c(size, seed, work_A, work_B, work_A_ascending, work_A_indep, work_B_ascending, work_B_indep, frac_A_ascending, work_A_indep, frac_B_ascending, work_B_indep),paste0("sub_clades/test_Mantel_sub_clade_simul_1_",size,"_mutualism_5_seed_",seed,".csv"), sep=";", quote=F )
    
    print(paste0("worked: ", seed))
    
    
  }
  
  mclapply(1:20, seed_parrallel, mc.cores = 20)
  
}



# results 

setwd("/Users/bperez/ownCloud/Recherche/These/ISYEB/signal_phylo/simul_1/")

results <- c()

for (size in c("A","B","C","D","E","F")){
  for (seed in 1:20){
    table <- read.table(paste0("sub_clades/test_Mantel_sub_clade_simul_1_",size,"_mutualism_5_seed_",seed,".csv"), sep=";", header=T)
    results <- rbind(results, as.vector(table$x))
  }}

results <- data.frame(results)
table(results[,3]) # 98/120
table(results[,4]) # 99/120

# ascending
table(results[,9]) # 98/120
table(results[,11]) # 99/120

table(results[,5]) # 98/120
table(results[,7]) # 99/120

# indep
table(results[,6]) # 2/120
table(results[,8]) # 0/120

table(results[,10]) # 2/120
table(results[,12]) # 0/120


##### Step 8: Check the correlation between degree and ecological distances ####

# on cluster
rm(list=ls())

library(vegan)
library(phytools)
library(ggplot2)
library(GUniFrac,lib.loc="/users/biodiv/bperez/packages/")
library(parallel)
require(R.utils)


setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1")


source("../script/function_phylo_signal_network.R")
dyn.load("../script/permute.so")

method="Jaccard_weighted"
method="Jaccard_binary"
method="GUniFrac"
method="UniFrac_unweighted"

list_networks <- as.vector(outer(c("simul_1_A","simul_1_B","simul_1_C","simul_1_D","simul_1_E", "simul_1_F"), c(paste0("neutral_seed_",1:100), paste0("mutualism_", as.vector(outer(1:6, 1:20, paste, sep="_seed_"))), paste0("antagonism_",as.vector(outer(1:9, 1:20, paste, sep="_seed_")))), paste, sep="_"))

results_signal <- c()

for (name in list_networks){
  
  network <- read.table(paste0("network_",name,".csv"), header=T,sep=";")
  
  colnames(network) <- seq(1:ncol(network))
  rownames(network) <- seq(1:nrow(network))
  
  tree_A <- read.tree(paste0("network_tree_guild_A_",name,".tre"))
  tree_B <- read.tree(paste0("network_tree_guild_B_",name,".tre"))
  
  network <- network[tree_B$tip.label,tree_A$tip.label]
  
  # Ecological matrix
  # binary Jaccard distances
  if (method=="Jaccard_binary"){
    jaccard_A <- as.matrix(vegan::vegdist(t(network), "jaccard", binary=T))
    jaccard_B <- as.matrix(vegan::vegdist(network, "jaccard", binary=T))
    eco_A <- jaccard_A
    eco_B <- jaccard_B
  }
  
  # quantitative Jaccard distances
  if (method=="Jaccard_weighted"){
    jaccard_A <- as.matrix(vegan::vegdist(t(network), "jaccard", binary=F))
    jaccard_B <- as.matrix(vegan::vegdist(network, "jaccard", binary=F))
    eco_A <- jaccard_A
    eco_B <- jaccard_B
  }
  
  # Unifrac (generalized UniFrac, with alpha=0.5)
  if (method=="GUniFrac"){
    unifrac_B <- GUniFrac::GUniFrac(network, tree = tree_A)
    unifrac_A <- GUniFrac::GUniFrac(t(network), tree = tree_B)
    index=2
    eco_A <- unifrac_A$unifracs[,,index]
    eco_B <- unifrac_B$unifracs[,,index]
  }
  
  # Unifrac (unweighted UniFrac)
  if (method=="UniFrac_unweighted"){
    unifrac_B <- GUniFrac::GUniFrac(network, tree = tree_A)
    unifrac_A <- GUniFrac::GUniFrac(t(network), tree = tree_B)
    index=4
    eco_A <- unifrac_A$unifracs[,,index]
    eco_B <- unifrac_B$unifracs[,,index]
  }
  
  # Degree matrix 
  network[network>0] <- 1 # makes binary
  
  # Mantel test (Pearson)
  mantel_A <- mantel(as.dist(eco_A) ~ dist(colSums(network)),  nperm = 10000, correlation="Pearson")
  mantel_B <- mantel(as.dist(eco_B) ~ dist(rowSums(network)),  nperm = 10000, correlation="Pearson")
  
  
  results_signal <- rbind(results_signal, c(name, mantel_A[1], mantel_A[2], mantel_A[3], mantel_B[1], mantel_B[2], mantel_B[3]))
  
}

results_signal <- data.frame(results_signal)

colnames(results_signal) <- c("name", "mantel_cor_A","pvalue_high_A","pvalue_low_A", "mantel_cor_B", "pvalue_high_B", "pvalue_low_B")

results_signal$size <- "A"
results_signal$size[grep(pattern = "_B_", results_signal$name)] <- "B"
results_signal$size[grep(pattern = "_C_", results_signal$name)] <- "C"
results_signal$size[grep(pattern = "_D_", results_signal$name)] <- "D"
results_signal$size[grep(pattern = "_E_", results_signal$name)] <- "E"
results_signal$size[grep(pattern = "_F_", results_signal$name)] <- "F"
results_signal$param <- "neutral"
results_signal$param[grep(pattern = "_antagonism_1_", results_signal$name)] <- "antagonism_1"
results_signal$param[grep(pattern = "_antagonism_2_", results_signal$name)] <- "antagonism_2"
results_signal$param[grep(pattern = "_antagonism_3_", results_signal$name)] <- "antagonism_3"
results_signal$param[grep(pattern = "_antagonism_4_", results_signal$name)] <- "antagonism_4"
results_signal$param[grep(pattern = "_antagonism_5_", results_signal$name)] <- "antagonism_5"
results_signal$param[grep(pattern = "_antagonism_6_", results_signal$name)] <- "antagonism_6"
results_signal$param[grep(pattern = "_antagonism_7_", results_signal$name)] <- "antagonism_7"
results_signal$param[grep(pattern = "_antagonism_8_", results_signal$name)] <- "antagonism_8"
results_signal$param[grep(pattern = "_antagonism_9_", results_signal$name)] <- "antagonism_9"
results_signal$param[grep(pattern = "_mutualism_1_", results_signal$name)] <- "mutualism_1"
results_signal$param[grep(pattern = "_mutualism_2_", results_signal$name)] <- "mutualism_2"
results_signal$param[grep(pattern = "_mutualism_3_", results_signal$name)] <- "mutualism_3"
results_signal$param[grep(pattern = "_mutualism_4_", results_signal$name)] <- "mutualism_4"
results_signal$param[grep(pattern = "_mutualism_5_", results_signal$name)] <- "mutualism_5"
results_signal$param[grep(pattern = "_mutualism_6_", results_signal$name)] <- "mutualism_6"

write.table(results_signal,paste0("results_correlation_degree_simul_1_",method,".csv"), quote=F, sep=";",row.names=F)



# Plot results 


rm(list=ls()) # on cluster

library(ggplot2)
library(vegan)
library(phytools)
library(ggplot2)
library(GUniFrac,lib.loc="/users/biodiv/bperez/packages/")
library(parallel)
require(R.utils)

setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1")

dyn.load("../script/permute.so")
source("../script/function_phylo_signal_network.R")

dir.create(file.path("./simulations/"), showWarnings = FALSE)



transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))


for (method in c("Jaccard_weighted", "Jaccard_binary", "GUniFrac", "UniFrac_unweighted")){
  
  results_all <- read.table(paste0("results_correlation_degree_simul_1_",method,".csv"), header=T,sep=";")
  
  
  results_all$param <- as.character(results_all$param)
  results_all$size <- as.character(results_all$size)
  results_all$mantel_cor_A <- as.numeric(as.character(results_all$mantel_cor_A))
  results_all$pvalue_high_A <- as.numeric(as.character(results_all$pvalue_high_A))
  results_all$mantel_cor_B <- as.numeric(as.character(results_all$mantel_cor_B))
  results_all$pvalue_high_B <- as.numeric(as.character(results_all$pvalue_high_B))
  
  
  print(method)
  print("A")
  print(length(intersect(which(results_all$pvalue_high_A<0.05), which(results_all$mantel_cor_A>0.10))))  
  
  print("B")
  print(length(intersect(which(results_all$pvalue_high_B<0.05), which(results_all$mantel_cor_B>0.10))))  
  
  # Plot 
  results_signal <- results_all
  
  
  # for tree A
  
  results_signal$size <- as.factor(results_signal$size)
  results_signal$size <-  factor(results_signal$size,levels(results_signal$size)[c(5,1,2,3,4,6)])
  
  ggplot(results_signal,aes(x=param,y=mantel_cor_A,fill=size))+xlab("Parameters")+ylab("Mantel correlation (A)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
    geom_hline(yintercept=0, linetype="dashed", color="#273746")+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    geom_boxplot(alpha=0.9,size=0.75, color="#273746")
  
  ggsave(filename=paste0("boxplot_correlation_degree_Mantel_tests_",method,"_size_correlation_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)
  
  
  # for tree B
  
  ggplot(results_signal,aes(x=param,y=mantel_cor_B,fill=size))+xlab("Parameters")+ylab("Mantel correlation (B)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
    geom_hline(yintercept=0, linetype="dashed", color="#273746")+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    geom_boxplot(alpha=0.9,size=0.75, color="#273746")
  
  ggsave(filename=paste0("boxplot_correlation_degree_Mantel_tests_",method,"_size_correlation_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)
  
  
}



##### Step 8-B: Check the correlation between degree and phylogenetic distances ####

# on cluster
rm(list=ls())

library(vegan)
library(phytools)
library(ggplot2)

library(GUniFrac,lib.loc="/users/biodiv/bperez/packages/")
library(parallel)
require(R.utils)


setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1")


source("../script/function_phylo_signal_network.R")
dyn.load("../script/permute.so")


list_networks <- as.vector(outer(c("simul_1_A","simul_1_B","simul_1_C","simul_1_D","simul_1_E", "simul_1_F"), c(paste0("neutral_seed_",1:100), paste0("mutualism_", as.vector(outer(1:6, 1:20, paste, sep="_seed_"))), paste0("antagonism_",as.vector(outer(1:9, 1:20, paste, sep="_seed_")))), paste, sep="_"))

results_signal <- c()

for (name in list_networks){
  
  network <- read.table(paste0("network_",name,".csv"), header=T,sep=";")
  
  colnames(network) <- seq(1:ncol(network))
  rownames(network) <- seq(1:nrow(network))
  
  tree_A <- read.tree(paste0("network_tree_guild_A_",name,".tre"))
  tree_B <- read.tree(paste0("network_tree_guild_B_",name,".tre"))
  
  network <- network[tree_B$tip.label,tree_A$tip.label]
  
  # Degree matrix 
  network[network>0] <- 1 # makes binary
  
  cophe_A <- cophenetic.phylo(tree_A)
  cophe_B <- cophenetic.phylo(tree_B)
  
  # Mantel test (Pearson)
  mantel_A <- mantel(as.dist(cophe_A) ~ dist(colSums(network)),  nperm = 10000, correlation="Pearson")
  mantel_B <- mantel(as.dist(cophe_B) ~ dist(rowSums(network)),  nperm = 10000, correlation="Pearson")
  
  
  results_signal <- rbind(results_signal, c(name, Ntip(tree_A), Ntip(tree_B),  mantel_A[1], mantel_A[2], mantel_A[3], mantel_B[1], mantel_B[2], mantel_B[3]))
  
}

results_signal <- data.frame(results_signal)

colnames(results_signal) <- c("name", "nb_A","nb_B","mantel_cor_A","pvalue_high_A","pvalue_low_A", "mantel_cor_B", "pvalue_high_B", "pvalue_low_B")

results_signal$size <- "A"
results_signal$size[grep(pattern = "_B_", results_signal$name)] <- "B"
results_signal$size[grep(pattern = "_C_", results_signal$name)] <- "C"
results_signal$size[grep(pattern = "_D_", results_signal$name)] <- "D"
results_signal$size[grep(pattern = "_E_", results_signal$name)] <- "E"
results_signal$size[grep(pattern = "_F_", results_signal$name)] <- "F"
results_signal$param <- "neutral"
results_signal$param[grep(pattern = "_antagonism_1_", results_signal$name)] <- "antagonism_1"
results_signal$param[grep(pattern = "_antagonism_2_", results_signal$name)] <- "antagonism_2"
results_signal$param[grep(pattern = "_antagonism_3_", results_signal$name)] <- "antagonism_3"
results_signal$param[grep(pattern = "_antagonism_4_", results_signal$name)] <- "antagonism_4"
results_signal$param[grep(pattern = "_antagonism_5_", results_signal$name)] <- "antagonism_5"
results_signal$param[grep(pattern = "_antagonism_6_", results_signal$name)] <- "antagonism_6"
results_signal$param[grep(pattern = "_antagonism_7_", results_signal$name)] <- "antagonism_7"
results_signal$param[grep(pattern = "_antagonism_8_", results_signal$name)] <- "antagonism_8"
results_signal$param[grep(pattern = "_antagonism_9_", results_signal$name)] <- "antagonism_9"
results_signal$param[grep(pattern = "_mutualism_1_", results_signal$name)] <- "mutualism_1"
results_signal$param[grep(pattern = "_mutualism_2_", results_signal$name)] <- "mutualism_2"
results_signal$param[grep(pattern = "_mutualism_3_", results_signal$name)] <- "mutualism_3"
results_signal$param[grep(pattern = "_mutualism_4_", results_signal$name)] <- "mutualism_4"
results_signal$param[grep(pattern = "_mutualism_5_", results_signal$name)] <- "mutualism_5"
results_signal$param[grep(pattern = "_mutualism_6_", results_signal$name)] <- "mutualism_6"

write.table(results_signal,paste0("results_correlation_degree_phylo_simul_1.csv"), quote=F, sep=";",row.names=F)


# Plot results 


rm(list=ls()) # on cluster

library(ggplot2)
library(vegan)
library(phytools)
library(ggplot2)

library(GUniFrac,lib.loc="/users/biodiv/bperez/packages/")
library(parallel)
require(R.utils)


setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1")

dyn.load("../script/permute.so")
source("../script/function_phylo_signal_network.R")

transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))



results_all <- read.table(paste0("results_correlation_degree_phylo_simul_1.csv"), header=T,sep=";")


results_all$param <- as.character(results_all$param)
results_all$size <- as.character(results_all$size)
results_all$mantel_cor_A <- as.numeric(as.character(results_all$mantel_cor_A))
results_all$pvalue_high_A <- as.numeric(as.character(results_all$pvalue_high_A))
results_all$mantel_cor_B <- as.numeric(as.character(results_all$mantel_cor_B))
results_all$pvalue_high_B <- as.numeric(as.character(results_all$pvalue_high_B))
results_all$nb_A <- as.numeric(as.character(results_all$nb_A))
results_all$nb_B <- as.numeric(as.character(results_all$nb_B))

print(method)
print("A")
print(length(intersect(which(results_all$pvalue_high_A<0.05), which(results_all$mantel_cor_A>0.10))))  

print("B")
print(length(intersect(which(results_all$pvalue_high_B<0.05), which(results_all$mantel_cor_B>0.10))))  

# Plot 
results_signal <- results_all


# for tree A

results_signal$size <- as.factor(results_signal$size)
results_signal$size <-  factor(results_signal$size,levels(results_signal$size)[c(5,1,2,3,4,6)])

ggplot(results_signal,aes(x=param,y=mantel_cor_A,fill=size))+xlab("Parameters")+ylab("Mantel correlation (A)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
  geom_hline(yintercept=0, linetype="dashed", color="#273746")+
  scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                              bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
  geom_boxplot(alpha=0.9,size=0.75, color="#273746")

ggsave(filename=paste0("boxplot_correlation_degree_phylo_simul_1_correlation_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)


# for tree B

ggplot(results_signal,aes(x=param,y=mantel_cor_B,fill=size))+xlab("Parameters")+ylab("Mantel correlation (B)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
  geom_hline(yintercept=0, linetype="dashed", color="#273746")+
  scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                              bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
  geom_boxplot(alpha=0.9,size=0.75, color="#273746")

ggsave(filename=paste0("boxplot_correlation_degree_phylo_simul_1_correlation_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 7)


library(dplyr)
brks <- c(0, 0.25, 0.5, 0.75, 1)

results_signal$size <- as.factor(results_signal$size)
results_signal$size <-  factor(results_signal$size,levels(results_signal$size)[c(5,1,2,3,4,6)])

results_signal$size_tot <- "<150"
results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>=150)] <- "150-250"
results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>250)] <- ">250"

results_signal$size_tot <- as.factor(results_signal$size_tot)
results_signal$size_tot <-  factor(results_signal$size_tot,levels(results_signal$size_tot)[c(1,3,2)])


# for tree A

results_signal$Inference <- "d) Not significant signal"
results_signal$Inference[which(results_signal$pvalue_high_A<0.05)] <- "e) Significant signal"
results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.05))] <- "f) Significant signal (R>0.05)"
results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.15))] <- "g) Significant signal (R>0.15)"
results_signal$Inference[which(results_signal$pvalue_low_A<0.05)] <- "c) Significant anti-signal"
results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.05)))] <- "b) Significant anti-signal (R<-0.05)"
results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.15)))] <- "a) Significant anti-signal (R<-0.15)"

colors <- c()
if ("a) Significant anti-signal (R<-0.15)" %in% results_signal$Inference) {colors <- c(colors, "#943126")}
if ("b) Significant anti-signal (R<-0.05)" %in% results_signal$Inference) {colors <- c(colors, "#e74c3c")}
if ("c) Significant anti-signal" %in% results_signal$Inference) {colors <- c(colors, "#f1948a")}
if ("d) Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
if ("e) Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
if ("f) Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
if ("g) Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}

results_signal_all <- results_signal %>% 
  group_by(param,Inference,size_tot) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))

results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
results_signal_all$size_tot <-  factor(results_signal_all$size_tot,levels(results_signal_all$size_tot)[c(1,3,2)])

ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
  xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
  scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
  scale_color_manual(values=colors)+
  scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                              bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
  facet_grid(~size_tot)

ggsave(filename=paste0("boxplot_correlation_degree_phylo_simul_1_size_correlation_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 11, height = 7)



##  Resume for clade B

results_signal$Inference <- "d) Not significant signal"
results_signal$Inference[which(results_signal$pvalue_high_B<0.05)] <- "e) Significant signal"
results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.05))] <- "f) Significant signal (R>0.05)"
results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.15))] <- "g) Significant signal (R>0.15)"
results_signal$Inference[which(results_signal$pvalue_low_B<0.05)] <- "c) Significant anti-signal"
results_signal$Inference[intersect(which(results_signal$pvalue_low_B<0.05), which(results_signal$mantel_cor_B<(-0.05)))] <- "b) Significant anti-signal (R<-0.05)"
results_signal$Inference[intersect(which(results_signal$pvalue_low_B<0.05), which(results_signal$mantel_cor_B<(-0.15)))] <- "a) Significant anti-signal (R<-0.15)"

#colors=c("#943126","#e74c3c","#f1948a","#f7dc6f","#7dcea0","#27ae60","#196f3d")
colors <- c()
if ("a) Significant anti-signal (R<-0.15)" %in% results_signal$Inference) {colors <- c(colors, "#943126")}
if ("b) Significant anti-signal (R<-0.05)" %in% results_signal$Inference) {colors <- c(colors, "#e74c3c")}
if ("c) Significant anti-signal" %in% results_signal$Inference) {colors <- c(colors, "#f1948a")}
if ("d) Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
if ("e) Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
if ("f) Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
if ("g) Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}


results_signal_all <- results_signal %>% 
  group_by(param,Inference,size_tot) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))

results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
results_signal_all$size_tot <-  factor(results_signal_all$size_tot,levels(results_signal_all$size_tot)[c(1,3,2)])

ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
  xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
  scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
  scale_color_manual(values=colors)+
  scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                              bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                              bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                              bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                              bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                              bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
  facet_grid(~size_tot)

ggsave(filename=paste0("boxplot_correlation_degree_phylo_simul_1_size_correlation_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 11, height = 7)




### Resume the power when successive Mantel tests ####


rm(list=ls()) # on computer

library(ggplot2)
library(vegan)
library(phytools)
library(ggplot2)

library(GUniFrac,lib.loc="/users/biodiv/bperez/packages/")
library(parallel)
require(R.utils)

setwd("/Users/bperez/ownCloud/Recherche/These/ISYEB/signal_phylo/simul_1/")

dyn.load("../script/permute.so")
source("../script/function_phylo_signal_network.R")

transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))



results_all <- read.table(paste0("results_correlation_degree_phylo_simul_1.csv"), header=T,sep=";")


results_all$param <- as.character(results_all$param)
results_all$size <- as.character(results_all$size)
results_all$mantel_cor_A <- as.numeric(as.character(results_all$mantel_cor_A))
results_all$pvalue_high_A <- as.numeric(as.character(results_all$pvalue_high_A))
results_all$mantel_cor_B <- as.numeric(as.character(results_all$mantel_cor_B))
results_all$pvalue_high_B <- as.numeric(as.character(results_all$pvalue_high_B))
results_all$nb_A <- as.numeric(as.character(results_all$nb_A))
results_all$nb_B <- as.numeric(as.character(results_all$nb_B))

results_all <- results_all[which(results_all$param!="neutral"),]
results_all <- results_all[which(results_all$param!="mutualism_1"),]
results_all <- results_all[which(results_all$param!="mutualism_2"),]
results_all <- results_all[which(results_all$param!="mutualism_4"),]

for (method in c("Jaccard_weighted","Jaccard_binary","GUniFrac","UniFrac_unweighted")){
  
  correlation= "Pearson"
  
  results_signal <- read.table(paste0("results_phylogenetic_signal_simul_1_",method,"_", correlation, ".csv"), header=T, sep=";")
  
  results_signal$mantel_cor_A <- as.numeric(results_signal$mantel_cor_A)
  results_signal$mantel_cor_B <- as.numeric(results_signal$mantel_cor_B)
  results_signal$pvalue_high_A <- as.numeric(results_signal$pvalue_high_A)
  results_signal$pvalue_high_B <- as.numeric(results_signal$pvalue_high_B)
  
  print(method)
  
  # do not keep the neutral (only mutualistic or antagonisitic)
  results_signal <- results_signal[which(results_signal$param!="neutral"),]
  results_signal <- results_signal[which(results_signal$param!="mutualism_1"),]
  results_signal <- results_signal[which(results_signal$param!="mutualism_2"),]
  results_signal <- results_signal[which(results_signal$param!="mutualism_4"),]
  
  print("Power regular Mantel test:")
  print(length(which(results_signal$pvalue_high_A<0.05))/nrow(results_signal))
  print(length(which(results_signal$pvalue_high_B<0.05))/nrow(results_signal))
  
  print("Power separate Mantel tests:")
  print(length(intersect(which(results_signal$pvalue_high_A<0.05), which(results_all$pvalue_high_A>0.05)))/nrow(results_all))
  print(length(intersect(which(results_signal$pvalue_high_B<0.05), which(results_all$pvalue_high_B>0.05)))/nrow(results_all))
  
  
  # Partial Mantel test (Part 8C)
  
  results_signal <- read.table(paste0("results_partial_Mantel_test_degree_simul_1_",method,".csv"), header=T, sep=";")
  
  results_signal$mantel_cor_A <- as.numeric(results_signal$mantel_cor_A)
  results_signal$mantel_cor_B <- as.numeric(results_signal$mantel_cor_B)
  results_signal$pvalue_high_A <- as.numeric(results_signal$pvalue_high_A)
  results_signal$pvalue_high_B <- as.numeric(results_signal$pvalue_high_B)
  
  # do not keep the neutral (only mutualistic or antagonisitic)
  results_signal <- results_signal[which(results_signal$param!="neutral"),]
  results_signal <- results_signal[which(results_signal$param!="mutualism_1"),]
  results_signal <- results_signal[which(results_signal$param!="mutualism_2"),]
  results_signal <- results_signal[which(results_signal$param!="mutualism_4"),]
  
  print("Power partial Mantel test:")
  
  print(length(which(results_signal$pvalue_high_A<0.05))/nrow(results_signal))
  print(length(which(results_signal$pvalue_high_B<0.05))/nrow(results_signal))
  
  print(" ")
  
}


##### Step 8-C: Partial Mantel tests with degree ####

# on cluster
rm(list=ls())

library(vegan)
library(phytools)
library(ggplot2)

library(GUniFrac,lib.loc="/users/biodiv/bperez/packages/")
library(parallel)
require(R.utils)
library(ecodist)


setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1")


source("../script/function_phylo_signal_network.R")
dyn.load("../script/permute.so")

method="Jaccard_weighted"
method="Jaccard_binary"
method="GUniFrac"
method="UniFrac_unweighted"

list_networks <- as.vector(outer(c("simul_1_A","simul_1_B","simul_1_C","simul_1_D","simul_1_E", "simul_1_F"), c(paste0("neutral_seed_",1:100), paste0("mutualism_", as.vector(outer(1:6, 1:20, paste, sep="_seed_"))), paste0("antagonism_",as.vector(outer(1:9, 1:20, paste, sep="_seed_")))), paste, sep="_"))

results_signal <- c()

for (name in list_networks){
  
  network <- read.table(paste0("network_",name,".csv"), header=T,sep=";")
  
  colnames(network) <- seq(1:ncol(network))
  rownames(network) <- seq(1:nrow(network))
  
  tree_A <- read.tree(paste0("network_tree_guild_A_",name,".tre"))
  tree_B <- read.tree(paste0("network_tree_guild_B_",name,".tre"))
  
  network <- network[tree_B$tip.label,tree_A$tip.label]
  
  # Ecological matrix
  # binary Jaccard distances
  if (method=="Jaccard_binary"){
    jaccard_A <- as.matrix(vegan::vegdist(t(network), "jaccard", binary=T))
    jaccard_B <- as.matrix(vegan::vegdist(network, "jaccard", binary=T))
    eco_A <- jaccard_A
    eco_B <- jaccard_B
  }
  
  # quantitative Jaccard distances
  if (method=="Jaccard_weighted"){
    jaccard_A <- as.matrix(vegan::vegdist(t(network), "jaccard", binary=F))
    jaccard_B <- as.matrix(vegan::vegdist(network, "jaccard", binary=F))
    eco_A <- jaccard_A
    eco_B <- jaccard_B
  }
  
  # Unifrac (generalized UniFrac, with alpha=0.5)
  if (method=="GUniFrac"){
    unifrac_B <- GUniFrac::GUniFrac(network, tree = tree_A)
    unifrac_A <- GUniFrac::GUniFrac(t(network), tree = tree_B)
    index=2
    eco_A <- unifrac_A$unifracs[,,index]
    eco_B <- unifrac_B$unifracs[,,index]
  }
  
  # Unifrac (unweighted UniFrac)
  if (method=="UniFrac_unweighted"){
    unifrac_B <- GUniFrac::GUniFrac(network, tree = tree_A)
    unifrac_A <- GUniFrac::GUniFrac(t(network), tree = tree_B)
    index=4
    eco_A <- unifrac_A$unifracs[,,index]
    eco_B <- unifrac_B$unifracs[,,index]
  }
  
  # Degree matrix 
  network[network>0] <- 1 # makes binary
  
  # cophenetic distances
  cophe_A <- cophenetic.phylo(tree_A)
  cophe_B <- cophenetic.phylo(tree_B)
  
  # Partial Mantel test (Pearson)
  mantel_A <- ecodist::mantel(as.dist(cophe_A) ~ as.dist(eco_A) + dist(colSums(network)),  nperm = 10000, mrank=F, nboot = 0) # Pearson
  mantel_B <- ecodist::mantel(as.dist(cophe_B) ~ as.dist(eco_B) + dist(rowSums(network)),  nperm = 10000, mrank=F, nboot = 0) # Pearson
  
  
  results_signal <- rbind(results_signal, c(name, Ntip(tree_A), Ntip(tree_B), mantel_A[1], mantel_A[2], mantel_A[3], mantel_B[1], mantel_B[2], mantel_B[3]))
  
}

results_signal <- data.frame(results_signal)

colnames(results_signal) <- c("name","nb_A","nb_B", "mantel_cor_A","pvalue_high_A","pvalue_low_A", "mantel_cor_B", "pvalue_high_B", "pvalue_low_B")

results_signal$size <- "A"
results_signal$size[grep(pattern = "_B_", results_signal$name)] <- "B"
results_signal$size[grep(pattern = "_C_", results_signal$name)] <- "C"
results_signal$size[grep(pattern = "_D_", results_signal$name)] <- "D"
results_signal$size[grep(pattern = "_E_", results_signal$name)] <- "E"
results_signal$size[grep(pattern = "_F_", results_signal$name)] <- "F"
results_signal$param <- "neutral"
results_signal$param[grep(pattern = "_antagonism_1_", results_signal$name)] <- "antagonism_1"
results_signal$param[grep(pattern = "_antagonism_2_", results_signal$name)] <- "antagonism_2"
results_signal$param[grep(pattern = "_antagonism_3_", results_signal$name)] <- "antagonism_3"
results_signal$param[grep(pattern = "_antagonism_4_", results_signal$name)] <- "antagonism_4"
results_signal$param[grep(pattern = "_antagonism_5_", results_signal$name)] <- "antagonism_5"
results_signal$param[grep(pattern = "_antagonism_6_", results_signal$name)] <- "antagonism_6"
results_signal$param[grep(pattern = "_antagonism_7_", results_signal$name)] <- "antagonism_7"
results_signal$param[grep(pattern = "_antagonism_8_", results_signal$name)] <- "antagonism_8"
results_signal$param[grep(pattern = "_antagonism_9_", results_signal$name)] <- "antagonism_9"
results_signal$param[grep(pattern = "_mutualism_1_", results_signal$name)] <- "mutualism_1"
results_signal$param[grep(pattern = "_mutualism_2_", results_signal$name)] <- "mutualism_2"
results_signal$param[grep(pattern = "_mutualism_3_", results_signal$name)] <- "mutualism_3"
results_signal$param[grep(pattern = "_mutualism_4_", results_signal$name)] <- "mutualism_4"
results_signal$param[grep(pattern = "_mutualism_5_", results_signal$name)] <- "mutualism_5"
results_signal$param[grep(pattern = "_mutualism_6_", results_signal$name)] <- "mutualism_6"

write.table(results_signal,paste0("results_partial_Mantel_test_degree_simul_1_",method,".csv"), quote=F, sep=";",row.names=F)



# Plot results 

# scp bperez@jord.biologie.ens.fr:/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1/results_partial_Mantel_test_degree_simul_1_*.csv  /Users/bperez/ownCloud/Recherche/These/ISYEB/signal_phylo/simul_1/

rm(list=ls()) # on computer

library(ggplot2)
library(vegan)
library(phytools)
library(ggplot2)

library(GUniFrac,lib.loc="/users/biodiv/bperez/packages/")
library(parallel)
require(R.utils)
library(dplyr)
brks <- c(0, 0.25, 0.5, 0.75, 1)

setwd("/Users/bperez/ownCloud/Recherche/These/ISYEB/signal_phylo/simul_1/")
#setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1")

dyn.load("../script/permute.so")
source("../script/function_phylo_signal_network.R")

dir.create(file.path("./simulations/"), showWarnings = FALSE)

transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))
correlation="Pearson"

for (method in c(  "UniFrac_unweighted",  "GUniFrac",  "Jaccard_binary", "Jaccard_weighted")){
  
  results_all <- read.table(paste0("results_partial_Mantel_test_degree_simul_1_",method,".csv"), header=T,sep=";")
  
  
  results_all$param <- as.character(results_all$param)
  results_all$size <- as.character(results_all$size)
  results_all$mantel_cor_A <- as.numeric(as.character(results_all$mantel_cor_A))
  results_all$pvalue_high_A <- as.numeric(as.character(results_all$pvalue_high_A))
  results_all$mantel_cor_B <- as.numeric(as.character(results_all$mantel_cor_B))
  results_all$pvalue_high_B <- as.numeric(as.character(results_all$pvalue_high_B))
  results_all$nb_A <- as.numeric(as.character(results_all$nb_A))
  results_all$nb_B <- as.numeric(as.character(results_all$nb_B))
  
  
  # Plot 
  results_signal <- results_all
  
  results_signal$size <- as.factor(results_signal$size)
  results_signal$size <-  factor(results_signal$size,levels(results_signal$size)[c(5,1,2,3,4,6)])
  
  results_signal$size_tot <- "<150"
  results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>=150)] <- "150-250"
  results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>250)] <- ">250"
  
  #results_signal$size_tot <- as.factor(results_signal$size_tot)
  #results_signal$size_tot <-  factor(results_signal$size_tot,levels(results_signal$size_tot)[c(1,3,2)])
  
  
  # for tree A
  
  results_signal$Inference <- "d) Not significant signal"
  results_signal$Inference[which(results_signal$pvalue_high_A<0.05)] <- "e) Significant signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.05))] <- "f) Significant signal (R>0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.15))] <- "g) Significant signal (R>0.15)"
  results_signal$Inference[which(results_signal$pvalue_low_A<0.05)] <- "c) Significant anti-signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.05)))] <- "b) Significant anti-signal (R<-0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.15)))] <- "a) Significant anti-signal (R<-0.15)"
  
  colors <- c()
  if ("a) Significant anti-signal (R<-0.15)" %in% results_signal$Inference) {colors <- c(colors, "#943126")}
  if ("b) Significant anti-signal (R<-0.05)" %in% results_signal$Inference) {colors <- c(colors, "#e74c3c")}
  if ("c) Significant anti-signal" %in% results_signal$Inference) {colors <- c(colors, "#f1948a")}
  if ("d) Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
  if ("e) Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
  if ("f) Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
  if ("g) Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}
  
  results_signal_all <- results_signal %>% 
    group_by(param,Inference,size_tot) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
  results_signal_all$size_tot <-  factor(results_signal_all$size_tot,levels(results_signal_all$size_tot)[c(1,3,2)])
  
  ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=colors)+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    facet_grid(~size_tot)
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_partial_Mantel_test_degree_simul_1_",method,"_", correlation,"_inference_size_tot_anticorrelation_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 11, height = 7)
  
  
  
  ##  Resume for clade B
  
  results_signal$Inference <- "d) Not significant signal"
  results_signal$Inference[which(results_signal$pvalue_high_B<0.05)] <- "e) Significant signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.05))] <- "f) Significant signal (R>0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.15))] <- "g) Significant signal (R>0.15)"
  results_signal$Inference[which(results_signal$pvalue_low_B<0.05)] <- "c) Significant anti-signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_low_B<0.05), which(results_signal$mantel_cor_B<(-0.05)))] <- "b) Significant anti-signal (R<-0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_low_B<0.05), which(results_signal$mantel_cor_B<(-0.15)))] <- "a) Significant anti-signal (R<-0.15)"
  
  #colors=c("#943126","#e74c3c","#f1948a","#f7dc6f","#7dcea0","#27ae60","#196f3d")
  colors <- c()
  if ("a) Significant anti-signal (R<-0.15)" %in% results_signal$Inference) {colors <- c(colors, "#943126")}
  if ("b) Significant anti-signal (R<-0.05)" %in% results_signal$Inference) {colors <- c(colors, "#e74c3c")}
  if ("c) Significant anti-signal" %in% results_signal$Inference) {colors <- c(colors, "#f1948a")}
  if ("d) Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
  if ("e) Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
  if ("f) Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
  if ("g) Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}
  
  
  results_signal_all <- results_signal %>% 
    group_by(param,Inference,size_tot) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
  results_signal_all$size_tot <-  factor(results_signal_all$size_tot,levels(results_signal_all$size_tot)[c(1,3,2)])
  
  ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=colors)+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    facet_grid(~size_tot)
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_partial_Mantel_test_degree_simul_1_",method,"_", correlation,"_inference_size_tot_anticorrelation_clade_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 11, height = 7)
  
}

# scp bperez@jord.biologie.ens.fr:/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1/*_partial_Mantel_test_degree_simul_1_*.pdf  /Users/bperez/ownCloud/Recherche/These/ISYEB/signal_phylo/simul_1/


##### Step 8-D: Partial Mantel tests with Phylogenetic Permutations ####



# on cluster
rm(list=ls())

library(vegan)
library(phytools)
library(ggplot2)

library(GUniFrac,lib.loc="/users/biodiv/bperez/packages/")
library(parallel)
require(R.utils)


setwd("/users/biodiv/bperez/data/network_project/bipartite/simulations/simul_1")

source("../script/evo_973_sm_phylomantel.r")
source("../script/function_phylo_signal_network.R")
dyn.load("../script/permute.so")


method="Jaccard_weighted"
method="Jaccard_binary"
method="GUniFrac"
method="UniFrac_unweighted"

list_networks <- as.vector(outer(c("simul_1_A","simul_1_B","simul_1_C","simul_1_D","simul_1_E", "simul_1_F"), c(paste0("neutral_seed_",1:100), paste0("mutualism_", as.vector(outer(1:6, 1:20, paste, sep="_seed_"))), paste0("antagonism_",as.vector(outer(1:9, 1:20, paste, sep="_seed_")))), paste, sep="_"))

results_signal <- c()

#for (name in list_networks){
compute_phylo_signal <- function(name, method, correlation, nperm){
  
  network <- read.table(paste0("network_",name,".csv"), header=T,sep=";")
  
  colnames(network) <- seq(1:ncol(network))
  rownames(network) <- seq(1:nrow(network))
  
  tree_A <- read.tree(paste0("network_tree_guild_A_",name,".tre"))
  tree_B <- read.tree(paste0("network_tree_guild_B_",name,".tre"))
  
  network <- network[tree_B$tip.label,tree_A$tip.label]
  
  # Ecological matrix
  # binary Jaccard distances
  if (method=="Jaccard_binary"){
    jaccard_A <- as.matrix(vegan::vegdist(t(network), "jaccard", binary=T))
    jaccard_B <- as.matrix(vegan::vegdist(network, "jaccard", binary=T))
    eco_A <- jaccard_A
    eco_B <- jaccard_B
  }
  
  # quantitative Jaccard distances
  if (method=="Jaccard_weighted"){
    jaccard_A <- as.matrix(vegan::vegdist(t(network), "jaccard", binary=F))
    jaccard_B <- as.matrix(vegan::vegdist(network, "jaccard", binary=F))
    eco_A <- jaccard_A
    eco_B <- jaccard_B
  }
  
  # Unifrac (generalized UniFrac, with alpha=0.5)
  if (method=="GUniFrac"){
    unifrac_B <- GUniFrac::GUniFrac(network, tree = tree_A)
    unifrac_A <- GUniFrac::GUniFrac(t(network), tree = tree_B)
    index=2
    eco_A <- unifrac_A$unifracs[,,index]
    eco_B <- unifrac_B$unifracs[,,index]
  }
  
  # Unifrac (unweighted UniFrac)
  if (method=="UniFrac_unweighted"){
    unifrac_B <- GUniFrac::GUniFrac(network, tree = tree_A)
    unifrac_A <- GUniFrac::GUniFrac(t(network), tree = tree_B)
    index=4
    eco_A <- unifrac_A$unifracs[,,index]
    eco_B <- unifrac_B$unifracs[,,index]
  }
  
  # Degree matrix 
  network[network>0] <- 1 # makes binary
  
  
  # cophenetic distances
  cophe_A <- cophenetic.phylo(tree_A)
  cophe_B <- cophenetic.phylo(tree_B)
  
  # Phylogenetic Mantel test (input matrices)
  
  PP_mantel_A <- phyloMantel(tree_A,m1=cophe_A,m2=eco_A,m3=as.matrix(dist(colSums(network))),nperm=10000)
  PP_mantel_B <- phyloMantel(tree_B,m1=cophe_B,m2=eco_B,m3=as.matrix(dist(rowSums(network))),nperm=10000)
  
  #results_signal <- rbind(results_signal, c(name, Ntip(tree_A), Ntip(tree_B), PP_mantel_A[1], PP_mantel_A[2], PP_mantel_A[3], PP_mantel_B[1], PP_mantel_B[2], PP_mantel_B[3]))
  return(c(name, Ntip(tree_A), Ntip(tree_B), PP_mantel_A[1], PP_mantel_A[2], PP_mantel_A[3], PP_mantel_B[1], PP_mantel_B[2], PP_mantel_B[3]))
  
}

results_signal <- matrix(unlist(mclapply(list_networks,compute_phylo_signal,mc.cores=20, mc.preschedule=F, method=method, correlation="Pearson", nperm=10000)),nrow=length(list_networks),byrow=T)
results_signal <- data.frame(results_signal)

colnames(results_signal) <- c("name","nb_A","nb_B", "mantel_cor_A","pvalue_high_A","pvalue_low_A", "mantel_cor_B", "pvalue_high_B", "pvalue_low_B")

results_signal$size <- "A"
results_signal$size[grep(pattern = "_B_", results_signal$name)] <- "B"
results_signal$size[grep(pattern = "_C_", results_signal$name)] <- "C"
results_signal$size[grep(pattern = "_D_", results_signal$name)] <- "D"
results_signal$size[grep(pattern = "_E_", results_signal$name)] <- "E"
results_signal$size[grep(pattern = "_F_", results_signal$name)] <- "F"
results_signal$param <- "neutral"
results_signal$param[grep(pattern = "_antagonism_1_", results_signal$name)] <- "antagonism_1"
results_signal$param[grep(pattern = "_antagonism_2_", results_signal$name)] <- "antagonism_2"
results_signal$param[grep(pattern = "_antagonism_3_", results_signal$name)] <- "antagonism_3"
results_signal$param[grep(pattern = "_antagonism_4_", results_signal$name)] <- "antagonism_4"
results_signal$param[grep(pattern = "_antagonism_5_", results_signal$name)] <- "antagonism_5"
results_signal$param[grep(pattern = "_antagonism_6_", results_signal$name)] <- "antagonism_6"
results_signal$param[grep(pattern = "_antagonism_7_", results_signal$name)] <- "antagonism_7"
results_signal$param[grep(pattern = "_antagonism_8_", results_signal$name)] <- "antagonism_8"
results_signal$param[grep(pattern = "_antagonism_9_", results_signal$name)] <- "antagonism_9"
results_signal$param[grep(pattern = "_mutualism_1_", results_signal$name)] <- "mutualism_1"
results_signal$param[grep(pattern = "_mutualism_2_", results_signal$name)] <- "mutualism_2"
results_signal$param[grep(pattern = "_mutualism_3_", results_signal$name)] <- "mutualism_3"
results_signal$param[grep(pattern = "_mutualism_4_", results_signal$name)] <- "mutualism_4"
results_signal$param[grep(pattern = "_mutualism_5_", results_signal$name)] <- "mutualism_5"
results_signal$param[grep(pattern = "_mutualism_6_", results_signal$name)] <- "mutualism_6"

write.table(results_signal,paste0("results_PP_partial_Mantel_test_degree_simul_1_",method,".csv"), quote=F, sep=";",row.names=F)



# Plot results 

rm(list=ls()) # on computer

library(ggplot2)
library(vegan)
library(phytools)
library(ggplot2)

library(GUniFrac,lib.loc="/users/biodiv/bperez/packages/")
library(parallel)
require(R.utils)
library(dplyr)
brks <- c(0, 0.25, 0.5, 0.75, 1)

setwd("/Users/bperez/ownCloud/Recherche/These/ISYEB/signal_phylo/simul_1/")

dyn.load("../script/permute.so")
source("../script/function_phylo_signal_network.R")

dir.create(file.path("./simulations/"), showWarnings = FALSE)



transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))
correlation="Pearson"

for (method in c(  "UniFrac_unweighted",  "GUniFrac",   "Jaccard_weighted", "Jaccard_binary")){
  
  results_all <- read.table(paste0("results_PP_partial_Mantel_test_degree_simul_1_",method,".csv"), header=T,sep=";")
  
  
  results_all$param <- as.character(results_all$param)
  results_all$size <- as.character(results_all$size)
  results_all$mantel_cor_A <- as.numeric(as.character(results_all$mantel_cor_A))
  results_all$pvalue_high_A <- as.numeric(as.character(results_all$pvalue_high_A))
  results_all$mantel_cor_B <- as.numeric(as.character(results_all$mantel_cor_B))
  results_all$pvalue_high_B <- as.numeric(as.character(results_all$pvalue_high_B))
  results_all$nb_A <- as.numeric(as.character(results_all$nb_A))
  results_all$nb_B <- as.numeric(as.character(results_all$nb_B))
  
  
  # Plot 
  results_signal <- results_all
  
  results_signal$size <- as.factor(results_signal$size)
  results_signal$size <-  factor(results_signal$size,levels(results_signal$size)[c(5,1,2,3,4,6)])
  
  results_signal$size_tot <- "<150"
  results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>=150)] <- "150-250"
  results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>250)] <- ">250"
  
  # for tree A
  
  results_signal$Inference <- "d) Not significant signal"
  results_signal$Inference[which(results_signal$pvalue_high_A<0.05)] <- "e) Significant signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.05))] <- "f) Significant signal (R>0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.15))] <- "g) Significant signal (R>0.15)"
  results_signal$Inference[which(results_signal$pvalue_low_A<0.05)] <- "c) Significant anti-signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.05)))] <- "b) Significant anti-signal (R<-0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.15)))] <- "a) Significant anti-signal (R<-0.15)"
  
  colors <- c()
  if ("a) Significant anti-signal (R<-0.15)" %in% results_signal$Inference) {colors <- c(colors, "#943126")}
  if ("b) Significant anti-signal (R<-0.05)" %in% results_signal$Inference) {colors <- c(colors, "#e74c3c")}
  if ("c) Significant anti-signal" %in% results_signal$Inference) {colors <- c(colors, "#f1948a")}
  if ("d) Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
  if ("e) Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
  if ("f) Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
  if ("g) Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}
  
  results_signal_all <- results_signal %>% 
    group_by(param,Inference,size_tot) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
  results_signal_all$size_tot <-  factor(results_signal_all$size_tot,levels(results_signal_all$size_tot)[c(1,3,2)])
  
  ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=colors)+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    facet_grid(~size_tot)
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_PP_partial_Mantel_test_degree_simul_1_",method,"_", correlation,"_inference_size_tot_anticorrelation_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 11, height = 7)
  
  
  
  ##  Resume for clade B
  
  results_signal$Inference <- "d) Not significant signal"
  results_signal$Inference[which(results_signal$pvalue_high_B<0.05)] <- "e) Significant signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.05))] <- "f) Significant signal (R>0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_B<0.05), which(results_signal$mantel_cor_B>0.15))] <- "g) Significant signal (R>0.15)"
  results_signal$Inference[which(results_signal$pvalue_low_B<0.05)] <- "c) Significant anti-signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_low_B<0.05), which(results_signal$mantel_cor_B<(-0.05)))] <- "b) Significant anti-signal (R<-0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_low_B<0.05), which(results_signal$mantel_cor_B<(-0.15)))] <- "a) Significant anti-signal (R<-0.15)"
  
  colors <- c()
  if ("a) Significant anti-signal (R<-0.15)" %in% results_signal$Inference) {colors <- c(colors, "#943126")}
  if ("b) Significant anti-signal (R<-0.05)" %in% results_signal$Inference) {colors <- c(colors, "#e74c3c")}
  if ("c) Significant anti-signal" %in% results_signal$Inference) {colors <- c(colors, "#f1948a")}
  if ("d) Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
  if ("e) Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
  if ("f) Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
  if ("g) Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}
  
  
  results_signal_all <- results_signal %>% 
    group_by(param,Inference,size_tot) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  results_signal_all$size_tot <- as.factor(results_signal_all$size_tot)
  results_signal_all$size_tot <-  factor(results_signal_all$size_tot,levels(results_signal_all$size_tot)[c(1,3,2)])
  
  ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=colors)+
    scale_x_discrete(labels = c(bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.1"),
                                bquote("   antagonism: \u03b1"[A]~"=-1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.1; \u03b1"[B]~"=0.01"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=1"),
                                bquote("   antagonism: \u03b1"[A]~"=-0.01; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=1"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=0.01; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.1"),
                                bquote("mutualism: \u03b1"[A]~"=1; \u03b1"[B]~"=0.01"),
                                bquote("mutualism: \u03b1"[A]~"=0.1; \u03b1"[B]~"=0.01"),
                                bquote("neutral: \u03b1"[A]~"=0; \u03b1"[B]~"=0")))+
    facet_grid(~size_tot)
  
  ggsave(filename=paste0("boxplot_phylogenetic_signal_PP_partial_Mantel_test_degree_simul_1_",method,"_", correlation,"_inference_size_tot_anticorrelation_clade_B.pdf"), device=cairo_pdf, plot = last_plot(), width = 11, height = 7)
  
}


##### Step 9A: Simulate phylogenetic signal in generalism using degree ####

rm(list=ls())

library(vegan)
library(phytools)
library(ggplot2)
library(GUniFrac)
library(parallel)
require(R.utils)

setwd("/Users/bperez/ownCloud/Recherche/These/ISYEB/signal_phylo/simul_3/")

source("../script/function_phylo_signal_network.R")
dyn.load("../script/permute.so")

library(mvMORPH)
library(reshape2)

seed=1

# "param_1" d_A=1 
# "param_2" d_A=0.5
# "param_3" d_A=0.05 
# "param_4" d_A=0


method="Jaccard_binary"

for (method in c("Jaccard_binary", "UniFrac_unweighted")){
  
  results_signal <- c()
  
  for (param in paste0("param_",0:4)){
    print(param)
    
    if (param=="param_0"){d_A=5  }  
    if (param=="param_1"){d_A=1  }
    if (param=="param_2"){d_A=0.5  }
    if (param=="param_3"){d_A=0.05}
    if (param=="param_4"){d_A=0  }
    
    for (seed in 1:100){
      
      print(seed)
      
      # simulate phylogenies
      
      set.seed(seed)
      nb_A <- floor(runif(1, 40, 150))
      nb_B <- floor(runif(1, 40, 150))
      
      tree_A <- phytools::pbtree(n=nb_A)
      tree_B <- phytools::pbtree(n=nb_B)
      
      # scale branch length
      tree_A$edge.length <- tree_A$edge.length/max(node.depth.edgelength(tree_A))
      tree_B$edge.length <- tree_B$edge.length/max(node.depth.edgelength(tree_B))
      
      # simulate OU with mean of 0 (theta=0), variance of 0.1 (noise of the brownian motion) and  a parameter of attraction toward 0, d (central tendency)
      
      if (d_A!="star"){
        if (d_A!="0"){
          data_A <- mvSIM(tree_A, param=list(sigma=0.1, alpha=d_A, ntraits=1, theta=0), model="OU1", nsim=1)
        }else{ # brownian motion
          data_A <- mvSIM(tree_A, param=list(sigma=0.1, ntraits=1, theta=0), model="BM1", nsim=1)
        }
      }else{
        star_tree <- di2multi(pbtree(n=Ntip(tree_A)), tol = 10000000)
        star_tree$edge.length <- rep(1,length(star_tree$edge.length))
        data_A <- mvSIM( star_tree, param=list(sigma=0.1, alpha=1, ntraits=1, theta=0), model="OU1", nsim=1)
      }
      
      data_A <- round( 1+(abs(data_A)-min(abs(data_A)))/(max(abs(data_A))-min(abs(data_A)))*(nb_B-1)) # new
      
      if (min(data_A)<1) { print("problem")}
      
      network <- matrix(0, nrow=Ntip(tree_B), ncol= Ntip(tree_A))
      colnames(network) <- tree_A$tip.label
      rownames(network) <- tree_B$tip.label
      
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
      
      # Ecological matrix
      # binary Jaccard distances
      if (method=="Jaccard_binary"){
        jaccard_A <- as.matrix(vegan::vegdist(t(network), "jaccard", binary=T))
        jaccard_B <- as.matrix(vegan::vegdist(network, "jaccard", binary=T))
        eco_A <- jaccard_A
        eco_B <- jaccard_B
      }
      
      # quantitative Jaccard distances
      if (method=="Jaccard_weighted"){
        jaccard_A <- as.matrix(vegan::vegdist(t(network), "jaccard", binary=F))
        jaccard_B <- as.matrix(vegan::vegdist(network, "jaccard", binary=F))
        eco_A <- jaccard_A
        eco_B <- jaccard_B
      }
      
      # Unifrac (generalized UniFrac, with alpha=0.5)
      if (method=="GUniFrac"){
        unifrac_B <- GUniFrac::GUniFrac(network, tree = tree_A)
        unifrac_A <- GUniFrac::GUniFrac(t(network), tree = tree_B)
        index=2
        eco_A <- unifrac_A$unifracs[,,index]
        eco_B <- unifrac_B$unifracs[,,index]
      }
      
      # Unifrac (unweighted UniFrac)
      if (method=="UniFrac_unweighted"){
        unifrac_B <- GUniFrac::GUniFrac(network, tree = tree_A)
        unifrac_A <- GUniFrac::GUniFrac(t(network), tree = tree_B)
        index=4
        eco_A <- unifrac_A$unifracs[,,index]
        eco_B <- unifrac_B$unifracs[,,index]
      }
      
      # Degree matrix 
      network_binary <- network
      network_binary[network_binary>0] <- 1 # makes binary
      
      # cophenetic distances
      cophe_A <- cophenetic.phylo(tree_A)
      cophe_B <- cophenetic.phylo(tree_B)
      
      # Mantel tests (Pearson)
      mantel_A_partial <- ecodist::mantel(as.dist(cophe_A) ~ as.dist(eco_A) + dist(colSums(network_binary)),  nperm = 10000, mrank=F, nboot = 0) # Pearson
      mantel_A_identity <- ecodist::mantel(as.dist(cophe_A) ~ as.dist(eco_A) ,  nperm = 10000, mrank=F, nboot = 0) # Pearson
      mantel_A_degree <- ecodist::mantel(as.dist(cophe_A) ~ dist(colSums(network_binary)),  nperm = 10000, mrank=F, nboot = 0) # Pearson
      
      
      # Mantel test (Pearson)
      mantel_A_eco <- ecodist::mantel(as.dist(eco_A) ~ dist(colSums(network_binary)),  nperm = 10000, mrank=F, nboot = 0) 
      
      results_signal <- rbind(results_signal, c(param, seed, d_A, Ntip(tree_A), Ntip(tree_B), length(which(network>0))/(ncol(network)*nrow(network)) , round(mantel_A_partial[1:3],4), round(mantel_A_identity[1:3],4), round(mantel_A_degree[1:3],4), round(mantel_A_eco[1:3],4))) #,  mantel_B_partial[1:3], mantel_B_identity[1:3], mantel_B_degree[1:3]))
      
    }}
  
  results_signal <- data.frame(results_signal)
  
  colnames(results_signal) <- c("param","seed","d_A", #"d_B", 
                                "nb_A","nb_B","connectance",
                                "partial_mantel_cor_A","partial_pvalue_high_A","partial_pvalue_low_A",
                                "identity_mantel_cor_A","identity_pvalue_high_A","identity_pvalue_low_A",
                                "degree_mantel_cor_A","degree_pvalue_high_A","degree_pvalue_low_A",
                                "eco_mantel_cor_A","eco_pvalue_high_A","eco_pvalue_low_A") 
  
  write.table(results_signal,paste0("results_generalism_partial_Mantel_test_degree_simul_3_new_simul_",method,".csv"), quote=F, sep=";",row.names=F)
  
  
  hist(as.numeric(as.character(results_signal$partial_pvalue_high_A)))
  hist(as.numeric(as.character(results_signal$identity_pvalue_high_A)))
  hist(as.numeric(as.character(results_signal$degree_pvalue_high_A)))
  
  print(table(as.numeric(as.character(results_signal$partial_pvalue_high_A))<0.05))
  print(table(as.numeric(as.character(results_signal$identity_pvalue_high_A))<0.05))
  print(table(as.numeric(as.character(results_signal$degree_pvalue_high_A))<0.05))
  
}



#### Plot results  #####


rm(list=ls()) # on computer

library(ggplot2)
library(vegan)
library(phytools)
library(ggplot2)
library(GUniFrac)
#library(GUniFrac,lib.loc="/users/biodiv/bperez/packages/")
library(parallel)
require(R.utils)
library(dplyr)
brks <- c(0, 0.25, 0.5, 0.75, 1)

setwd("/Users/bperez/ownCloud/Recherche/These/ISYEB/signal_phylo/simul_3/")

dyn.load("../script/permute.so")
source("../script/function_phylo_signal_network.R")


transparent_theme_no <- theme(panel.grid = element_blank(),
                              axis.line = element_line("black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              plot.background = element_rect(fill = "transparent",colour = NA),
                              legend.key=element_blank(),legend.background=element_blank(),
                              axis.text.x = element_text(angle = 60, hjust = 1,size=8) ,
                              axis.title=element_text(face="bold",family="Helvetica",size=10),axis.text=element_text(face="bold"), legend.text = element_text(face="bold"), legend.title = element_text(face="bold"))
correlation="Pearson"

method="GUniFrac"

for (method in c(  "UniFrac_unweighted", "Jaccard_binary")){ #"GUniFrac",   "Jaccard_weighted",
  
  results_all <- read.table(paste0("results_generalism_partial_Mantel_test_degree_simul_3_new_simul_",method,".csv"), header=T,sep=";")
  
  
  results_all$param <- as.character(results_all$param)
  results_all$seed <- as.character(results_all$seed)
  results_all$nb_A <- as.numeric(as.character(results_all$nb_A))
  results_all$nb_B <- as.numeric(as.character(results_all$nb_B))
  results_all$d_A <- (as.character(results_all$d_A))
  
  results_all$partial_mantel_cor_A <- as.numeric(as.character(results_all$partial_mantel_cor_A))
  results_all$partial_pvalue_high_A <- as.numeric(as.character(results_all$partial_pvalue_high_A))
  results_all$partial_pvalue_low_A <- as.numeric(as.character(results_all$partial_pvalue_low_A))
  results_all$identity_mantel_cor_A <- as.numeric(as.character(results_all$identity_mantel_cor_A))
  results_all$identity_pvalue_high_A <- as.numeric(as.character(results_all$identity_pvalue_high_A))
  results_all$identity_pvalue_low_A <- as.numeric(as.character(results_all$identity_pvalue_low_A))
  results_all$degree_mantel_cor_A <- as.numeric(as.character(results_all$degree_mantel_cor_A))
  results_all$degree_pvalue_high_A <- as.numeric(as.character(results_all$degree_pvalue_high_A))
  results_all$degree_pvalue_low_A <- as.numeric(as.character(results_all$degree_pvalue_low_A))
  
  results_all$eco_mantel_cor_A <- as.numeric(as.character(results_all$eco_mantel_cor_A))
  results_all$eco_pvalue_high_A <- as.numeric(as.character(results_all$eco_pvalue_high_A))
  results_all$eco_pvalue_low_A <- as.numeric(as.character(results_all$eco_pvalue_low_A))
  
  
  # Plot 
  results_signal <- results_all
  
  results_signal$size_tot <- "<150"
  results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>=150)] <- "150-200"
  results_signal$size_tot[which(results_signal$nb_A+results_signal$nb_B>200)] <- ">200"
  
  results_signal$size_tot <- as.factor(results_signal$size_tot)
  results_signal$size_tot <-  factor(results_signal$size_tot,levels(results_signal$size_tot)[c(1,3,2)])
  
  
  test="partial_Mantel_test"
  
  # for tree A
  results_signal$mantel_cor_A <- results_signal$partial_mantel_cor_A
  results_signal$pvalue_high_A <- results_signal$partial_pvalue_high_A
  results_signal$pvalue_low_A <- results_signal$partial_pvalue_low_A
  
  results_signal$Inference <- "d) Not significant signal"
  results_signal$Inference[which(results_signal$pvalue_high_A<0.05)] <- "e) Significant signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.05))] <- "f) Significant signal (R>0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.15))] <- "g) Significant signal (R>0.15)"
  results_signal$Inference[which(results_signal$pvalue_low_A<0.05)] <- "c) Significant anti-signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.05)))] <- "b) Significant anti-signal (R<-0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.15)))] <- "a) Significant anti-signal (R<-0.15)"
  
  colors <- c()
  if ("a) Significant anti-signal (R<-0.15)" %in% results_signal$Inference) {colors <- c(colors, "#943126")}
  if ("b) Significant anti-signal (R<-0.05)" %in% results_signal$Inference) {colors <- c(colors, "#e74c3c")}
  if ("c) Significant anti-signal" %in% results_signal$Inference) {colors <- c(colors, "#f1948a")}
  if ("d) Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
  if ("e) Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
  if ("f) Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
  if ("g) Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}
  
  results_signal_all <- results_signal %>% 
    group_by(param,Inference,size_tot) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=colors)+
    scale_x_discrete(labels = c(bquote(" a"[A]~"=5"),
                                bquote(" a"[A]~"=1"),
                                bquote(" a"[A]~"=0.5"),
                                bquote(" a"[A]~"=0.05"),
                                bquote(" a"[A]~"=0")))+
    facet_grid(~size_tot)
  
  ggsave(filename=paste0("boxplot_new_simul_phylogenetic_signal_generalism_",test,"_simul_3_",method,"_", correlation,"_inference_size_tot_anticorrelation_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 5)
  
  
  test="simple_Mantel_test_identity"
  
  # for tree A
  results_signal$mantel_cor_A <- results_signal$identity_mantel_cor_A
  results_signal$pvalue_high_A <- results_signal$identity_pvalue_high_A
  results_signal$pvalue_low_A <- results_signal$identity_pvalue_low_A
  
  results_signal$Inference <- "d) Not significant signal"
  results_signal$Inference[which(results_signal$pvalue_high_A<0.05)] <- "e) Significant signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.05))] <- "f) Significant signal (R>0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.15))] <- "g) Significant signal (R>0.15)"
  results_signal$Inference[which(results_signal$pvalue_low_A<0.05)] <- "c) Significant anti-signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.05)))] <- "b) Significant anti-signal (R<-0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.15)))] <- "a) Significant anti-signal (R<-0.15)"
  
  colors <- c()
  if ("a) Significant anti-signal (R<-0.15)" %in% results_signal$Inference) {colors <- c(colors, "#943126")}
  if ("b) Significant anti-signal (R<-0.05)" %in% results_signal$Inference) {colors <- c(colors, "#e74c3c")}
  if ("c) Significant anti-signal" %in% results_signal$Inference) {colors <- c(colors, "#f1948a")}
  if ("d) Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
  if ("e) Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
  if ("f) Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
  if ("g) Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}
  
  results_signal_all <- results_signal %>% 
    group_by(param,Inference,size_tot) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=colors)+
    scale_x_discrete(labels = c(bquote(" a"[A]~"=5"),
                                bquote(" a"[A]~"=1"),
                                bquote(" a"[A]~"=0.5"),
                                bquote(" a"[A]~"=0.05"),
                                bquote(" a"[A]~"=0")))+
    facet_grid(~size_tot)
  
  ggsave(filename=paste0("boxplot_new_simul_phylogenetic_signal_generalism_",test,"_simul_3_",method,"_", correlation,"_inference_size_tot_anticorrelation_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 5)
  
  
  
  test="simple_Mantel_test_degree"
  
  # for tree A
  results_signal$mantel_cor_A <- results_signal$degree_mantel_cor_A
  results_signal$pvalue_high_A <- results_signal$degree_pvalue_high_A
  results_signal$pvalue_low_A <- results_signal$degree_pvalue_low_A
  
  results_signal$Inference <- "d) Not significant signal"
  results_signal$Inference[which(results_signal$pvalue_high_A<0.05)] <- "e) Significant signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.05))] <- "f) Significant signal (R>0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.15))] <- "g) Significant signal (R>0.15)"
  results_signal$Inference[which(results_signal$pvalue_low_A<0.05)] <- "c) Significant anti-signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.05)))] <- "b) Significant anti-signal (R<-0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.15)))] <- "a) Significant anti-signal (R<-0.15)"
  
  colors <- c()
  if ("a) Significant anti-signal (R<-0.15)" %in% results_signal$Inference) {colors <- c(colors, "#943126")}
  if ("b) Significant anti-signal (R<-0.05)" %in% results_signal$Inference) {colors <- c(colors, "#e74c3c")}
  if ("c) Significant anti-signal" %in% results_signal$Inference) {colors <- c(colors, "#f1948a")}
  if ("d) Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
  if ("e) Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
  if ("f) Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
  if ("g) Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}
  
  results_signal_all <- results_signal %>% 
    group_by(param,Inference,size_tot) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=colors)+
    scale_x_discrete(labels = c(bquote(" a"[A]~"=5"),
                                bquote(" a"[A]~"=1"),
                                bquote(" a"[A]~"=0.5"),
                                bquote(" a"[A]~"=0.05"),
                                bquote(" a"[A]~"=0")))+
    facet_grid(~size_tot)
  
  ggsave(filename=paste0("boxplot_new_simul_phylogenetic_signal_generalism_",test,"_simul_3_",method,"_", correlation,"_inference_size_tot_anticorrelation_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 5)
  
  
  
  #  plot eco
  
  test="correlation_eco_degree_Mantel_tests"
  
  ggplot(results_signal,aes(x=param,y=eco_mantel_cor_A,fill=size_tot))+xlab("Parameters")+ylab("Mantel correlation (A)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
    geom_hline(yintercept=0, linetype="dashed", color="#273746")+
    scale_x_discrete(labels = c(bquote(" a"[A]~"=5"),
                                bquote(" a"[A]~"=1"),
                                bquote(" a"[A]~"=0.5"),
                                bquote(" a"[A]~"=0.05"),
                                bquote(" a"[A]~"=0")))+
    geom_boxplot(alpha=0.9,size=0.75, color="#273746")
  
  ggsave(filename=paste0("boxplot_new_simul_phylogenetic_signal_generalism_",test,"_simul_3_",method,"_", correlation,"_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 5)
  
  
  # for tree A
  results_signal$mantel_cor_A <- results_signal$eco_mantel_cor_A
  results_signal$pvalue_high_A <- results_signal$eco_pvalue_high_A
  results_signal$pvalue_low_A <- results_signal$eco_pvalue_low_A
  
  results_signal$Inference <- "d) Not significant signal"
  results_signal$Inference[which(results_signal$pvalue_high_A<0.05)] <- "e) Significant signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.05))] <- "f) Significant signal (R>0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_high_A<0.05), which(results_signal$mantel_cor_A>0.15))] <- "g) Significant signal (R>0.15)"
  results_signal$Inference[which(results_signal$pvalue_low_A<0.05)] <- "c) Significant anti-signal"
  results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.05)))] <- "b) Significant anti-signal (R<-0.05)"
  results_signal$Inference[intersect(which(results_signal$pvalue_low_A<0.05), which(results_signal$mantel_cor_A<(-0.15)))] <- "a) Significant anti-signal (R<-0.15)"
  
  colors <- c()
  if ("a) Significant anti-signal (R<-0.15)" %in% results_signal$Inference) {colors <- c(colors, "#943126")}
  if ("b) Significant anti-signal (R<-0.05)" %in% results_signal$Inference) {colors <- c(colors, "#e74c3c")}
  if ("c) Significant anti-signal" %in% results_signal$Inference) {colors <- c(colors, "#f1948a")}
  if ("d) Not significant signal" %in% results_signal$Inference) {colors <- c(colors, "#f7dc6f")}
  if ("e) Significant signal" %in% results_signal$Inference) {colors <- c(colors, "#7dcea0")}
  if ("f) Significant signal (R>0.05)" %in% results_signal$Inference) {colors <- c(colors, "#27ae60")}
  if ("g) Significant signal (R>0.15)" %in% results_signal$Inference) {colors <- c(colors, "#196f3d")}
  
  results_signal_all <- results_signal %>% 
    group_by(param,Inference,size_tot) %>% 
    summarise(count=n()) %>% 
    mutate(perc=count/sum(count))
  
  ggplot(results_signal_all,aes(x=param, y = count))+  geom_bar(position="fill", stat="identity", alpha=0.7, aes(fill=Inference,color=Inference))+
    xlab("Parameters")+ylab("Percentage of simulated networks") +labs(title=" ")+scale_fill_manual(values=colors)+transparent_theme_no + 
    scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
    scale_color_manual(values=colors)+
    scale_x_discrete(labels = c(bquote(" a"[A]~"=5"),
                                bquote(" a"[A]~"=1"),
                                bquote(" a"[A]~"=0.5"),
                                bquote(" a"[A]~"=0.05"),
                                bquote(" a"[A]~"=0")))+
    facet_grid(~size_tot)
  
  ggsave(filename=paste0("boxplot_new_simul_phylogenetic_signal_generalism_",test,"_simul_3_",method,"_", correlation,"_inference_size_tot_anticorrelation_clade_A.pdf"), device=cairo_pdf, plot = last_plot(), width = 8, height = 5)
  
  
  
  
  print(method)
  
  results_all_A <- results_all[results_all$d_A!="star",]
  
  results_all_A_restricted <- results_all_A[which(results_all_A$degree_pvalue_high_A<0.05),]
  
  print("type I error (simple Mantel test):")
  print(length(which(results_all_A$identity_pvalue_high_A<0.05))/nrow(results_all_A))
  print(length(which(results_all_A_restricted$identity_pvalue_high_A<0.05))/nrow(results_all_A_restricted))
  
  print("type I error (partial Mantel test):")
  print(length(which(results_all_A$partial_pvalue_high_A<0.05))/nrow(results_all_A))
  print(length(which(results_all_A_restricted$partial_pvalue_high_A<0.05))/nrow(results_all_A_restricted))
  
  print("type I error (successive Mantel tests):")
  print(length(intersect(which(results_all_A$identity_pvalue_high_A<0.05), which(results_all_A$degree_pvalue_high_A>0.05)))/nrow(results_all_A))
  print(length(intersect(which(results_all_A_restricted$identity_pvalue_high_A<0.05), which(results_all_A_restricted$degree_pvalue_high_A>0.05)))/nrow(results_all_A_restricted))
  
  
}


# Plot sizes and connectance

results_all <- read.table(paste0("results_generalism_partial_Mantel_test_degree_simul_3_new_simul_UniFrac_unweighted.csv"), header=T,sep=";")


results_all$connectance <- as.numeric(as.character(results_all$connectance))
results_all$nb_A <- as.numeric(as.character(results_all$nb_A))
results_all$nb_B <- as.numeric(as.character(results_all$nb_B))



# Plot 

ggplot(results_all,aes(x=param,y=nb_A))+xlab("Parameters")+ylab("Number of species (A)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
  ylim(0, max(results_all$nb_A))+
  scale_x_discrete(labels = c(bquote(" a"[A]~"=5"),
                              bquote(" a"[A]~"=1"),
                              bquote(" a"[A]~"=0.5"),
                              bquote(" a"[A]~"=0.05"),
                              bquote(" a"[A]~"=0")))+
  geom_boxplot(alpha=0.9,size=0.75, color="#273746", fill="#eb984e")

ggsave(filename=paste0("boxplot_new_simul_size_A_simul_3.pdf"), device=cairo_pdf, plot = last_plot(), width = 5, height = 4)

ggplot(results_all,aes(x=param,y=nb_B))+xlab("Parameters")+ylab("Number of species (B)") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
  ylim(0, max(results_all$nb_B))+
  scale_x_discrete(labels = c(bquote(" a"[A]~"=5"),
                              bquote(" a"[A]~"=1"),
                              bquote(" a"[A]~"=0.5"),
                              bquote(" a"[A]~"=0.05"),
                              bquote(" a"[A]~"=0")))+
  geom_boxplot(alpha=0.9,size=0.75, color="#273746", fill="#eb984e")

ggsave(filename=paste0("boxplot_new_simul_size_B_simul_3.pdf"), device=cairo_pdf, plot = last_plot(), width = 5, height = 4)


# Connectance
ggplot(results_all,aes(x=param,y=connectance))+xlab("Parameters")+ylab("Connectance") +labs(title=" ")+scale_fill_manual(values=c("#27ae60","#7dcea0","#f4d03f","#f5b041","#e67e22","#ba4a00"))+transparent_theme_no + 
  #geom_hline(yintercept=0, linetype="dashed", color="#273746")+
  scale_x_discrete(labels = c(bquote(" a"[A]~"=5"),
                              bquote(" a"[A]~"=1"),
                              bquote(" a"[A]~"=0.5"),
                              bquote(" a"[A]~"=0.05"),
                              bquote(" a"[A]~"=0")))+
  geom_boxplot(alpha=0.9,size=0.75, color="#273746", fill="#eb984e")

ggsave(filename=paste0("boxplot_new_simul_connectance_simul_3.pdf"), device=cairo_pdf, plot = last_plot(), width = 5, height = 4)






##### Step 10: Empirical application on orchid network from La Réunion ####


rm(list=ls())


setwd("/Users/bperez/ownCloud/Recherche/These/ISYEB/signal_phylo/orchids_Reunion/")

library(RPANDA)
data(mycorrhizal_network)

tree_Fungi <- mycorrhizal_network[[3]]
tree_Orchids_calibrated <- mycorrhizal_network[[2]]
network <- mycorrhizal_network[[1]]

# sub-clade
tree_Orchids_calibrated <- ladderize(tree_Orchids_calibrated)
results_sub_clades <- phylosignal_sub_network(network, tree_Orchids_calibrated, tree_Fungi, method = "GUniFrac", nperm = 100000, correlation = "Pearson", minimum=10)
pdf("test_phylogenetic_signal_orcids_species_levels_GUniFrac.pdf", width=10, height=14)
plot_phylosignal_sub_network(tree_Orchids_calibrated, results_sub_clades, legend=TRUE, show.tip.label=TRUE, where="bottomleft")
dev.off()

pdf("test_phylogenetic_signal_orcids_species_levels_GUniFrac_plot.pdf", width=9, height=12)
plot_phylosignal_sub_network(tree_Orchids_calibrated, results_sub_clades, legend=TRUE, show.tip.label=F, where="bottomleft")
axisPhylo(las = 1)
dev.off()


pdf("plot_tree.pdf", width=10.53, height=14)
plot(tree_Orchids_calibrated, show.tip.label=F,edge.width = 3)
dev.off()

# Jaccard
results_sub_clades <- phylosignal_sub_network(network, tree_Orchids_calibrated, tree_Fungi, method = "Jaccard_weighted", nperm = 100000, correlation = "Pearson", minimum=10)
pdf("test_phylogenetic_signal_orcids_species_levels_Jaccard.pdf", width=10, height=14)
plot_phylosignal_sub_network(tree_Orchids_calibrated, results_sub_clades, legend=TRUE, show.tip.label=TRUE, where="bottomleft")
dev.off()

pdf("test_phylogenetic_signal_orcids_species_levels_Jaccard_plot.pdf", width=9, height=12)
plot_phylosignal_sub_network(tree_Orchids_calibrated, results_sub_clades, legend=TRUE, show.tip.label=F, where="bottomleft")
dev.off()



trait_info <- interactions[,c(1,2)]
trait_info <- trait_info[!duplicated(trait_info$Orchid.species),]
trait_info_all <- trait_info$Type
names(trait_info_all) <- trait_info$Orchid.species
plot(tree_Orchids_calibrated)

pdf("tree_state_epiphytic.pdf", width=10, height=14)
tree_Orchids_calibrated_plot <- ladderize(tree_Orchids_calibrated)
trait_info_all <- trait_info_all[tree_Orchids_calibrated_plot$tip.label]
tree_Orchids_calibrated_plot$tip.label <- paste0("   ", tree_Orchids_calibrated_plot$tip.label)
plot(tree_Orchids_calibrated_plot)
tiplabels(pie=to.matrix(trait_info_all,sort(unique(trait_info_all))),piecol=c("#1e8449","#935116"),cex=0.3)
dev.off()



# Correct for uncertainty in Angraecinae
tree_Orchids_polytomies <- tree_Orchids_calibrated

# Add polytomies in Angraecum

graft_tree <- di2multi(pbtree(n=length(grep("Angraecum", tree_Orchids_polytomies$tip.label))),tol=500)
graft_tree$tip.label <- tree_Orchids_polytomies$tip.label[grep("Angraecum", tree_Orchids_polytomies$tip.label)]
graft_tree$edge.length <- graft_tree$edge.length/graft_tree$edge.length * 0.15
tree_Orchids_polytomies$tip.label[grep("Angraecum", tree_Orchids_polytomies$tip.label)] <- paste0(tree_Orchids_polytomies$tip.label[grep("Angraecum", tree_Orchids_polytomies$tip.label)],"_duplicated")
tree_Orchids_polytomies <- bind.tree(tree_Orchids_polytomies, graft_tree, where=which(tree_Orchids_polytomies$tip.label=="Angraecum_eburneum_duplicated"),position = 0.10)
tree_Orchids_polytomies <- drop.tip(tree_Orchids_polytomies,tip=grep("_duplicated", tree_Orchids_polytomies$tip.label))
tree_Orchids_polytomies$edge.length[which(tree_Orchids_polytomies$edge[,2]==getMRCA(tree_Orchids_polytomies, tip=tree_Orchids_polytomies$tip.label[grep("Angraecum", tree_Orchids_polytomies$tip.label)]))] <- tree_Orchids_polytomies$edge.length[which(tree_Orchids_polytomies$edge[,2]==getMRCA(tree_Orchids_polytomies, tip=tree_Orchids_polytomies$tip.label[grep("Angraecum", tree_Orchids_polytomies$tip.label)]))] -0.05

write.tree(tree_Orchids_polytomies, "tree_orchids_species_level_polytomies_Angraecum.tre")


# sub-clade with Angraecum as polytomy
results_sub_clades_poly <- phylosignal_sub_network(network, tree_Orchids_polytomies, tree_Fungi, method = "GUniFrac", nperm = 100000, correlation = "Pearson", minimum=10)
pdf("test_phylogenetic_signal_orcids_species_levels_polytomies_Angraecum.pdf", width=10, height=14)
plot_phylosignal_sub_network(tree_Orchids_polytomies, results_sub_clades_poly, legend=TRUE, show.tip.label=TRUE, where="bottomleft")
dev.off()



# sub-clade with Angraecum when subsampling only 15 Angraecineae species 
set.seed(15)
discarded_species <- sample(size=34-15, tree_Angraecinae$tip.label)
network_subsampled <- network[,which(!colnames(network) %in% discarded_species)]
tree_Orchids_polytomies_subsampled <- drop.tip(tree_Orchids_polytomies, tip=discarded_species)
results_sub_clades_subsampled <- phylosignal_sub_network(network_subsampled, tree_Orchids_polytomies_subsampled, tree_Fungi, method = "GUniFrac", nperm = 100000, correlation = "Pearson", minimum=10)
pdf("test_phylogenetic_signal_orcids_species_levels_subsampled_Angraecum.pdf", width=10, height=14)
plot_phylosignal_sub_network(tree_Orchids_polytomies_subsampled, results_sub_clades_subsampled, legend=TRUE, show.tip.label=TRUE, where="bottomleft")
dev.off()

set.seed(5)
discarded_species <- sample(size=34-10, tree_Angraecinae$tip.label)
network_subsampled <- network[,which(!colnames(network) %in% discarded_species)]
tree_Orchids_polytomies_subsampled <- drop.tip(tree_Orchids_polytomies, tip=discarded_species)
results_sub_clades_subsampled <- phylosignal_sub_network(network_subsampled, tree_Orchids_polytomies_subsampled, tree_Fungi, method = "GUniFrac", nperm = 100000, correlation = "Pearson", minimum=10)
pdf("test_phylogenetic_signal_orcids_species_levels_subsampled_Angraecum_10.pdf", width=10, height=14)
plot_phylosignal_sub_network(tree_Orchids_polytomies_subsampled, results_sub_clades_subsampled, legend=TRUE, show.tip.label=TRUE, where="bottomleft")
dev.off()


# Apply PBLM 
phylosignal_network(network, tree_Orchids_calibrated, tree_Fungi, method = "PBLM")
# significant for fungi 

