######## Dependances ###############################################



######## Simulation ###############################################

# nx,ny matrix size
# NG number of time steps
# dSpace size of dispersal kernel
# D dimention of trait space 
# iniP, iniH initial value of traits
# muP, muH mutation rates
# effect standard deviation of mutation kernel
# alphaP, alphaH effect of trait difference on fitness
# nP, nH number of dead individuals per time step
# oneByOne are all traits mutated at the same time? (or several/one multidimentional trait)
# death is the fitness acting on birth or on death?
# timeStep used for time execution pb. A list of matrices is created instead of on giant matrix
# verbose the time is printed every verbose time steps
# f fitness function type
# R used if f==4
# thin the state is recorded every thin /!\ must be tested with new version
# P, H used to continue one precedent run: traits of the individuals at the end of the precedent run


model.spatial=function(nx,ny=nx,NG,D=1,muP,muH,alphaP=0,alphaH=0,iniP=0,iniH=0,nP=1,nH=1,rP=1,rH=1,f=3,effect=1,
                       verbose=100,timeStep=floor(sqrt(NG/(nH+nP))),R=1,oneByOne=F,thin=1,P=NULL,H=NULL,fitness_death=FALSE){
  
  N=nx*ny     # number of individuals
  
  parentsP=1:N
  parentsH=1:N
  pmut_act=rep(0,N)
  hmut_act=rep(0,N)
  changeP=c()
  changeH=c()
  
  
  # choice of the fitness function for birth
  if (fitness_death==FALSE){
  f=function(X,x,r,alpha,pos){
    .Call("fitnessFunction", X=as.numeric(t(X)), x=as.numeric(x),r=as.numeric(r),alpha=as.numeric(alpha), Ncol=as.integer(N), D=as.integer(D))
  }
  }else{
    f=function(X,x,r,alpha,pos){rep(1,N)}
  }
  
  # choice of the fitness function for death
  if (fitness_death==TRUE){
  f_death=function(X,x,r,alpha){
    .Call("fitnessFunctionDeath", X=as.numeric(t(X)), x=as.numeric(t(x)),r=as.numeric(r),alpha=as.numeric(alpha), Ncol=as.integer(N), D=as.integer(D))
  }
  }else{
    f_death=function(X,x,r,alpha){rep(1,N)}
  }
  
  # trait initial values
  if(is.null(P)) P=matrix(iniP,nrow=D,ncol=N)
  if(is.null(H)) H=matrix(iniH,nrow=D,ncol=N)
  
  # creation of lists to help with the execution time
  Phist=list(a=lapply(1:D,function(i){Matrix(0,nrow=1,ncol=N,sparse=T)}))
  Hhist=list(a=lapply(1:D,function(i){Matrix(0,nrow=1,ncol=N,sparse=T)}))
  for(j in 1:D){
    Phist[[1]][[j]][1,]=P[j,]
    Hhist[[1]][[j]][1,]=H[j,]
  }
  Pmut=list(a=Matrix(0,nrow=1,ncol=N,sparse=T))
  Hmut=Pmut
  Pgen=Pmut
  Hgen=Pmut
  
  
  listInd=1
  time=0
  t=0     # number of time steps executed in the loop
  
  # beginning of the main loop
  for(u in 1:max(1,ceiling(NG/(timeStep*thin)))){
    
    NTime=ceiling((min(u*timeStep*thin,NG)-time)/thin) # number of time steps in one loop run
    
    t=0
    
    # creation of the matrices in which we record what happens 
    phist=lapply(1:D,function(e){Matrix(0,nrow=NTime,ncol=N,sparse=T)})  # trait of P when there is a birth/death
    hhist=phist         # trait of H when there is a birth/death
    pmut=phist[[1]]     # number of mutations for P
    hmut=phist[[1]]     # number of mutations for H
    pgen=phist[[1]]     # parents for P
    hgen=phist[[1]]     # parents for H
    
    
    # begining of the internal loop _ execution between two timeStep
    while(time<min(u*timeStep*thin,NG)){
      
      # initialisation
      if((t) %% thin ==0){
        parentsP=1:N
        parentsH=1:N
        pmut_act=rep(0,N)
        hmut_act=rep(0,N)
        changeP=c()
        changeH=c()
      }
      
      t=t+1
      time=time+1

      for(k in 1:nH){
        
        prob=f_death(H,P,rH,alphaH)
        if(all(prob==0)) { prob[all(prob==0)] <- 1 }
        deadH=sample(1:N,1,prob = prob)
        
        p=P[,deadH]
        prob=f(H,p,rH,alphaH,deadH)
        if(all(prob==0)) { prob[all(prob==0)] <- 1 }
        newH=sample(1:N,1,prob = prob)
        
        mutH=(runif(1) < muH)     # mutation ?
        H[,deadH]=H[,newH]        # new H trait values
        changeH=unique(c(changeH,deadH,newH))    # which H changed since the last recording
        parentsH[deadH]=parentsH[newH]           # parents of new Hs
        hmut_act[deadH]=hmut_act[parentsH[newH]]+mutH  # number of mutations for H since last recording
        
        # modification of the trait by mutation for H
        if(sum(mutH)>0){

          H[,deadH]=H[,deadH]+rnorm(D,mean = 0,sd=effect)
        }
      }
      
      for(k in 1:nP){
      
        prob=f_death(P,H,rP,alphaP)
        if(all(prob==0)) { prob[all(prob==0)] <- 1 }
        deadP=sample(1:N,1,prob = prob)
        
        h=H[,deadP]
        prob=f(P,h,rP,alphaP,deadP)
        if(all(prob==0)) {prob[all(prob==0)] <- 1}
        newP=sample(1:N,1,prob = prob)
        
        mutP=(runif(1) < muP)     # mutation ?
        P[,deadP]=P[,newP]        # new H trait values
        changeP=unique(c(changeP,deadP,newP))    # which H changed since the last recording
        parentsP[deadP]=parentsP[newP]           # parents of new Hs
        pmut_act[deadP]=pmut_act[parentsP[newP]]+mutP  # number of mutations for H since last recording
        
        # modification of the trait by mutation for H
        if(sum(mutP)>0){
          P[,deadP]=P[,deadP]+rnorm(D,mean = 0,sd=effect)
        }
      }
      
      # recording
      if((t) %% thin ==0){
        for(i in 1:D){
          phist[[i]][t/thin,changeP]=P[i,changeP]
          hhist[[i]][t/thin,changeH]=H[i,changeH]
        }
        
        pgen[t/thin,changeP]=parentsP[changeP]
        hgen[t/thin,changeH]=parentsH[changeH]
        
        pmut[t/thin,]=pmut_act
        hmut[t/thin,]=hmut_act
      }
      
    }
    
    Phist[[listInd]]=phist
    Hhist[[listInd]]=hhist
    Pmut[[listInd]]=pmut
    Hmut[[listInd]]=hmut
    Pgen[[listInd]]=pgen
    Hgen[[listInd]]=hgen
    
    listInd=listInd+1
  }
  
  return(list(Pgenealogy=Pgen,Hgenealogy=Hgen,xP=Phist,xH=Hhist,P=P,H=H,Pmut=Pmut,Hmut=Hmut,iniP=iniP,iniH=iniH,thin.factor=thin))
}



######## Construct the genealogy ###################################

# out an object created by the simulation function
# treeP, treeH former genealogy the tree have to be graphted on
# timeStep used for time execution pb


make.gen=function(out,time.step=nrow(out$Pmut), treeP=NULL, treeH=NULL, verbose=T){
  
  N=ncol(out$Pmut[[1]])                                                  # number of individuals
  NG=sum(sapply(1:length(out$Pmut),function(i){nrow(out$Pmut[[i]])}))    # total time steps number
  D=length(out$xP[[1]])                                                  # dimention of trait space

  # auxiliary function creating a genalogy for one type
  aux=function(P.trait,Gen,X,Mut,ini,thin,tree){
    
    edgeP=matrix(ncol=2)                       # the edges of the genealogy tree
    xP=lapply(1:D,function(i){c()})            # trait values at the beginning of the edges
    currentP=1:N                               # individuals at current time that have descendants in the present
    node=N+1                                   # current node name (1:N are the tips)
    nodes=1:N                                  # current nodes
    nodes_age=rep(NG+1,N)                      # current nodes age
    edge.length=c()                            # length of the edges of the genealogy
    nMut=rep(0,N)                              # number of mutation between a node and its parent node
    t=NG                                       # actual time
    
    listInd=length(Mut)+1                      # indice of the current matrix in the list
    time=NG                                    # time at the first row of the current matrix
    newInd=NULL                                # who was born at time t
    
    # begining of main loop
    while (t>0 & length(currentP)>1){
    # ie we are not yet at the root and there is more than one individual with descendant in the present
      
      listInd=listInd-1
      time.step=nrow(Mut[[listInd]])
      time=time-time.step
      if(time<0){
        time.step=time.step+time
        time=0
      }
      
      x=lapply(1:D,function(i){X[[listInd]][[i]]})  # trait values
      gen=Gen[[listInd]]                            # parents
      mut=Mut[[listInd]]                            # mutation number
      
      # internal loop
      while(t>time){
      # ie we are not at the first row of the matrix yet

        if(is.null(newInd)){
          newInd=which(gen[t-time,currentP]>0)
        }
        
        # if there was a birth at time t
        if (sum(newInd)>0){
          newCurrentP=currentP
          newNodes=nodes
          fathers=gen[t-time,currentP[newInd]]    # parents of the new individuals
          
          for(i in 1:length(newInd)){
            father=fathers[i]
            newCurrentP[newInd[i]]=father
            nMut[nodes[newInd[i]]]=nMut[nodes[newInd[i]]]+(mut[t-time,currentP])[newInd[i]]
            
            if(sum(fathers[1:i]==father)==1){
            # ie the individual do not have the same parent than another born at this time already treated
              
              if(!(father %in% currentP[newInd]) & father %in% currentP){
              # ie the parent is in currentP (else there is a pb !) and did not die at time t 
                
                # we add a new node 
                newNodes[newInd[i]]=node
                newNodes[currentP==father]=node
                nodes_age=c(nodes_age,t)
                nMut=c(nMut,0)
                
                # and the corresponding offspring edges
                # edge 1
                edgeP=rbind(edgeP,c(node,nodes[newInd[i]]))
                for(j in 1:D){
                  xP[[j]]=c(xP[[j]],(x[[j]][t-time,currentP])[newInd[i]])
                }
                edge.length=c(edge.length,nodes_age[nodes[newInd[i]]]-t)

                # edge 2
                edgeP=rbind(edgeP,c(node,nodes[currentP==father]))
                for(j in 1:D){
                  xP[[j]]=c(xP[[j]],x[[j]][t-time,father])
                }
                edge.length=c(edge.length,nodes_age[nodes[currentP==father]]-t)
                
                node=node+1
 
              }else{
              # ie the parent died at time t
                
                if(sum(fathers==father)>1){
                # ie there is another current individual born from the same partent at time t
                  
                  # we add a new node
                  newNodes[newInd[i]]=node
                  nodes_age=c(nodes_age,t)
                  nMut=c(nMut,0)
                  
                  # and the first of its offspring edges
                  edgeP=rbind(edgeP,c(node,nodes[newInd[i]]))
                  for(j in 1:D){
                    xP[[j]]=c(xP[[j]],(x[[j]][t-time,currentP])[newInd[i]])
                  }
                  edge.length=c(edge.length,nodes_age[nodes[newInd[i]]]-t)

                  node=node+1
                }
                # else we don't do anything, the current individual identity is only changed in current P (already done)
              }
            }else{
            # ie there is an individual born at time t (and already seen) with the same parent
              
              if(!(father %in% currentP[newInd]) & father %in% currentP){
              # ie the parent did not die at time t and is in currentP (always, else there is a pb),  we select the parent node
                Node=newNodes[currentP==father]
              }else{
              # ie the parent died at time t, we select the parent node
                Node=(newNodes[newCurrentP==father])[1]
              }
              
              #and add the corresponding edge
              edgeP=rbind(edgeP,c(Node,nodes[newInd[i]]))
              for(j in 1:D){
                xP[[j]]=c(xP[[j]],(x[[j]][t-time,currentP])[newInd[i]])
              }
              newNodes[newInd[i]]=Node
              edge.length=c(edge.length,nodes_age[nodes[newInd[i]]]-t)
            }
          }

          currentP=unique(newCurrentP)
          nodes=unique(newNodes)
        }
        
        # select next time
        if(t==(time+1)){
        # ie the matrix has no line left
          t=t-1
          newInd=NULL
        }else{
          newInd=which(gen[t-time-1,currentP]>0)
          if(sum(newInd)>0 ){
          # ie there is a change in the next time step
            t=t-1
          }else{
            # we look for the next change
            newInd=NULL
            if(t>1 & length(currentP)>1){
              mat=which(gen[1:(t-time-1),currentP]>0,arr.ind = T)
              if(length(mat)==0) {
                t=time
              }else{
                t=max(mat[,1])+time
                ind=mat[mat[,1]==(t-time),2]
              }
            }else{
              t=time
            }}}
      } # end internal loop
      } # end main loop
    
    edgeP=edgeP[-1,]   # edgeP was initialized with a row of NA, we suppress it
    
    if(length(currentP)>1){
    # ie all the individual did not coalesce
      
      if(is.null(tree)){
        for(i in 1:length(currentP)){
          edgeP=rbind(edgeP,c(node,nodes[i]))
          for(j in 1:D){xP[[j]]=c(xP[[j]],ini)}
          edge.length=c(edge.length,nodes_age[nodes[i]])
        }
        
      }else{
        # we graft the new tree on the former tree
        M=max(max(edgeP),N)
        newmut=rep(0,max(tree$edge))
        for(j in 1:nrow(tree$edge)){ newmut[tree$edge[j,2]]=tree$nMut[j]}
        nMut=c(nMut,newmut)
        tree$edge=tree$edge+M
        edgeP=rbind(edgeP,tree$edge)
        edge.length=c(edge.length,tree$edge.length/thin)
        for(j in 1:D) xP[[j]]=c(xP[[j]],tree$x[[j]])
        for (i in 1:length(currentP)){
          ind=which(tree$tip.label==currentP[i])+M
          ind.edge=which(edgeP[,2]==ind)
          nMut[nodes[i]]=nMut[nodes[i]]+nMut[edgeP[ind.edge,2]]
          edgeP[ind.edge,2]=nodes[i]
          edge.length[ind.edge]=edge.length[ind.edge]+nodes_age[nodes[i]]-1
        }
      }
    }
    
    # create a phylo object
    if(is.null(tree) | length(currentP)==1){
      treeP=list(Nnode=max(edgeP)-N,edge=edgeP,tip.label=1:N,edge.length=(edge.length)*thin)
    }else{
      treeP=list(Nnode=length(unique(as.vector(edgeP)))-(2*N-length(currentP)),edge=edgeP,tip.label=c(1:(max(edgeP))),edge.length=(edge.length)*thin)
    }
    class(treeP)="phylo"
    X=rep(0,nrow(edgeP))
    nMut=nMut[treeP$edge[,2]]
    traits=xP
    traits[[D+1]]=nMut
 
    # make the tree well conformed and supress extinct lineages
    treeP=rigth.order(treeP,traits)
    if(!is.null(tree) & (2*N-length(currentP))>N){
      extinct=treeP$tree$tip.label[treeP$tree$tip.label>N]
      treeP=prune.with.traits(treeP$tree,extinct,treeP$trait,extensiveTrait = D+1)}

    # add the traits in the phylo object
    xP=treeP$trait[1:D]
    nMut=treeP$trait[[D+1]]
    x.tip=matrix(0,ncol=N,nrow=D)
    for(j in 1:D){
      x.tip[j,]=P.trait[j,treeP$tree$tip.label]
    }

    treeP=treeP$tree
    treeP$x=xP
    treeP$ini=max(node.depth.edgelength(treeP))
    treeP$nMut=nMut
    treeP$x.tip=x.tip
    
    return(treeP)
  }
  
  P=aux(out$P,out$Pgenealogy,out$xP,out$Pmut,out$iniP,out$thin.factor,treeP)

  H=aux(out$H,out$Hgenealogy,out$xH,out$Hmut,out$iniH,out$thin.factor,treeH)
  return(list(P=P,H=H))
}



######## Create a species tree from the genealogy #################

# genealogy an object created with make.gen
# threshold number of mutation to be in different species
# distanceP, distanceH distance (ie nb of mutations) matrix between the individual
# if NULL it is computed within the function


define.species=function(genealogy,threshold=1,distanceH=NULL,distanceP=NULL, verbose=T, monophyly=TRUE, seed=NULL){
  
  if (monophyly==TRUE){
  
  D=length(genealogy$P$x)       # dimension of trait space

  # auxiliary function to define the species for one type
  aux = function(gen,distance){
      N=length(gen$tip.label)    # number of individuals
    
      # first we define the genetic types of the individuals
    
      if(threshold==1){
          type=rep(1,N)             
          new.type=2
          for(i in 1:nrow(gen$edge)){
            if(gen$nMut[i]>0){
              if(gen$edge[i,2]<=N){
                type[gen$tip.label[gen$edge[i,2]]]=new.type
                new.type=new.type+1
              }else{
                type[as.integer(extract.clade(gen,gen$edge[i,2])$tip.label)]=new.type}
              new.type=new.type+1
          }
        }
     }else{ #threshold >1
        if(is.null(distance)){
          # we compute the distance matrix if it is not given as an argument
          distance=matrix(0,nrow=N,ncol=N)
          for(i in 1:nrow(gen$edge)){
            if(gen$nMut[i]>0){
              #add the number of mutation between the individuals separated by this edge
              if(gen$edge[i,2]<=N){
                tip=as.integer(gen$tip.label[gen$edge[i,2]])
                distance[tip,(1:N)[-tip]]=distance[tip,(1:N)[-tip]]+gen$nMut[i]
                distance[(1:N)[-tip],tip]=distance[(1:N)[-tip],tip]+gen$nMut[i]
              }else{
                tip=as.integer(extract.clade(gen,gen$edge[i,2])$tip.label)
                distance[tip,(1:N)[-tip]]=distance[tip,(1:N)[-tip]]+gen$nMut[i]
                distance[(1:N)[-tip],tip]=distance[(1:N)[-tip],tip]+gen$nMut[i]
              }
            }
          }
          distance=distance[as.integer(gen$tip.label),as.integer(gen$tip.label)]
        }
        
        # we use the distance matrix to define the genetic type of individuals
        phylo=gen
        phylo$edge.length=phylo$nMut
        type=1:N
        for(i in 1:N){
          ind=((1:i)[distance[i,1:i]<threshold])[1]
          type[phylo$tip.label[i]]=type[phylo$tip.label[ind]]

        }
    }
    
    # species are then defined as monophyletic types
    species=rep(1,N)        # the species of each tip
    newSpecies=2
    nextNode=N+1
    
    # test of the tips descendig from a particular node _ main loop
    while(length(nextNode)>0){
      i=nextNode[1]                                       # candidate node
      nextNode=nextNode[-1]                               # update candidate nodes vector
      
      a=extract.clade(gen,i)                              # subtree with root i
      ty=as.integer(gen$tip.label)
      ty=ty[which(!(ty %in% as.integer(a$tip.label)))]
      ty=type[ty]                                         # type of tips not in a
      
      if(length(intersect(type[as.integer(a$tip.label)],ty))==0){
      # seeing how nextNode is built, this should always be the case... Test needed
        n=length(a$tip.label)                             # number of tips in the subtree
        offspring=gen$edge[gen$edge[,1]==i,2]             # offspring nodes of i in the main tree
        offspringSubTree=a$edge[a$edge[,1]==(n+1),2]      # offspring node of i in the subtree
        
        if(length(offspringSubTree)>1){
        # but it should always be the case... Test needed
          
          for(j in 1:(length(offspringSubTree))){
            
            if(offspringSubTree[j]<=n){
            # ie offspringSubTree[j] is a tip
              
              ty=as.integer(a$tip.label)
              ty=ty[which(!(ty %in% as.integer(a$tip.label[offspringSubTree[j]])))]
              ty=type[ty]                                  # type of tips in a that are not offspring[j]
              if(!(type[as.integer(a$tip.label[offspringSubTree[j]])] %in% ty)){
                # offspringSubTree[j] is a species
                species[as.integer(a$tip.label[offspringSubTree[j]])]=newSpecies
                newSpecies=newSpecies+1
              }
              
            }else{
            # ie offspringSubTree[j] is an internal node
              a1=extract.clade(a,offspringSubTree[j])
              ty=as.integer(a$tip.label)                  
              ty=ty[which(!(ty %in% as.integer(a1$tip.label)))]
              ty=type[ty]                                # type of tips in a but not in a1
              if(length(intersect(type[as.integer(a1$tip.label)],ty))==0){
                # the tips in a1 are in a different species than the other tips in a
                nextNode=c(nextNode,offspring[j])        # new candidate node
                species[as.integer(a1$tip.label)]=newSpecies
                newSpecies=newSpecies+1
              }
            }
          }
        }else{print("error : one node has only one offspring... ???")}
      }else{print(paste("error : this node should not be in next node :",i))}
    }   # end of main loop
    
    # rename the species so that there is no gap in the names
    M=rep(0,max(species))
    name=1
    for(i in 1:length(species)){
      if(M[species[i]]==0){
        M[species[i]]=name
        name=name+1
      }
      species[i]=M[species[i]]
    }
    
    return(species)
  }
  
  # auxiliary function to build the species tree from the genealogy and the species of each tip
  make.phylo=function(gen,spec){
    
    N=length(gen$tip.label)                                                                             # number of individuals
    abundance=sapply(unique(spec),function(i){sum(spec==i)})                                            # abundance of each species
    mean.trait=sapply(unique(spec),function(i){sapply(1:D,function(e){mean(gen$x.tip[e,spec[gen$tip.label]==i])})})    # mean trait of each species
    if(D == 1){mean.trait = matrix(mean.trait,nrow = D)}
    
    if(max(spec)==1){return(list(abundance=abundance,mean.trait=mean.trait))}
    
    keep.tip=c()                                                                                        # which tips will be in the species tree
    keep.tip=sapply(unique(spec),function(i){vect=(1:N)[spec[gen$tip.label]==i];
      if(length(vect)==1){return(vect)};
      sample(vect,1)})

    tree=prune.with.traits(gen,gen$tip.label[-keep.tip],gen$x)               # suppress the other tips
    tree$abundance=abundance[spec[as.integer(tree$tree$tip.label)]]
    tree$mean.trait=matrix(mean.trait[,spec[as.integer(tree$tree$tip.label)]], nrow = D, byrow = F)
    tree$tree$tip.label=spec[as.integer(tree$tree$tip.label)]
    return(tree)
  }
  
  P=aux(genealogy$P,distanceP)
  H=aux(genealogy$H,distanceH)
  Pphylo=make.phylo(genealogy$P,P)
  Hphylo=make.phylo(genealogy$H,H)
  
  } else {   # Monophyly == False
    
    # set.seed(seed) # test on 27/01/2021 to solve problem non reproducibility
    
    D=length(genealogy$P$x)       # dimention of trait space
    
    # auxiliary function to define the species for one type
    aux = function(gen,distance){
      N=length(gen$tip.label)    # number of individuals
      
      # first we define the genetic types of the individuals
      
      if(threshold==1){
        type=rep(1,N)             
        new.type=2
        for(i in 1:nrow(gen$edge)){
          if(gen$nMut[i]>0){
            if(gen$edge[i,2]<=N){
              type[gen$tip.label[gen$edge[i,2]]]=new.type
              new.type=new.type+1
            }else{
              type[as.integer(extract.clade(gen,gen$edge[i,2])$tip.label)]=new.type}
            new.type=new.type+1
          }
        }
      }else{ #threshold >1
        if(is.null(distance)){
          # we compute the distance matrix if it is not given as an argument
          distance=matrix(0,nrow=N,ncol=N)
          for(i in 1:nrow(gen$edge)){
            if(gen$nMut[i]>0){
              #add the number of mutation between the individuals separated by this edge
              if(gen$edge[i,2]<=N){
                tip=as.integer(gen$tip.label[gen$edge[i,2]])
                distance[tip,(1:N)[-tip]]=distance[tip,(1:N)[-tip]]+gen$nMut[i]
                distance[(1:N)[-tip],tip]=distance[(1:N)[-tip],tip]+gen$nMut[i]
              }else{
                tip=as.integer(extract.clade(gen,gen$edge[i,2])$tip.label)
                distance[tip,(1:N)[-tip]]=distance[tip,(1:N)[-tip]]+gen$nMut[i]
                distance[(1:N)[-tip],tip]=distance[(1:N)[-tip],tip]+gen$nMut[i]
              }
            }
          }
          distance=distance[as.integer(gen$tip.label),as.integer(gen$tip.label)]
        }
        
        # we use the distance matrix to define the genetic type of individuals
        phylo=gen
        phylo$edge.length=phylo$nMut
        type=1:N
        for(i in 1:N){
          ind=((1:i)[distance[i,1:i]<threshold])[1]
          type[phylo$tip.label[i]]=type[phylo$tip.label[ind]]

        }
      }
      
      # rename the species (=type) so that there is no gap in the names
      M=rep(0,max(type))
      name=1
      for(i in 1:length(type)){
        if(M[type[i]]==0){
          M[type[i]]=name
          name=name+1
        }
        type[i]=M[type[i]]
      }
      
      return(type)
    }
    
    # auxiliary function to build the species tree from the genealogy and the species of each tip
    make.phylo=function(gen,spec){
      
      N=length(gen$tip.label)  
      list_species <- unique(spec)
      unique_type <- c()
      pruned_gen <- gen
      pruned_gen$tip.label <- as.character(pruned_gen$tip.label)
      pruned_gen <- multi2di(pruned_gen)
      
      
      if (!is.null(seed)) { # test on 27/01/2021
        set.seed(seed)
        sampled_tips <- sample(1:N)} else {sampled_tips <- 1:N}
      
      for (i in sampled_tips){
        if (spec[i] %in% unique_type){
            pruned_gen <- drop.tip(pruned_gen,tip=as.character(i), trim.internal = TRUE)
        }else{unique_type <- c(unique_type, spec[i])}
      }
      
      pruned_gen$tip.label <- spec[as.numeric(pruned_gen$tip.label)]
      
      abundance=sapply(unique(spec),function(i){sum(spec==i)})   # abundance of each species
      mean.trait=sapply(unique(spec),function(i){sapply(1:D,function(e){mean(gen$x.tip[e,spec[gen$tip.label]==i])})})    # mean trait of each species
      if(D == 1){mean.trait = matrix(mean.trait,nrow = D)}
      if(max(spec)==1){return(list(abundance=abundance,mean.trait=mean.trait))}
    
      
      tree=list(tree=pruned_gen) 
      tree$abundance=abundance[tree$tree$tip.label]
      tree$mean.trait=matrix(mean.trait[,tree$tree$tip.label], nrow = D, byrow = F)
      # tree$mean.trait and tree$abundance order like tree$tree$tip.label
  
      return(tree)
    }

    P=aux(genealogy$P,distanceP)
    H=aux(genealogy$H,distanceH)
    Pphylo=make.phylo(genealogy$P,P)
    Hphylo=make.phylo(genealogy$H,H)

  }
  
  return(list(P=P,H=H,Pphylo=Pphylo,Hphylo=Hphylo))
}



######## Make a well conformed tree ################################

# phy a phylo object
# trait a list of vector of trait corresponding to each edge


rigth.order=function(phy,trait){
  
  root=unique(phy$edge[,1])[which(sapply(unique(phy$edge[,1]),function(x){!(x %in% phy$edge[,2])}))]  # root of the phylogeny
  nextNode=c(root)                                  # next node
  order=rep(0,phy$Nnode+length(phy$tip.label))      # rank of each node in the tree traversal
  isTip=c()                                         # rank of the tips in the tree traversal
  i=1                                               # current rank
  label=rep(0,length(phy$tip.label))                # tip labels
  
  # fill order, isTip and label by traversing the tree (Deep First Search ?)
  while(length(nextNode)>0){
    order[nextNode[1]]=i
    offspring=phy$edge[phy$edge[,1]==nextNode[1],2]
    
    if(length(offspring)==0) {
      # then nextNode[1] is a tip
      isTip=c(isTip,i)
      label[length(isTip)]=phy$tip.label[nextNode[1]]
    }
    
    nextNode=c(offspring,nextNode[-1])            # update next node
    i=i+1
  }
  
  # rename the nodes
  phy$edge[,1]=order[phy$edge[,1]]
  phy$edge[,2]=order[phy$edge[,2]]
  phy$tip.label=label[1:length(isTip)]   # if the length was not correct in the first place or tips where badly labeled
  
  # reorder the edges
  order=order(phy$edge[,2])
  phy$edge=phy$edge[order,]
  
  # reorder trait values
  if(! is.null(trait)){
    for(i in 1:length(trait)){
      trait[[i]]=(trait[[i]])[order]
    }
  }
  
  # reorder edge lengths
  phy$edge.length=phy$edge.length[order]
  
  # rename the nodes according to their nature (tip <= ntip, internal nodes > ntip)
  iTip=1                               # current tip name
  iNode=length(phy$tip.label)+1        # current node name
  newNames=c()                         # vector of new names
  
  for(i in 1:(phy$Nnode+length(phy$tip.label))){
    if(i %in% isTip) {
      newNames=c(newNames,iTip)
      iTip=iTip+1
    }else{
      newNames=c(newNames,iNode)
      iNode=iNode+1
    }
  }
  
  phy$edge=matrix(newNames[phy$edge],ncol=2)   # change the names in edge

  return(list(tree=phy,trait=trait))
}



######## Compute number of mutation between tips ###################

# gen an object obtained with make.gen


compute.dist=function(gen, verbose=T){
  
  N=length(gen$tip.label)               # number of tips
  distance=matrix(0,nrow=N,ncol=N)      # initiate distance matrix
  
  #main loop
  for(i in 1:nrow(gen$edge)){
    
    if(gen$nMut[i]>0){
      #add the number of mutation between the individuals separated by this edge

      if(gen$edge[i,2]<=N){
        tip=as.integer(gen$tip.label[gen$edge[i,2]])
        distance[tip,(1:N)[-tip]]=distance[tip,(1:N)[-tip]]+gen$nMut[i]
        distance[(1:N)[-tip],tip]=distance[(1:N)[-tip],tip]+gen$nMut[i]
      }else{
        tip=as.integer(extract.clade(gen,gen$edge[i,2])$tip.label)
        distance[tip,(1:N)[-tip]]=distance[tip,(1:N)[-tip]]+gen$nMut[i]
        distance[(1:N)[-tip],tip]=distance[(1:N)[-tip],tip]+gen$nMut[i]
      }
    }
  }  # end main loop
  
  distance=distance[as.integer(gen$tip.label),as.integer(gen$tip.label)]
  return(distance)
}



######## Remove taxa when there are traits #########################

# phy a phylo object
# nodes the labels of the tips to suppress 
# traits a list of vector of traits corresponding to the edges of the tree
# additiveTrait vector of indice of extensive traits (eg branch length, nmut...)


prune.with.traits=function(phy,nodes,traits,extensiveTrait=c()){
  
  D=length(traits)                   # dimension of trait space
  obj=list(tree=phy,trait=traits)    # returned object
  edges=c()                          # removed edges
  if(length(extensiveTrait)>0){
    intensiveTrait=(1:D)[-extensiveTrait]
  }else{
    intensiveTrait=1:D
  }
  
  if(length(nodes)>0){
    # main loop
    for(i in 1:length(nodes)){
      edge=which(obj$tree$edge[,2]==which(obj$tree$tip.label==nodes[i]))   # edge leading to nodes[i]
      
      while(edge %in% edges){
      # ie edge has already been suppressed (all the sister clades of nodes[i] have been suppressed)
        # we have to suppress the parent edge
        edge=which(obj$tree$edge[,2]==obj$tree$edge[edge,1])
      }
      
      sisterEdge=which(obj$tree$edge[,1]==obj$tree$edge[edge,1])  # sister edges of edge (including edge)
      sisterEdge=sisterEdge[which(!(sisterEdge %in% edges))]      # remaining sister edges of edge
      
      if(length(sisterEdge)<3){
        # in that case we suppress edge and its sister edge
        edges=c(edges,sisterEdge)
        parent=which(obj$tree$edge[,2]==obj$tree$edge[edge,1])
        
        if(length(extensiveTrait)>0){
          # for an extensive trait we have to sum the trait of the two combined edges
          for(j in extensiveTrait){
            obj$trait[[j]][parent]=obj$trait[[j]][sisterEdge[sisterEdge!=edge]]+obj$trait[[j]][parent]
          }
        }
        
      }else{
        # in that case we only have to suppress edge, no edges are combined
        edges=c(edge,edges)
      }
    }  # end main loop

    for(j in 1:D){
      obj$trait[[j]]=obj$trait[[j]][-edges]
    }
    
    obj$tree=drop.tip(obj$tree,which(obj$tree$tip.label %in% nodes))
  }
  return(obj)
}



######## Plot the result of the model ############################

# gen an oject returned by make.gen
# spec the corresponding object returned by define.species
#


plot.model=function(gen,spec,trait.id,lwdgen=1,lwdsp=lwdgen,scale=NULL){
  
  # cut the plot window
  layout(mat = matrix(c(1,1,8,2,2,3,3,8,4,4,5,7,7,7,6),ncol=5,byrow = T),
         widths = c(3,2,1.5,2,3), heights = c(3,3,2))
  
  # define the scale range
  if(is.null(scale)) {
    scale=range(c(gen$P$x[[trait.id]],gen$H$x[[trait.id]]))
  }else{
    scale=range(c(scale,gen$P$x[[trait.id]],gen$H$x[[trait.id]]))
  }
  if(scale[1]==scale[2]){scale[2]=scale[1]+1}
  
  # plot the genealogies
  #for P
  plot.with.trait(gen$P,gen$P$x[[trait.id]],lwd=lwdgen,scale=scale)
  title(main=paste("P",gen$P$ini))
  # for H
  plot.with.trait(gen$H,gen$H$x[[trait.id]],lwd=lwdgen,scale=scale)
  title(main=paste("H",gen$H$ini))
  
  # plot the species trees
  # for P
  if(!is.null(spec$Pphylo$tree)){
    plot.with.trait(spec$Pphylo$tree,spec$Pphylo$trait[[trait.id]],lwd=lwdsp,scale=scale)
    title(main=length(spec$Pphylo$tree$tip.label))
  }else{plot(1,spec$Pphylo$mean.trait[trait.id])}
  # for H
  if(!is.null(spec$Hphylo$tree)){
    plot.with.trait(spec$Hphylo$tree,spec$Hphylo$trait[[trait.id]],lwd=lwdsp,scale=scale)
    title(main=length(spec$Hphylo$tree$tip.label))
  }else{plot(1,spec$Hphylo$mean.trait[trait.id])}
  
  # plot the trait densities
  densplot(mcmc(gen$P$x.tip[trait.id,]))
  densplot(mcmc(gen$H$x.tip[trait.id,]))
  
  # plot the traits with their correlation 
  plot(gen$P$x.tip[trait.id,],sapply(gen$P$tip.label,function(i){gen$H$x.tip[trait.id,gen$H$tip.label==i]}),xlab="P trait",ylab="H trait")
  lines(c(scale[1],scale[2]),c(scale[1],scale[2]),type = 'l',col="red")
  cor=try(cor(gen$P$x.tip[trait.id,],sapply(gen$P$tip.label,function(i){gen$H$x.tip[trait.id,gen$H$tip.label==i]})))
  if (inherits(cor,"try_error")){
    title(main="_") 
  }else{
    title(main=cor)  
  }
  
  # add legend
  plot.new()
  rasterImage(as.raster(colorRampPalette(rev(c("steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4")))( 100 )),0,0,0.5,1)
  text(x=0.7, y = seq(0,1,l=5), labels = round(seq(scale[1],scale[2],l=5),digits = 2))
  axis(side=4,at=seq(0,1,l=5),pos=0.5,labels = F)
  
}



######## Plot a phylogeny with colored branches ##############################

plot.with.trait=function(phylo,rate,scale=NULL,lwd=1,direction="rightwards"){
  
  # define scale
  scale=range(c(scale,rate))
  if(scale[1]==scale[2]){scale[2]=scale[1]+1}
  
  # define color palette
  Colors = colorRampPalette(c("steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4"))( 100 ) 
  
  # define branch colors
  col = round( (rate - min(scale)) / diff(range(scale))*99   )+1
  
  plot(phylo, edge.color = Colors[col], show.tip.label = F,edge.width =lwd,direction=direction)
}



######## Plot the colored spatial matrix ################################

spatial.plot=function(out,trait.id,scale=NULL,nx=NULL, sort_trait = F){
  if(is.null(nx))   nx=sqrt(length(out$P[1,]))
  scale=range(c(scale,out$P[trait.id,],out$H[trait.id,]))
  if(scale[1]==scale[2]){
    scale[1]=scale[1]-0.5
    scale[2]=scale[2]+0.5
  }
  MP=max(out$P[trait.id,])
  mP=min(out$P[trait.id,])
  MH=max(out$H[trait.id,])
  mH=min(out$H[trait.id,])
  Colors = colorRampPalette(c("steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4"))( 100 ) 
  par(mfrow=c(1,2))
  if (sort_trait){
    order_P = order(out$P[trait.id,])
    image(matrix(out$P[trait.id,][order_P],ncol=nx,byrow = T),col=Colors[(round( (mP - min(scale)) / diff(range(scale))*99   )+1):(round( (MP - min(scale)) / diff(range(scale))*99   )+1)],
          main="",axes=F)
    image(matrix(out$H[trait.id,][order_P],ncol=nx,byrow = T),col=Colors[(round( (mH - min(scale)) / diff(range(scale))*99   )+1):(round( (MH - min(scale)) / diff(range(scale))*99   )+1)]
          ,main="",axes=F)
  }else{
    image(matrix(out$P[trait.id,],ncol=nx,byrow = T),col=Colors[(round( (mP - min(scale)) / diff(range(scale))*99   )+1):(round( (MP - min(scale)) / diff(range(scale))*99   )+1)],
          main="",axes=F)
    image(matrix(out$H[trait.id,],ncol=nx,byrow = T),col=Colors[(round( (mH - min(scale)) / diff(range(scale))*99   )+1):(round( (MH - min(scale)) / diff(range(scale))*99   )+1)]
          ,main="",axes=F)
  }
}



######## Build the network ##################

build.network=function(xP, xH, gen, spec, alpha, seuil=0, rP=2, rH=2 , f=3){
  
  if(f==1){
    f=function(xP,xH,alpha,r){(1+sum((xP-xH)^2))^(-alpha)}             # the default
  }else if(f==2){
    f=function(xP,xH,alpha,r){if(alpha<0){exp(sum(((xP-xH)*alpha)^2))}else{0.00001+exp(-sum(((xP-xH)*alpha)^2))}}
  }else if(f==3){
    f=function(xP,xH,alpha,r){if(alpha<0){1/(r-1)+1-exp(-sum(((xP-xH)*alpha)^2))}else{1/(r-1)+exp(-sum(((xP-xH)*alpha)^2))}}    
  }else if(f==4){
    f=function(xP,xH,alpha,r){if(sum((xP-xH)^2)<R){exp(alpha)}else{1}}
  }
  
  P=spec$P                                         # what Pspecies is each P individual in ?
  H=spec$H                                         # what Hspecies is each H individual in ?
  N=length(P)                                      # number of individual in each clade
  nSP=max(P)                                       # number of Pspecies
  nSH=max(H)                                       # number of Hspecies
  link=Matrix(0,nrow = nSP,ncol = nSH,sparse = T)  # initialize the link matrix
  finfP=1/(rP-1)+(alphaP<0)
  finfH=1/(rH-1)+(alphaH<0)
  
  for(i in 1:N){
    fP=abs(f(xP[,i],xH[,i],alphaP,rP)-finfP)
    fH=abs(f(xH[,i],xP[,i],alphaH,rH)-finfH)
    if(max(fP,fH)>=seuil){
      link[P[i],H[i]]=link[P[i],H[i]]+1
    }
  }
  
  return(link)
}



######## Plot the network ##################################

plot.network=function(link,spec,trait.id=1,method="bipartite",order=T,scale=c()){

  if(method=="bipartite"){
    Colors = colorRampPalette(c("steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4"))( 100 )
    scale=range(c(scale,spec$Pphylo$mean.trait[trait.id,],spec$Hphylo$mean.trait[trait.id,]))
    col.P=rep(0,length(spec$Pphylo$tree$tip.label))
    col.H=rep(0,length(spec$Hphylo$tree$tip.label))
    col.P[spec$Pphylo$tree$tip.label]=Colors[(round( (spec$Pphylo$mean.trait[trait.id,] - min(scale)) / diff(range(scale))*99   )+1)]
    col.H[spec$Hphylo$tree$tip.label]=Colors[(round( (spec$Hphylo$mean.trait[trait.id,] - min(scale)) / diff(range(scale))*99   )+1)]
    try(plotweb(as.matrix(link),col.low=col.P,col.high = col.H,col.interaction = "lightgray",empty=T))
  }else{
    Mat=as.matrix(link)
    if(order) Mat=Mat[order(rowSums(Mat),decreasing = T),order(colSums(Mat),decreasing = F)]
    image(log(Mat+1), axes = FALSE, col = grey(seq(1, 0.2, length = 256)))
    title(xlab = "P", ylab="H",outer = F,line=1)
  }
}



######## Plot the results with the network ###############

plot.model.network=function(gen,spec,trait.id, link,out,lwdgen=1,lwdsp=lwdgen,scale=NULL,nx=NULL,cor=F,network.method="bipartite",spatial=T, sort_trait = F){
  
  # cut the plot window
  if(spatial){
    layout(mat = matrix(c(1,1,8,2,2,2,2,3,3,8,4,4,4,4,5,7,7,7,6,6,6,9,9,9,10,10,11,11),ncol=7,byrow = T),
         widths = c(3,2,1.5,2,0.5,0.5,2), heights = c(3,3,2.5,4))
  }else{
    layout(mat = matrix(c(1,1,8,2,2,2,2,3,3,8,4,4,4,4,5,7,7,7,6,6,6,9,9,9,9,9,9,9),ncol=7,byrow = T),
           widths = c(3,2,1.5,2,0.5,0.5,2), heights = c(3,3,2.5,4))
  }
  
  # define the scale range
  if(is.null(scale)) {
    scale=range(c(gen$P$x[[trait.id]],gen$H$x[[trait.id]]))
  }else{
    scale=range(c(scale,gen$P$x[[trait.id]],gen$H$x[[trait.id]]))
  }
  if(scale[1]==scale[2]){scale[2]=scale[1]+1}
  
  # plot the genealogies
  #for P
  plot.with.trait(gen$P,gen$P$x[[trait.id]],lwd=lwdgen,scale=scale)
  title(main=paste("P",gen$P$ini))
  # for H
  plot.with.trait(gen$H,gen$H$x[[trait.id]],lwd=lwdgen,scale=scale)
  title(main=paste("H",gen$H$ini))
  
  # plot the species trees
  # for P
  if(!is.null(spec$Pphylo$tree)){
    plot.with.trait(spec$Pphylo$tree,spec$Pphylo$trait[[trait.id]],lwd=lwdsp,scale=scale)
    title(main=length(spec$Pphylo$tree$tip.label))
  }else{plot(1,spec$Pphylo$mean.trait[trait.id])}
  # for H
  if(!is.null(spec$Hphylo$tree)){
    plot.with.trait(spec$Hphylo$tree,spec$Hphylo$trait[[trait.id]],lwd=lwdsp,scale=scale)
    title(main=length(spec$Hphylo$tree$tip.label))
  }else{plot(1,spec$Hphylo$mean.trait[trait.id])}
  
  # plot the trait densities
  densplot(mcmc(gen$P$x.tip[trait.id,]))
  densplot(mcmc(gen$H$x.tip[trait.id,]))
  
  # plot the traits with their correlation 
  if(cor){
    plot(gen$P$x.tip[trait.id,],sapply(gen$P$tip.label,function(i){gen$H$x.tip[trait.id,gen$H$tip.label==i]}),xlab="P trait",ylab="H trait")
    lines(c(scale[1],scale[2]),c(scale[1],scale[2]),type = 'l',col="red")
  }else{
    if(is.null(spec)){
      plot(1,1)
    }else{
    plot.species(spec,cex=2,xlab=c(),ylab=c(), net = link)}
  }
  cor=try(cor(gen$P$x.tip[trait.id,],sapply(gen$P$tip.label,function(i){gen$H$x.tip[trait.id,gen$H$tip.label==i]})))
  if (inherits(cor,"try_error")){
    title(main="_") 
  }else{
    title(main=cor)  
  }
  
  # add legend
  plot.new()
  rasterImage(as.raster(colorRampPalette(rev(c("steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4")))( 100 )),0,0,0.5,1)
  text(x=0.7, y = seq(0,1,l=5), labels = round(seq(scale[1],scale[2],l=5),digits = 2))
  axis(side=4,at=seq(0,1,l=5),pos=0.5,labels = F)
  
  # plot.new()
  t=try(plot.network(link,spec,trait.id,method=network.method))
  if(inherits(t,"try-error")) plot.new()
  
  if(spatial){
    if(is.null(nx))   nx=sqrt(length(out$P[1,]))
    scale=range(c(scale,out$P[trait.id,],out$H[trait.id,]))
    MP=max(out$P[trait.id,])
    mP=min(out$P[trait.id,])
    MH=max(out$H[trait.id,])
    mH=min(out$H[trait.id,])
    Colors = colorRampPalette(c("steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4"))( 100 ) 
    if (sort_trait){
      order_P = order(out$P[trait.id,])
      image(matrix(out$P[trait.id,][order_P],ncol=nx,byrow = T),col=Colors[(round( (mP - min(scale)) / diff(range(scale))*99   )+1):(round( (MP - min(scale)) / diff(range(scale))*99   )+1)],
            main="",axes=F)
      image(matrix(out$H[trait.id,][order_P],ncol=nx,byrow = T),col=Colors[(round( (mH - min(scale)) / diff(range(scale))*99   )+1):(round( (MH - min(scale)) / diff(range(scale))*99   )+1)]
            ,main="",axes=F)
    }else{
      image(matrix(out$P[trait.id,],ncol=nx,byrow = T),col=Colors[(round( (mP - min(scale)) / diff(range(scale))*99   )+1):(round( (MP - min(scale)) / diff(range(scale))*99   )+1)],
            main="",axes=F)
      image(matrix(out$H[trait.id,],ncol=nx,byrow = T),col=Colors[(round( (mH - min(scale)) / diff(range(scale))*99   )+1):(round( (MH - min(scale)) / diff(range(scale))*99   )+1)]
            ,main="",axes=F)
    }
    }
  
  
}



######## Plot species in the phenotype space ##############



plot.species=function(spec,trait.id=1:3,net = NULL,...){
  D = max(1,nrow(spec$Pphylo$mean.trait))
  trait.id[trait.id>D] = D
  gr = "gray50"
  Colors = colorRampPalette(c("steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4"))( 100 )
  scale=range(c(spec$Pphylo$mean.trait[trait.id[3],],spec$Hphylo$mean.trait[trait.id[3],]))
  col.P=Colors[(round( (spec$Pphylo$mean.trait[trait.id[3],] - min(scale)) / diff(range(scale))*99   )+1)]
  col.H=Colors[(round( (spec$Hphylo$mean.trait[trait.id[3],] - min(scale)) / diff(range(scale))*99   )+1)]
  
  scale_x=range(c(spec$Pphylo$mean.trait[trait.id[1],],spec$Hphylo$mean.trait[trait.id[1],]))
  scale_y=range(c(spec$Pphylo$mean.trait[trait.id[2],],spec$Hphylo$mean.trait[trait.id[2],]))
  plot(c(),c(),xlim=scale_x,ylim=scale_y)
  if(!is.null(net)){
    for (i in 1:nrow(net)){
      for (j in 1:ncol(net)){
        if (net[i,j]>0){
          lines(c(spec$Pphylo$mean.trait[trait.id[1],i],spec$Hphylo$mean.trait[trait.id[1],j]),
                c(spec$Pphylo$mean.trait[trait.id[2],i],spec$Hphylo$mean.trait[trait.id[2],j]),
                col = gr)
        }
      }
    }
  }
  points(spec$Pphylo$mean.trait[trait.id[1],],spec$Pphylo$mean.trait[trait.id[2],],pch=16,col=col.P,...)
  points(spec$Hphylo$mean.trait[trait.id[1],],spec$Hphylo$mean.trait[trait.id[2],],pch=18,col=col.H,...)
}