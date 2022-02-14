
#library('treebase')

proba_split_non_norm=function(beta,n,i){
  if(beta>-2){
    return(gamma(beta+i+1)*gamma(beta+n-i+1)/(gamma(i+1)*gamma(n-i+1)))
  }else{
    return(1*(i==1)+1*(i==(n-1)))
  }
}

proba_split_non_norm_int=function(beta,n,i,bound,silent=T){
  if(beta>-2){
    b=bound
    continue=TRUE
    aux=function(x){
      (x^(i+beta))*((1-x)^(n-i+beta))
    }
    while(continue & b>0){
    t=try(choose(n,i)*integrate(aux,lower=10^(-1*b),upper=1-10^(-1*b))$value, silent=T)
    if(inherits(t, "try-error")){
      if (silent==F){
      print(paste("beta =", beta,"n =",n,"i =",i,"b =",b))}
      b=b-1
    }else{
      continue=FALSE
      }
    }
    return(t)
  }else{
    return(1*(i==1)+1*(i==(n-1)))
  }
}

norm_cte=function(beta,n){
  i=1:(n-1)
  aux=function(j){
    return(proba_split_non_norm(beta,n,i=j))
  }
  return(sum(sapply(i,aux)))
}

norm_cte_int=function(beta,n,silent=T){
  i=1:(n-1)
  aux=function(j){
    return(proba_split_non_norm_int(beta,n,i=j,10,silent))
  }
  return(sum(sapply(i,aux)))
}

proba_split=function(beta,n,i){
  if(beta<(-2)){
    return(-Inf)
  }else{
  return(proba_split_non_norm(beta,n,i)/norm_cte(beta,n))
  }
}

proba_split_int=function(beta,n,i,silent=T){
  a=n
  b=i
  c=beta
  if(beta<(-2)){
    return(-Inf)
  }else{
    return(proba_split_non_norm_int(beta,n,i,10,silent)/norm_cte_int(beta,n,silent))
  }
}

split=function(arbre){
  n=arbre$Nnode+1
  rep=matrix(0,3,2*n-1)
  rep[1,1:n]=rep(1,n)
  for (i in (2*n-2):1){
#     a=arbre$edge[i,1]
#     b=arbre$edge[i,2]
#     a1=rep[1,arbre$edge[i,1]]
#     a2=rep[1,arbre$edge[i,2]]
#     rep[1,arbre$edge[i,1]]=a1+a2
    rep[1,arbre$edge[i,1]]=rep[1,arbre$edge[i,1]]+rep[1,arbre$edge[i,2]]
#     if (rep[2,arbre$edge[i,1]]==0){
#       rep[2,arbre$edge[i,1]]=rep[1,arbre$edge[i,2]]
#     }else{
#       rep[3,arbre$edge[i,1]]=rep[1,arbre$edge[i,2]]
#     }
      if (rep[3,arbre$edge[i,1]]==0){
        rep[3,arbre$edge[i,1]]=rep[1,arbre$edge[i,2]]
      }else{
        rep[2,arbre$edge[i,1]]=rep[1,arbre$edge[i,2]]
      }
  }
  return(rep)
}

show.splits=function (arbre,n,text, adj = c(0.5, 0.5), frame = "rect", pch = NULL, 
                         thermo = NULL, pie = NULL, piecol = NULL, col = "black", 
                         bg = "lightblue", horiz = FALSE, width = NULL, height = NULL, 
                         ...) 
{
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  splits=split(arbre)[1,(lastPP$Ntip + 1):length(lastPP$xx)]
  nodes <- (lastPP$Ntip + 1):length(lastPP$xx)
  XX <- lastPP$xx[nodes]
  YY <- lastPP$yy[nodes]
  BOTHlabels(text, splits, XX, YY, adj, frame, pch, thermo, pie, 
             piecol, col, bg, horiz, width, height, ...)
}

proba_arbre=function(beta,splits){
  n=ncol(splits)
  aux=function(i){
    return(proba_split(beta,splits[1,i],splits[2,i]))
  }
  i=((n+3)/2):n
  prod(sapply(i,aux))
}

proba_arbre_int=function(beta,splits,silent=T,log=F){
  b=beta
  n=ncol(splits)
  if(log){
    aux=function(i){
      return(log(proba_split_int(beta,splits[1,i],splits[2,i],silent)))
    }
    i=((n+3)/2):n
    sum(sapply(i,aux))
  }
  else{
    aux=function(i){
      return(proba_split_int(beta,splits[1,i],splits[2,i],silent))
    }
    i=((n+3)/2):n
    prod(sapply(i,aux))
  }
}



log_proba_arbre=function(beta,splits){
  n=ncol(splits)
  aux=function(i){
    return(log(proba_split(beta,splits[1,i],splits[2,i])))
  }
  i=((n+3)/2):n
  sum(sapply(i,aux))
}

mle_beta=function(arbre){
  prit("c'est reparti")
  splits=split(arbre)
  aux=function(beta){
    proba_arbre(beta,splits)
  }
  optimize(aux, c(-2,10), maximum = TRUE)$maximum
}

mle_beta_int=function(arbre,bmax,silent=T,log=T){
  splits=split(arbre)
  aux=function(beta){
    proba_arbre_int(beta,splits,silent,log)
  }
  optimize(aux, c(-2,bmax), maximum = TRUE)$maximum
}

log_mle_beta=function(arbre,bmax=10){
  splits=split(arbre)
  aux=function(beta){
    log_proba_arbre(beta,splits)
  }
  optimize(aux, c(-2,bmax), maximum = TRUE)$maximum
}

plot_beta=function(arbre,bmin,bmax,pas,log=T){
  x=seq(bmin,bmax,pas)
  spl=split(arbre)
  aux=function(b){
    proba_arbre_int(b,spl,log=log)
  }
  y=sapply(x,aux)
  #print(y)
    plot(y~x,type='l')
}

#######################test##########################################
# a=arbres.avec.tps[[37]]
# plot_beta(a,525,535,0.0001)
# 
# plot_beta(a,-2,2,0.01,T)
# proba_split(-2,7,1)
# undebug(mle_beta)
# load("arbres_20_1000")
# arbre=arbre.final[[1082]]
# #str(subtrees(arbre)[[5]])
# split(arbre)
# plot(arbre)
# mle_beta(arbre)
# undebug(proba_arbre)
# proba_arbre(-1,split(arbre))
