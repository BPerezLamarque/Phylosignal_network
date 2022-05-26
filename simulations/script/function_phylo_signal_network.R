pblm <- function(assocs,tree1=NULL,tree2=NULL,covars1=NULL,covars2=NULL,bootstrap=FALSE,nreps=10,maxit=10000,pstart=c(.5,.5)){
  
  # Make a vector of associations
  A<-as.matrix(as.vector(as.matrix(assocs)))
  data.vecs<-A
  
  #numbers of species and interactions
  nassocs<-length(A)
  nspp1<-dim(assocs)[1]
  nspp2<-dim(assocs)[2]
  sppnames1<-rownames(assocs)
  sppnames2<-colnames(assocs)
  #make names of species pairs
  pairnames=NULL  # make a vector of pairwise comparison names
  for (o in 1:(nspp2)) {
    for (u in 1:nspp1) {
      pairnames<-c(pairnames,paste(sppnames2[o],sppnames1[u],sep="-"))
    }
  }
  
  #Clean Covariates
  #If the covariate applies to both sets, then it should be in the matrix of the longer set 
  covnames<-NULL
  C1covs<-NULL
  if(is.null(covars1)) {
    C1<-NULL } else {
      if(is.null(dim(covars1))) {
        C1<-matrix(covars1,nspp1,nspp2,byrow=FALSE)
        if(is.factor(covars1)) {
          C1<-as.matrix(as.vector(C1))
          C1covs<-cbind(C1covs,C1)
          C1<-as.matrix(model.matrix(~as.factor(C1)-1)[,-1])
          colnames(C1)<-paste(rep("covar1",length(levels(covars1))-1),levels(covars1)[-1],sep="-")
        } else {
          C1<-as.matrix(as.vector(C1))
          C1covs<-cbind(C1covs,C1)
        }
        covnames<-c(covnames,"covar1")
      } else {
        C1<-NULL
        for(i in 1:dim(covars1)[2]) {
          C1hold<-matrix(covars1[,i],nspp1,nspp2,byrow=FALSE)
          if(is.factor(covars1[,i])) {
            C1hold<-as.matrix(as.vector(C1hold))
            C1covs<-cbind(C1covs,C1hold)
            C1hold<-as.matrix(model.matrix(~as.factor(C1hold)-1)[,-1])
            colnames(C1hold)<-paste(rep(colnames(covars1)[i],length(levels(covars1[,i]))-1),levels(covars1[,i])[-1],sep="-")
            C1<-cbind(C1,C1hold)
          } else { 
            C1hold<-as.matrix(as.vector(C1hold))
            C1covs<-cbind(C1covs,C1hold)
            colnames(C1hold)<-colnames(covars1)[i]
            C1<-cbind(C1,C1hold)
          }
          covnames<-c(covnames,colnames(covars1)[i])
        }
      }
      data.vecs<-cbind(data.vecs,C1covs)
    }
  
  C2covs<-NULL
  if(is.null(covars2)) {
    C2<-NULL } else {
      if(is.null(dim(covars2))) {
        C2<-matrix(covars2,nspp1,nspp2,byrow=TRUE)
        if(is.factor(covars2)) {
          C2<-as.matrix(as.vector(C2))
          C2covs<-cbind(C2covs,C2)
          C2<-as.matrix(model.matrix(~as.factor(C2)-1)[,-1])
          colnames(C2)<-paste(rep("covar2",length(levels(covars2))-1),levels(covars2)[-1],sep="-")
        } else {
          C2<-as.matrix(as.vector(C2))
          C2covs<-cbind(C2covs,C2)
        }
        covnames<-c(covnames,"covar2")
      } else {
        C2<-NULL
        for(i in 1:dim(covars2)[2]) {
          C2hold<-matrix(covars2[,i],nspp1,nspp2,byrow=TRUE)
          if(is.factor(covars2[,i])) {
            C2hold<-as.matrix(as.vector(C2hold))
            C2covs<-cbind(C2covs,C2hold)
            C2hold<-as.matrix(model.matrix(~as.factor(C2hold)-1)[,-1])
            colnames(C2hold)<-paste(rep(colnames(covars2)[i],length(levels(covars2[,i]))-1),levels(covars2[,i])[-1],sep="-")
            C2<-cbind(C2,C2hold)
          } else { 
            C2hold<-as.matrix(as.vector(C2hold))
            C2covs<-cbind(C2covs,C2hold)
            colnames(C2hold)<-colnames(covars2)[i]
            C2<-cbind(C2,C2hold)
          }
          covnames<-c(covnames,colnames(covars2)[i])
        }
      }
      data.vecs<-cbind(data.vecs,C2covs)
    }
  
  
  # Make U, the combined matrix of covariates 
  U<-NULL
  if(is.null(C1) & is.null(C2)) {
    U<-rep(1,length(A))
  } else {
    if(is.null(C1)) {
      U<-rep(1,length(A))
    } else {
      U<-cbind(rep(1,length(A)),C1)
    }
    
    if(is.null(C2)) {
      U<-U
    } else {
      U<-cbind(U,C2)
    }
  }
  
  # Begin to organize output
  if(is.null(dim(U)))
  {
    data.vecs<-data.frame(A)
    colnames(data.vecs)<-"associations"
  } else {    
    colnames(data.vecs)<-c("associations", covnames)
  }
  rownames(data.vecs)<-pairnames
  
  ######
  # Calculate Star Regression Coefficients
  #calculate for the star (assuming no phylogenetic correlation)
  astar<-solve((t(U)%*%U),(t(U)%*%A))
  MSETotal<-cov(A)
  s2aStar<-as.vector(MSETotal)*chol2inv(chol((t(U)%*%U)))
  
  sdaStar<-t(diag(s2aStar)^(.5))
  approxCFstar<-rbind(t(astar)-1.96%*%sdaStar, t(astar), t(astar)+1.96%*%sdaStar)
  Pstar<-U%*%astar
  Estar<-A-Pstar
  MSEStar<-cov(matrix(Estar))
  
  #######
  if(is.null(tree1) | is.null(tree2)) {
    coefs<-approxCFstar
    rownames(coefs)<-c("lower CI 95%","estimate","upper CI 95%")
    colnames(coefs)<-paste("star",c("intercept",colnames(U)[-1]),sep="-")
    MSEs<-cbind(data.frame(MSETotal),data.frame(MSEStar))
    Pstar<-data.frame(Pstar)
    colnames(Pstar)<-"star"
    Estar<-data.frame(Estar)
    colnames(Estar)<-"star"
    output<-list(MSE=MSEs,signal.strength=NULL,coefficients=data.frame(t(coefs)),CI.boot=NULL,variates=data.frame(data.vecs),residuals=Estar,predicted=Pstar,bootvalues=NULL,Vfull=NULL)
    class(output)<-"pblm"
    return(output)
    
  } else {
    
    #tree1 is the phylogeny for the rows
    #tree2 is the phylogeny for the columns
    
    #Clean Trees
    if(is(tree1)[1] %in% c("phylo","chronos")) {
      if(is.null(tree1$edge.length)){tree1<-compute.brlen(tree1, 1)}  #If phylo has no given branch lengths
      V1<-vcv.phylo(tree1,corr=TRUE)
      V1<-V1[rownames(assocs),rownames(assocs)]
    } else {
      V1<-tree1[rownames(assocs),rownames(assocs)]
    }    
    
    if(is(tree2)[1] %in% c("phylo","chronos")) {
      if(is.null(tree2$edge.length)){tree2<-compute.brlen(tree2, 1)}  #If phylo has no given branch lengths
      V2<-vcv.phylo(tree2,corr=TRUE)
      V2<-V2[colnames(assocs),colnames(assocs)]
    } else {
      V2<-tree2[colnames(assocs),colnames(assocs)]
    }
    
    #Calculate Regression Coefficents for the base (assuming strict brownian motion evolution, ds=1)
    V1<-as.matrix(V1)
    V2<-as.matrix(V2)
    
    detV1 <- det(V1)
    if (detV1==0) {detV1=1e-300}
    detV2 <- det(V2) 
    if (detV2==0) {detV2=1e-300}
    
    V1<-V1/detV1^(1/nspp1)   # scale covariance matrices (this reduces numerical problems caused by
    V2<-V2/detV2^(1/nspp2)   # determinants going to infinity or zero)
    
    # invert and then kronecker product
    invV1 <- chol2inv(chol(V1))
    invV2 <- chol2inv(chol(V2))
    invV <- kronecker(invV2,invV1)  
    
    abase<-solve((t(U)%*%invV%*%U),((t(U)%*%invV%*%A)))  
    MSEBase<-(t(A-U%*%abase)%*%invV%*%(A-U%*%abase))/(nassocs-1)  
    s2abase<-as.vector(MSEBase)*chol2inv(chol(t(U)%*%invV%*%U))
    sdabase<-t(diag(s2abase)^(.5))
    approxCFbase<-rbind(t(abase)-1.96%*%sdabase, t(abase), t(abase)+1.96%*%sdabase)
    Pbase<-t(t(U%*%abase)%*%invV)
    Ebase<-A-Pbase
    
    ###################
    # Full EGLS estimates of phylogenetic signal
    ##################
    initV1<-V1
    initV2<-V2
    
    # tau = tau_i + tau_j where tau_i equals the node to tip distance
    tau1<-matrix(diag(initV1),nspp1,nspp1) + matrix(diag(initV1),nspp1,nspp1)-2*initV1
    tau2<-matrix(diag(initV2),nspp2,nspp2) + matrix(diag(initV2),nspp2,nspp2)-2*initV2
    
    # The workhorse function to estimate ds
    pegls<-function(parameters){
      d1<-abs(parameters[1])
      d2<-abs(parameters[2])
      
      V1<-(d1^tau1)*(1-d1^(2*initV1))/(1-d1^2)
      V2<-(d2^tau2)*(1-d2^(2*initV2))/(1-d2^2)
      
      if (d1==1){V1 <- initV1} # Brownian motion
      if (d2==1){V2 <- initV2} # Brownian motion
      
      detV1 <- det(V1)
      
      if (detV1==0) {detV1=1e-300}
      detV2 <- det(V2) 
      if (detV2==0) {detV2=1e-300}
      
      V1<-V1/detV1^(1/nspp1)
      V2<-V2/detV2^(1/nspp2)
      
      if (all(!is.infinite(V1))&all(!is.infinite(V2))){
        
        if (all(V1<1e10)&all(V2<1e10)){
          
          invV1 <- chol2inv(chol(V1))
          invV2 <- chol2inv(chol(V2))
          invV <- kronecker(invV2,invV1)  
          
          a<-solve((t(U)%*%invV%*%U),((t(U)%*%invV%*%A))) 
          E<-(A-U%*%a)
          MSE=t(E)%*%invV%*%E/(nassocs-1)
          return(MSE)
          
        }else{return(1e300)}
      }else{return(1e300)}
    }
    # estimate d1 and d2 by minimizing MSE
    est<-optim(pstart,pegls,control=list(maxit=maxit, reltol=1e-4, abstol=1e-4))      
    
    MSEFull<-est$value
    d1<-abs(est$par[1])
    d2<-abs(est$par[2])
    
    # Calculate EGLS coef w estimated ds 
    V1<-(d1^tau1)*(1-d1^(2*initV1))/(1-d1^2)
    V2<-(d2^tau2)*(1-d2^(2*initV2))/(1-d2^2)
    
    detV1 <- det(V1)
    if (detV1==0) {detV1=1e-300}
    detV2 <- det(V2) 
    if (detV2==0) {detV2=1e-300}
    
    V1<-V1/detV1^(1/nspp1)
    V2<-V2/detV2^(1/nspp2)
    
    
    # test 27/12/20: V needed for bootstraping only
    V <- kronecker(V2, V1)
    
    
    invV1 <- chol2inv(chol(V1))
    invV2 <- chol2inv(chol(V2))
    invV <- kronecker(invV2,invV1)  
    aFull<-solve((t(U)%*%invV%*%U),((t(U)%*%invV%*%A)))
    s2aFull<-as.vector(MSEFull)*chol2inv(chol(t(U)%*%invV%*%U))
    sdaFull<-t(diag(s2aFull)^(.5))
    approxCFfull<-rbind(t(aFull)-1.96%*%sdaFull, t(aFull), t(aFull)+1.96%*%sdaFull)
    Pfull<-t(t(U%*%aFull)%*%invV)
    Efull<-A-Pfull
    
    ########################################
    
    #organize output
    coefs<-cbind(approxCFfull,approxCFstar,approxCFbase)
    rownames(coefs)<-c("approx lower CI 95%","estimate","approx upper CI 95%")
    colnames(coefs)<-c(paste("full",c("intercept",colnames(U)[-1]),sep="-"),paste("star",c("intercept",colnames(U)[-1]),sep="-"),paste("base",c("intercept",colnames(U)[-1]),sep="-"))
    coefs<-t(coefs)
    CI.boot<-NULL
    MSEs<-cbind(data.frame(MSETotal),data.frame(MSEFull), data.frame(MSEStar), data.frame(MSEBase))
    residuals<-cbind(data.frame(Efull),data.frame(Estar),data.frame(Ebase))
    predicted<-cbind(data.frame(Pfull),data.frame(Pstar),data.frame(Pbase))
    rownames(residuals)<-pairnames
    rownames(predicted)<-pairnames
    colnames(predicted)<-c("full","star","base")
    colnames(residuals)<-c("full","star","base")
    phylocovs=list(V1=V1,V2=V2)
    
    ################
    #bootstrap CIs
    if (bootstrap) {
      Vtrue<-V
      Atrue<-A
      atrue<-aFull
      dtrue<-c(d1,d2)
      ehold<-eigen(Vtrue,symmetric=TRUE)
      L<-ehold$vectors[,nassocs:1]    #A or L
      G<-sort(ehold$values)      #D
      iG<-diag(G^-.5)    #iD
      
      # Construct Y = TT*A so that 
      # E{(Y-b)*(Y-b)'} = E{(TT*A-b)*(TT*A-b)'}
      #				  = T*V*T'
      #				  = I
      
      TT <- iG%*%t(L)
      Y <- TT%*%Atrue
      Z <- TT%*%U
      
      res <- (Y-Z%*%atrue)	# residuals in orthogonalized space
      
      # test on 27/12/20
      #invT <- chol2inv(chol(TT)) # Error in chol.default(TT) : the leading minor of order 2 is not positive definite
      invT <- qr.solve(TT)
      
      bootlist=NULL
      for (i in 1:nreps) {
        randindex<-sample(1:nassocs,replace=TRUE)	# vector of random indices
        #randindex=1:nassocs					# retain order
        YY<-Z%*%atrue+res[randindex]	# create new values of Y with random residuals
        A<-invT%*%YY	# back-transformed data
        pstart<-dtrue+c(0,.1)
        estRand<-optim(pstart,pegls,control=list(maxit=maxit, reltol=1e-4, abstol=1e-4))      
        
        MSEFullrand<-estRand$value
        d1rand<-abs(estRand$par[1])
        d2rand<-abs(estRand$par[2])
        
        # Calculate EGLS coef w estimated ds 
        V1<-(d1rand^tau1)*(1-d1rand^(2*initV1))/(1-d1rand^2)
        V2<-(d2rand^tau2)*(1-d2rand^(2*initV2))/(1-d2rand^2)
        
        detV1 <- det(V1)
        if (detV1==0) {detV1=1e-300}
        detV2 <- det(V2) 
        if (detV2==0) {detV2=1e-300}
        
        V1<-V1/detV1^(1/nspp1)
        V2<-V2/detV2^(1/nspp2)
        
        invV1 <- chol2inv(chol(V1))
        invV2 <- chol2inv(chol(V2))
        invV <- kronecker(invV2,invV1)  
        arand<-solve((t(U)%*%invV%*%U),((t(U)%*%invV%*%A)))
        
        bootlist<-rbind(bootlist,c(d1rand, d2rand, t(arand)))
      }
      nr<-dim(bootlist)[1]
      nc<-dim(bootlist)[2]
      
      #Calculate bootstrapped CIs
      alpha<-0.05  # alpha is always 0.05, but could change here
      conf<-NULL
      for(j in 1:nc){
        bootconf<-quantile(bootlist[,j],probs = c(alpha/2, 1-alpha/2))
        conf<-rbind(conf,c(bootconf[1],bootconf[2]))
      }
      signal.strength<-data.frame(cbind(conf[1:2,1],dtrue,conf[1:2,2]))
      rownames(signal.strength)<-c("d1","d2")
      colnames(signal.strength)<-c("booted lower CI 95%","estimate","booted upper CI 95%")
      
      #organize output
      CI.boot<-conf
      rownames(CI.boot)<-c("d1","d2","intercept",colnames(U)[-1])
      colnames(CI.boot)<-c("booted lower CI 95%","booted upper CI 95%")
      colnames(bootlist)<-c("d1","d2","intercept",colnames(U)[-1])
      output<-list(MSE=MSEs,signal.strength=signal.strength,coefficients=data.frame(coefs),CI.boot=CI.boot,variates=data.frame(data.vecs),predicted=predicted,residuals=residuals,bootvalues=bootlist,phylocovs=phylocovs)
      class(output)<-"pblm"
      return(output)
      
    } else {
      ########
      # If bootstrapping not performed
      
      conf<-matrix(NA,2,2)
      signal.strength<-data.frame(cbind(conf[1,],c(d1,d2),conf[2,]))
      rownames(signal.strength)<-c("d1","d2")
      colnames(signal.strength)<-c("booted lower CI 95%","estimate","booted upper CI 95%")
      output<-list(MSE=MSEs,signal.strength=signal.strength,coefficients=data.frame(coefs),CI.boot=NULL,variates=data.frame(data.vecs),predicted=predicted,residuals=residuals,bootvalues=NULL,phylocovs=phylocovs)
      class(output)<-"pblm"
      return(output)
    }
  }                                                                                                       
}

mantel_test <- function(formula = formula(data), data = sys.parent(), correlation = "Pearson", nperm = 1000) {
  
  if (!correlation %in% c("Pearson", "Spearman", "Kendall")) {stop("\"correlation\" must be among 'Pearson', 'Spearman', or 'Kendall'.")}
  if (!is.numeric(nperm)) {stop("Please provide a numeric number of permutations (\"nperm\").")}
  
  m <- match.call(expand.dots = FALSE)
  m2 <- match(c("formula", "data"), names(m), nomatch = 0)
  m <- m[c(1, m2)]
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  m <- as.matrix(m)
  n <- (1 + sqrt(1 + 8 * nrow(m)))/2
  if (abs(n - round(n)) > 1e-12) 
    stop("Matrix not square.")
  n <- round(n)
  if (ncol(m) < 2) 
    stop("Not enough data.")
  
  ymat <- as.vector(m[, 1])
  xmat <- as.vector(m[, 2])
  
  if (correlation %in% c("Spearman", "Kendall")){
    ymat <- rank(ymat)
    xmat <- rank(xmat)
  }
  
  ycor <- ymat
  xcor <- xmat
  
  if (correlation=="Pearson"){mantelr <- cor(xcor, ycor, use = "everything", method = "pearson")}
  if (correlation=="Spearman"){mantelr <- cor(xcor, ycor, use = "everything", method = "spearman")}
  if (correlation=="Kendall"){mantelr <- cor(xcor, ycor, use = "everything", method = "kendall")}
  
  xmat <- full(xmat)
  ymat <- full(ymat)
  xmat <- xmat[col(xmat) > row(xmat)]
  ymat <- ymat[col(ymat) > row(ymat)]
  if (nperm > 0) {
    zstats <- numeric(nperm)
    tmat <- matrix(0, n, n)
    rarray <- rep(0, n)
    ncor <- length(xmat)
    w1 <- sum(xmat)/ncor
    w2 <- sum(xmat^2)
    w2 <- sqrt(w2/ncor - w1^2)
    xmat <- (xmat - w1)/w2
    w1 <- sum(ymat)/ncor
    w2 <- sum(ymat^2)
    w2 <- sqrt(w2/ncor - w1^2)
    ymat <- (ymat - w1)/w2
    
    if (correlation %in% c("Pearson", "Spearman")){  # sum of the cross products
      cresults <- .C("permute", as.double(xmat), as.double(ymat), 
                     as.integer(n), as.integer(length(xmat)), as.integer(nperm), 
                     zstats = as.double(zstats), as.double(as.vector(tmat)), 
                     as.integer(rarray),
                     PACKAGE = "RPANDA")
    }
    
    if (correlation=="Kendall"){
      cresults <- .C("permuteKendall", as.double(xmat), as.double(ymat),
                     as.integer(n), as.integer(length(xmat)), as.integer(nperm), 
                     zstats = as.double(zstats), as.double(as.vector(tmat)), 
                     as.integer(rarray),
                     PACKAGE = "RPANDA")
    }
    
    zstats <- cresults$zstats
    pval1 <- min(c(length(which(zstats >= zstats[1]))/nperm,1))
    pval2 <- min(c(length(which(zstats <= zstats[1]))/nperm,1))
    pval3 <- length(which(abs(zstats) >= abs(zstats[1])))/nperm
  } else {
    pval1 <- 0
    pval2 <- 0
    pval3 <- 0
  }
  
  c(mantelr = mantelr, pval1 = pval1, pval2 = pval2, pval3 = pval3)
}


full <- function (v) {
  n <- (1 + sqrt(1 + 8 * length(v)))/2
  if (abs(n - round(n)) > 1e-07) 
    stop("Matrix not square.")
  n <- round(n)
  full <- matrix(0, n, n)
  full[lower.tri(full)] <- v
  full2 <- t(full)
  diag(full2) <- 0
  full + full2
}



mantel_test_marginal <- function(network, tree_A, tree_B, method="Jaccard_binary", nperm=1000, correlation="Pearson"){
  
  if (!correlation %in% c("Pearson", "Spearman")) {stop("\"correlation\" must be among 'Pearson' or 'Spearman'.")}
  if (!is.numeric(nperm)) {stop("Please provide a numeric number of permutations (\"nperm\").")}
  
  compute_eco_dist <- function(network){
    # binary Jaccard distances
    if (method=="Jaccard_binary"){
      jaccard_A <- as.matrix(vegan::vegdist(t(network), "jaccard", binary=T))
      eco_A <- jaccard_A
    }
    
    # quantitative Jaccard distances
    if (method=="Jaccard_weighted"){
      jaccard_A <- as.matrix(vegan::vegdist(t(network), "jaccard", binary=F))
      eco_A <- jaccard_A
    }
    
    
    # Bray-Curtis dissimilarity 
    if (method=="Bray-Curtis"){
      bray_A <- as.matrix(vegan::vegdist(t(network), "bray", binary=F))
      eco_A <- bray_A
    }
    
    # Unifrac (generalized UniFrac, with alpha=0.5)
    if (method=="GUniFrac"){
      unifrac_A <- GUniFrac::GUniFrac(t(network), tree = tree_B, alpha=c(0.5))
      index=1
      eco_A <- unifrac_A$unifracs[,,index]
    }
    
    # Unifrac (unweighted UniFrac)
    if (method=="UniFrac_unweighted"){
      unifrac_A <- GUniFrac::GUniFrac(t(network), tree = tree_B, alpha=c(0.5))
      index=2
      eco_A <- unifrac_A$unifracs[,,index]
    }
    return(eco_A)
  }
  
  eco_A <- compute_eco_dist(network)
  
  # Perform Mantel test:
  
  # cophenetic distances
  cophe_A <- cophenetic.phylo(tree_A)

  nb_A <- ncol(network)
  nb_B <- nrow(network)
  
  results <- c(as.integer(nb_A), as.integer(nb_B), NA, NA, NA, NA, NA, NA)
  names(results) <-  c("nb_A","nb_B","mantel_cor_A","pvalue_upper_A","pvalue_lower_A", "mantel_cor_B", "pvalue_upper_B", "pvalue_lower_B")
  
  if (length(unique(as.vector(cophe_A)))<3) {
    print("The phylogenetic distance matrix is composed of only 2 different values (because of polytomies?).")
    return(results)}
  if (length(unique(as.vector(eco_A)))<3) {
    print("The ecological distance matrix is composed of only 1 value (identical patterns of interactions across species?).")
    return(results)}
  
  
  if (correlation=="Pearson") {correlation="pearson"}
  if (correlation=="Spearman") {correlation="spearman"}
  
  original_correlation <- cor(as.vector(as.dist(eco_A)), as.vector(as.dist(cophe_A)), use = "everything", method = correlation)
  
  
  # Make randomizations

  random_correlation <- vector(mode = "numeric", length = nperm)
  vector_cophe_A <- as.vector(as.dist(cophe_A))
  
  for (i in 1:nperm){
    rand_network <- network
    
    for (k in 1:nb_A){
      rand_network[,k] <- sample(rand_network[,k])
    }
    
    eco_A_rand <- compute_eco_dist(rand_network)
    
    random_correlation[i] <- cor(as.vector(as.dist(eco_A_rand)), vector_cophe_A, use = "everything", method = correlation)
  
  }
  
  results <- c(original_correlation, min(c(length(which(c(random_correlation,original_correlation)>=original_correlation))/nperm,1)), min(c(length(which(c(random_correlation,original_correlation)<=original_correlation))/nperm, 1)))

  return(results)
}



phylosignal_network <- function(network, tree_A, tree_B=NULL, method = "Jaccard_weighted", nperm = 10000, correlation = "Pearson", only_A = FALSE, permutation ="shuffle"){
  
  if (is.null(tree_B)) {only_A <- TRUE} 
  
  if (!inherits(tree_A, "phylo")) {stop("object \"tree_A\" is not of class \"phylo\".")}
  if (!is.null(tree_B)) {if (!inherits(tree_B, "phylo")) {stop("object \"tree_B\" is not of class \"phylo\".")}}
  
  if (is.null(method)) {stop("Please provide a \"method\" to compute phylogenetic signals among 'Jaccard_weighted', 'Jaccard_binary', 'Bray-Curtis', 'GUniFrac', 'UniFrac_unweighted', 'PBLM', 'PBLM_binary', and 'degree'.")}
  if (method %in% c("GUniFrac", "UniFrac_unweighted", "PBLM", "PBLM_binary")) {if (is.null(tree_B)) stop("Please provide a phylogenetic tree \"tree_B\" for guild B.")}
  if (!method %in% c("Jaccard_weighted","Jaccard_binary", "Bray-Curtis", "GUniFrac", "UniFrac_unweighted", "PBLM", "PBLM_binary", "degree")) {stop("Please provide a \"method\" to compute phylogenetic signals among 'Jaccard_weighted', 'Jaccard_binary', 'Bray-Curtis', 'GUniFrac', 'UniFrac_unweighted', 'PBLM', 'PBLM_binary', and 'degree'.")}
  
  
  if (!permutation %in% c("shuffle","nbpartners")) {stop("Please provide a type of \"permutation\" among 'shuffle' and 'nbpartners'.")}
  if (permutation!="shuffle") {if (method %in% c("PBLM", "PBLM_binary", "degree")) stop("The argument \"permutation\" is not used for this method.")}
  
  
  if (!correlation %in% c("Pearson", "Spearman", "Kendall")) {stop("Please pick a \"correlation\" among Pearson, Spearman, and Kendall.")}
  
  if (nrow(network)<2){stop("Please provide a \"network\" with at least 2 species in clade B.")}
  if (ncol(network)<2){stop("Please provide a \"network\" with at least 2 species in clade A.")}
  
  
  # Only keep species with at least 1 interaction
  network <- network[rowSums(network)>0,]
  network <- network[,colSums(network)>0]
  
  # A in columns and B in rows
  nb_A <- ncol(network)
  nb_B <- nrow(network)
  names(nb_A) <- "nb_A"
  names(nb_B) <- "nb_B"
  
  # Check names
  if (all(is.null(colnames(network)))|all(is.null(rownames(network)))) {stop("Please provide a \"network\" with row names and columns names matching the species names.")}
  
  if (!all(colnames(network) %in% tree_A$tip.label)){stop("Please provide a \"tree_A\" for all the species in clade A (the columns of the intercation network).")}
  if (only_A==FALSE) { if (!all(rownames(network) %in% tree_B$tip.label)){stop("Please provide a \"tree_B\" for all the species in clade B (the rows of the intercation network).")}}
  
  tree_A <- drop.tip(tree_A,tip=tree_A$tip.label[which(!tree_A$tip.label %in% colnames(network))])
  if (only_A==FALSE) { tree_B <- drop.tip(tree_B,tip=tree_B$tip.label[which(!tree_B$tip.label %in% rownames(network))])}
  
  
  if (!is.rooted(tree_A)){tree_A <- midpoint.root(tree_A) }
  if (only_A==FALSE) { if (!is.rooted(tree_B)){tree_A <- midpoint.root(tree_B) }}
  
  if (only_A==TRUE) { 
    network <- network[1:nrow(network),tree_A$tip.label]
  } else {
    network <- network[tree_B$tip.label,tree_A$tip.label]
  }
  
  # Mantel tests
  if (!method %in% c("PBLM_binary","PBLM")){
    
    if (permutation=="shuffle"){
      
      # binary Jaccard distances
      if (method=="Jaccard_binary"){
        jaccard_A <- as.matrix(vegan::vegdist(t(network), "jaccard", binary=T))
        if (only_A==FALSE) jaccard_B <- as.matrix(vegan::vegdist(network, "jaccard", binary=T))
        eco_A <- jaccard_A
        if (only_A==FALSE) eco_B <- jaccard_B
      }
      
      # quantitative Jaccard distances
      if (method=="Jaccard_weighted"){
        jaccard_A <- as.matrix(vegan::vegdist(t(network), "jaccard", binary=F))
        if (only_A==FALSE) jaccard_B <- as.matrix(vegan::vegdist(network, "jaccard", binary=F))
        eco_A <- jaccard_A
        if (only_A==FALSE) eco_B <- jaccard_B
      }
      
      
      # Bray-Curtis dissimilarity 
      if (method=="Bray-Curtis"){
        bray_A <- as.matrix(vegan::vegdist(t(network), "bray", binary=F))
        if (only_A==FALSE) bray_B <- as.matrix(vegan::vegdist(network, "bray", binary=F))
        eco_A <- bray_A
        if (only_A==FALSE) eco_B <- bray_B
      }
      
      # Unifrac (generalized UniFrac, with alpha=0.5)
      if (method=="GUniFrac"){
        if (only_A==FALSE) unifrac_B <- GUniFrac::GUniFrac(network, tree = tree_A, alpha=c(0.5))
        unifrac_A <- GUniFrac::GUniFrac(t(network), tree = tree_B, alpha=c(0.5))
        index=1
        eco_A <- unifrac_A$unifracs[,,index]
        if (only_A==FALSE) eco_B <- unifrac_B$unifracs[,,index]
      }
      
      # Unifrac (unweighted UniFrac)
      if (method=="UniFrac_unweighted"){
        if (only_A==FALSE) unifrac_B <- GUniFrac::GUniFrac(network, tree = tree_A, alpha=c(0.5))
        unifrac_A <- GUniFrac::GUniFrac(t(network), tree = tree_B, alpha=c(0.5))
        index=2
        eco_A <- unifrac_A$unifracs[,,index]
        if (only_A==FALSE) eco_B <- unifrac_B$unifracs[,,index]
      }
      
      # Degree
      if (method=="degree"){
        network_binary <- network
        network_binary[network_binary>0] <- 1
        
        eco_A <- as.matrix(dist(colSums(network_binary)))
        if (only_A==FALSE) eco_B <- as.matrix(dist(rowSums(network_binary)))
      }
      
      
      # Perform Mantel test:
      
      # cophenetic distances
      cophe_A <- cophenetic.phylo(tree_A)
      if (only_A==FALSE) cophe_B <- cophenetic.phylo(tree_B)
      
      results <- c(as.integer(nb_A), as.integer(nb_B), NA, NA, NA, NA, NA, NA)
      names(results) <-  c("nb_A","nb_B","mantel_cor_A","pvalue_upper_A","pvalue_lower_A", "mantel_cor_B", "pvalue_upper_B", "pvalue_lower_B")
      
      if (length(unique(as.vector(cophe_A)))<3) {
        print("The phylogenetic distance matrix of guild A is composed of only 2 different values (because of polytomies?).")
        return(results)}
      if (only_A==FALSE) {if (length(unique(as.vector(cophe_B)))<3) {
        print("The phylogenetic distance matrix of guild B only is composed of only 2 different values (because of polytomies?).")
        return(results)}}
      if (length(unique(as.vector(eco_A)))<3) {
        print("The ecological distance matrix of guild A is composed of only 1 value (identical patterns of interactions across species?).")
        return(results)}
      if (only_A==FALSE) {if (length(unique(as.vector(eco_B)))<3) {
        print("The ecological distance matrix of guild B is composed of only 1 value (identical patterns of interactions across species?).")
        return(results)}}
      
      
      if (correlation=="Pearson"){
        mantel_A <- RPANDA::mantel_test(as.dist(eco_A) ~ as.dist(cophe_A),  nperm = nperm, correlation="Pearson")
        if (only_A==FALSE) mantel_B <- RPANDA::mantel_test(as.dist(eco_B) ~ as.dist(cophe_B),  nperm = nperm, correlation="Pearson")
      }
      
      if (correlation=="Spearman"){
        mantel_A <- RPANDA::mantel_test(as.dist(eco_A) ~ as.dist(cophe_A),  nperm = nperm, correlation="Spearman")
        if (only_A==FALSE) mantel_B <- RPANDA::mantel_test(as.dist(eco_B) ~ as.dist(cophe_B),  nperm = nperm, correlation="Spearman")
      }
      
      if (correlation=="Kendall"){
        mantel_A <- RPANDA::mantel_test(as.dist(eco_A) ~ as.dist(cophe_A),  nperm = nperm, correlation="Kendall")
        if (only_A==FALSE) mantel_B <- RPANDA::mantel_test(as.dist(eco_B) ~ as.dist(cophe_B),  nperm = nperm, correlation="Kendall")
      }
      
      if (only_A==TRUE) mantel_B <- c(NA, NA, NA)
      results <- c(as.integer(nb_A), as.integer(nb_B), mantel_A[1], mantel_A[2], mantel_A[3], mantel_B[1], mantel_B[2], mantel_B[3])
      names(results) <-  c("nb_A","nb_B","mantel_cor_A","pvalue_upper_A","pvalue_lower_A", "mantel_cor_B", "pvalue_upper_B", "pvalue_lower_B")
      
    }
    
    
    if (permutation=="nbpartners"){
      mantel_A <- #RPANDA::
        mantel_test_marginal(network, tree_A, tree_B, method, nperm, correlation)
      if (only_A==FALSE) {mantel_B <- #RPANDA::
        mantel_test_marginal(t(network), tree_B, tree_A, method, nperm, correlation)
      }else{mantel_B <- c(NA, NA, NA)}
    }
    
    results <- c(as.integer(nb_A), as.integer(nb_B), mantel_A[1], mantel_A[2], mantel_A[3], mantel_B[1], mantel_B[2], mantel_B[3])
    names(results) <-  c("nb_A","nb_B","mantel_cor_A","pvalue_upper_A","pvalue_lower_A", "mantel_cor_B", "pvalue_upper_B", "pvalue_lower_B")
    
    return(results)
  }
  
  # PBLM (non binary)
  if ((method=="PBLM")&(only_A==FALSE)){
    model_pblm <- R.utils::withTimeout(pblm(assocs=network, tree1=tree_B, tree2=tree_A, bootstrap=F, nreps=0), timeout = 60*60*24, onTimeout = "silent")
    
    if (!is.null(model_pblm)) {
      results <- c(as.integer(nb_A), as.integer(nb_B), model_pblm$signal.strength$estimate[2], model_pblm$signal.strength$estimate[1], model_pblm$MSE )
      names(results) <- c("nb_A", "nb_B", "dA", "dB", "MSETotal", "MSEFull", "MSEStar", "MSEBase")
      return(unlist(results))
    }else{
      results <- c(as.integer(nb_A), as.integer(nb_B), NA, NA, NA, NA, NA, NA)
      names(results) <- c("nb_A", "nb_B", "dA", "dB", "MSETotal", "MSEFull", "MSEStar", "MSEBase")
      return(unlist(results))}
  }
  
  # PBLM binary
  if ((method=="PBLM_binary")&(only_A==FALSE)){
    network_binary <- network
    network_binary[network_binary>0] <- 1
    
    model_pblm <- R.utils::withTimeout(pblm(assocs=network_binary, tree1=tree_B, tree2=tree_A, bootstrap=F, nreps=0), timeout = 60*60*24, onTimeout = "silent")
    
    if (!is.null(model_pblm)) {
      results <- c(as.integer(nb_A), as.integer(nb_B), model_pblm$signal.strength$estimate[2], model_pblm$signal.strength$estimate[1], model_pblm$MSE )
      names(results) <- c("nb_A", "nb_B", "dA", "dB", "MSETotal", "MSEFull", "MSEStar", "MSEBase")
      return(unlist(results))
    }else{
      results <- c(as.integer(nb_A), as.integer(nb_B), NA, NA, NA, NA, NA, NA)
      names(results) <- c("nb_A", "nb_B", "dA", "dB", "MSETotal", "MSEFull", "MSEStar", "MSEBase")
      return(unlist(results))}
  }
}



phylosignal_sub_network <- function(network, tree_A, tree_B=NULL, method = "Jaccard_weighted", 
                                    nperm = 1000, correlation = "Pearson", minimum=10, degree=FALSE, permutation ="shuffle"){
  
  host_tree <- tree_A
  symbiont_tree <- tree_B
  
  if (!inherits(host_tree, "phylo")) {stop("object \"tree_A\" is not of class \"phylo\".")}
  if (!is.null(tree_B)) {if (!inherits(symbiont_tree, "phylo")) {stop("object \"tree_B\" is not of class \"phylo\".")}}
  if (!method %in% c("Jaccard_weighted","Jaccard_binary", "Bray-Curtis", "GUniFrac", "UniFrac_unweighted")) {stop("Please provide a \"method\" to compute phylogenetic signals.")}
  
  if (all(is.null(colnames(network)))|all(is.null(rownames(network)))) {stop("Please provide a network with row names and columns names matching the species names.")}
  
  if (method %in% c("GUniFrac", "UniFrac_unweighted")) {if (is.null(tree_B)) stop("Please provide a phylogenetic tree \"tree_B\" for guild B.")}
  
  if (!correlation %in% c("Pearson", "Spearman", "Kendall")) {stop("Please pick a \"correlation\" among Pearson, Spearman, and Kendall.")}
  
  if (nrow(network)<2){stop("Please provide a \"network\" with at least 2 species in clade B.")}
  if (ncol(network)<2){stop("Please provide a \"network\" with at least 2 species in clade A.")}
  
  if (minimum<2){stop("The minimal number of descending species (\"minimum\") must be with at least of 2 (or even larger!).")}
  
  if (!permutation %in% c("shuffle","nbpartners")) {stop("Please provide a type of \"permutation\" among 'shuffle' and 'nbpartners'.")}

  # only keep species having at least one interaction
  network <- network[rowSums(network)>0,]
  network <- network[,colSums(network)>0]
  host_tree <- drop.tip(host_tree, tip=host_tree$tip.label[!host_tree$tip.label %in% colnames(network)])
  if (!is.null(symbiont_tree)){
    symbiont_tree <- drop.tip(symbiont_tree, tip=symbiont_tree$tip.label[!symbiont_tree$tip.label %in% rownames(network)])
    network <- network[symbiont_tree$tip.label,host_tree$tip.label]
  }else{
    network <- network[,host_tree$tip.label]
  }
  
  # check a second time (in case of a missing species) 
  network <- network[rowSums(network)>0,]
  network <- network[,colSums(network)>0]
  host_tree <- drop.tip(host_tree, tip=host_tree$tip.label[!host_tree$tip.label %in% colnames(network)])
  if (!is.null(symbiont_tree)){
    symbiont_tree <- drop.tip(symbiont_tree, tip=symbiont_tree$tip.label[!symbiont_tree$tip.label %in% rownames(network)])
    network <- network[symbiont_tree$tip.label,host_tree$tip.label]
  }else{
    network <- network[,host_tree$tip.label]
  }
  
  
  set.seed(1)
  
  nb_sub_clades <- 0
  results_sub_clades <- c()
  for (i in sort(unique(host_tree$edge[,1]))){  # include root  and can be non binary
    sub_host_tree <- extract.clade(host_tree, i)
    if (Ntip(sub_host_tree)>=minimum){
      sub_network <- network[,sub_host_tree$tip.label]
      sub_network <- sub_network[which(rowSums(sub_network)>0),,drop=F]
      if (!is.null(symbiont_tree)){
        sub_symbiont_tree <- drop.tip(symbiont_tree, tip= symbiont_tree$tip.label[!symbiont_tree$tip.label %in% rownames(sub_network)])
      }else{sub_symbiont_tree <- NULL}
      
      if (nrow(sub_network)>1){
        nb_sub_clades <- nb_sub_clades+1
        mantel_test <- phylosignal_network(sub_network, sub_host_tree, sub_symbiont_tree, method = method, nperm = nperm, correlation = correlation, permutation = permutation)
        
        if (degree==TRUE){
          mantel_degree <- rep("NA", 5)
          tryCatch({
            mantel_degree <- phylosignal_network(sub_network, sub_host_tree, sub_symbiont_tree, method = "degree", nperm = nperm, correlation = correlation)
          }, error=function(e){cat("clade ",i,": ", conditionMessage(e), "\n")})
          results_sub_clades <- rbind(results_sub_clades, c(i, mantel_test[1:5],NA,NA, mantel_degree[3:5] ))
        }else{
          results_sub_clades <- rbind(results_sub_clades, c(i, mantel_test[1:5],NA,NA))
        }
        
      }
    }
  }
  if (degree==TRUE){
    colnames(results_sub_clades) <- c("node", "nb_A", "nb_B", "mantel_cor", "pvalue_upper", "pvalue_lower", "pvalue_upper_corrected","pvalue_lower_corrected", "degree_mantel_cor", "degree_pvalue_upper", "degree_pvalue_lower") 
  }else{
    colnames(results_sub_clades) <- c("node", "nb_A", "nb_B", "mantel_cor", "pvalue_upper", "pvalue_lower", "pvalue_upper_corrected","pvalue_lower_corrected") 
  }
  results_sub_clades <- data.frame(results_sub_clades, stringsAsFactors = F)
  results_sub_clades$nb_A <- as.integer(as.numeric(results_sub_clades$nb_A))
  results_sub_clades$nb_B <- as.integer(as.numeric(results_sub_clades$nb_B))
  results_sub_clades$mantel_cor <- as.numeric(results_sub_clades$mantel_cor)
  results_sub_clades$pvalue_upper <- as.numeric(results_sub_clades$pvalue_upper)
  results_sub_clades$pvalue_lower <- as.numeric(results_sub_clades$pvalue_lower)
  
  results_sub_clades$pvalue_upper_corrected <- results_sub_clades$pvalue_upper*nb_sub_clades
  results_sub_clades$pvalue_lower_corrected <- results_sub_clades$pvalue_lower*nb_sub_clades
  results_sub_clades$pvalue_upper_corrected[results_sub_clades$pvalue_upper_corrected>1] <- 1
  results_sub_clades$pvalue_lower_corrected[results_sub_clades$pvalue_lower_corrected>1] <- 1
  
  if (degree==TRUE){
    results_sub_clades$degree_mantel_cor <- as.numeric(results_sub_clades$degree_mantel_cor)
    results_sub_clades$degree_pvalue_upper <- as.numeric(results_sub_clades$degree_pvalue_upper)
    results_sub_clades$degree_pvalue_lower <- as.numeric(results_sub_clades$degree_pvalue_lower)
    
    results_sub_clades$degree_pvalue_upper_corrected <- results_sub_clades$degree_pvalue_upper*nb_sub_clades
    results_sub_clades$degree_pvalue_lower_corrected <- results_sub_clades$degree_pvalue_lower*nb_sub_clades
    results_sub_clades$degree_pvalue_upper_corrected[results_sub_clades$degree_pvalue_upper_corrected>1] <- 1
    results_sub_clades$degree_pvalue_lower_corrected[results_sub_clades$degree_pvalue_lower_corrected>1] <- 1
  }
  
  return(results_sub_clades)
}


plot_phylosignal_sub_network <- function(tree_A, results_sub_clades, network=NULL, legend=TRUE, show.tip.label=FALSE, where="bottomleft"){
  
  set.seed(1)
  host_tree <- tree_A
  
  if (!inherits(host_tree, "phylo")) {stop("object \"tree_A\" is not of class \"phylo\".")}
  
  
  if (!is.null(network)){
    network <- network[rowSums(network)>0,]
    network <- network[,colSums(network)>0]
    host_tree <- drop.tip(host_tree, tip=host_tree$tip.label[!host_tree$tip.label %in% colnames(network)])
    network <- network[,host_tree$tip.label]
  }
  
  
  if ((Ntip(host_tree)+1)!=results_sub_clades$node[1]){
    stop("object \"tree_A\" contains more node than \"results_sub_clades\". Remove these nodes from \"tree_A\" before plotting.")
  }
  
  plot(host_tree, show.tip.label=show.tip.label)
  # significant and R>=0.05
  nodes=results_sub_clades$node[intersect(which(results_sub_clades$pvalue_upper_corrected<=0.05),which(results_sub_clades$mantel_cor>=0.5))]
  nodelabels(node=nodes,pie=rep(1,length(nodes)),piecol="#78281f",cex=0.5)
  
  # significant and 0.3<R<0.5
  nodes=results_sub_clades$node[intersect(which(results_sub_clades$pvalue_upper_corrected<=0.05),intersect(which(results_sub_clades$mantel_cor<0.5), which(results_sub_clades$mantel_cor>=0.3)))]
  nodelabels(node=nodes,pie=rep(1,length(nodes)),piecol="#b03a2e",cex=0.5)
  
  # significant and 0.1<R<0.3
  nodes=results_sub_clades$node[intersect(which(results_sub_clades$pvalue_upper_corrected<=0.05),intersect(which(results_sub_clades$mantel_cor<0.3), which(results_sub_clades$mantel_cor>=0.1)))]
  nodelabels(node=nodes,pie=rep(1,length(nodes)),piecol="#ec7063",cex=0.5)
  
  # significant and 0.1>R
  nodes=results_sub_clades$node[intersect(which(results_sub_clades$pvalue_upper_corrected<=0.05),which(results_sub_clades$mantel_cor<0.1))]
  nodelabels(node=nodes,pie=rep(1,length(nodes)),piecol="#f5b7b1",cex=0.5)
  
  # not significant 
  nodes=results_sub_clades$node[which(results_sub_clades$pvalue_upper_corrected>0.05)]
  nodelabels(node=nodes,pie=rep(1,length(nodes)),piecol="#aed6f1",cex=0.5)
  
  # too small
  nodes=((Ntip(host_tree)+2):(2*Ntip(host_tree)-1))[!(Ntip(host_tree)+2):(2*Ntip(host_tree)-1) %in% results_sub_clades$node]
  nodelabels(node=nodes,pie=rep(1,length(nodes)),piecol="#5d6d7e",cex=0.2)
  
  # add legend
  if (legend==TRUE) {legend(where, legend = c("0.5 < R","0.3 < R < 0.5","0.1 < R < 0.3","R < 0.1","Not significant","Small sub-clade"), pch = 19, 
                            col = c("#78281f","#b03a2e","#ec7063","#f5b7b1","#aed6f1","#5d6d7e"), bg="transparent", bty = "n", cex=0.89, y.intersp=0.9)}
  
}

