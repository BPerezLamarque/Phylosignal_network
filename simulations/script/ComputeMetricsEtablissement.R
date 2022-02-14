

computeModulesInSilence=function(M,steps,deep=F){
  capture.output(y<-computeModules(M,steps = steps, deep = deep,forceLPA = T));return(y)
}

compute.metrics=function(output,genealogy,distanceH,distanceP,threshold,rP,rH,min_interaction_strength=0.1,module=T){
  species=define.species(genealogy,threshold,distanceH,distanceP,verbose=F)
  build_network=T
  dist_phylo=F
  res=list(species=species)
  nP=1
  nH=1
    
  if( is.null(species$Hphylo$tree)){
    if(is.null(species$Pphylo$tree)){
      build_network=F
    }
    gammaH=NA
    betaH=NA
    spectraH=NA
    asymmetryH=NA
    peackednessH=NA
    }else{
      if(!is.null(species$Pphylo$tree)){
        dist_phylo=T}
    tmp=multi2di(species$Hphylo$tree)
    gammaH=gammaStat(tmp)
    betaH=mle_beta_int(tmp,bmax=100)
    nH=ncol(species$Hphylo$mean.trait)
    spectraH=spectR(species$Hphylo$tree)
    asymmetryH=spectraH$asymmetry
    peackednessH=spectraH$peakedness1
    }
  res=c(res,list(gammaH=gammaH,betaH=betaH,asymmetryH=asymmetryH,peackednessH=peackednessH))
  
  if(! is.null(species$Pphylo$tree)){
    tmp=multi2di(species$Pphylo$tree)
    gammaP=gammaStat(tmp)
    betaP=mle_beta_int(tmp,bmax=100)
    nP=ncol(species$Pphylo$mean.trait)
    spectraP=spectR(species$Pphylo$tree)
    asymmetryP=spectraP$asymmetry
    peackednessP=spectraP$peakedness1
  }else{
    gammaP=NA
    betaP=NA
    spectraP=NA
    asymmetryP=NA
    peackednessP=NA
  }
  res=c(res,list(gammaP=gammaP,betaP=betaP,asymmetryP=asymmetryP,peackednessP=peackednessP))
  

  if(build_network){
    network=build.network(output$P,output$H,p,species,rP=rP,rH=rH,seuil=min_interaction_strength)
    M=as.matrix(network)
    rownames(M)=paste("P",1:nrow(M))
    colnames(M)=paste("H",1:ncol(M))
    if(nP>1 & nH>1){
      NODF=try(nested(M,method = "NODF2"))
      if(inherits(NODF,"try-error")) NODF=NA
      weightedNODF=try(nested(M,method = "weighted NODF"))
      if(inherits(weightedNODF,"try-error")) NODF=NA
      if(module){
      modules=try(computeModulesInSilence(M,steps = 1e7))
      if(inherits(modules,"try-error")){
        modularity=NA
        nModules=NA
      }else if (is.null(modules)){
        modularity=NA
        nModules=NA
      }else{
        modularity=slot(modules,"likelihood")
        nModules=sum(slot(modules,"modules")[,1]==1)
      }
      res=c(res,list(modularity=modularity,nModules=nModules))}
      res=c(res,list(NODF=NODF,weightedNODF=weightedNODF))
      }

    res=c(res,list(network=network))
  }else{
    modularity=NA
    nModules=NA
    NODF=NA
    weightedNODF=NA
        res=c(res,list(modularity=modularity,nModules=nModules,NODF=NODF,weightedNODF=weightedNODF))

    res=c(res,list(network=1))
  }
  
  if(dist_phylo){
    distance_PH=JSDtree(list(species$Pphylo$tree,species$Hphylo$tree))[1,2]
    res=c(res,list(distance_PH=distance_PH))
  }else{
    res=c(res,list(distance_PH=NA))
  }
  
  res=c(res,list(nH=nH,nP=nP))
  return(res)
  
}

add.iteration=function(thresholds,output,genealogy,nx,ny,NG,dSpace,traitDimention,muH,muP,nH,verbose=2*NG,
                       nP,alphaP,alphaH,rP,rH,effect,totalTime,iniP = 0,timeStep = 20,min_interaction_strength=0.1,module=T){
  if(is.null(output)){
    output=model.spatial(nx=nx,ny=ny,NG=NG,dSpace=dSpace,D=traitDimention,muH=muH,muP=muP,nH=nH,
                         nP=nP,alphaP =alphaP,alphaH =alphaH,thin=1,rP=rP,rH=rH,
                         effect = effect,verbose = verbose,iniP = iniP,
                         oneByOne = F,timeStep = 20)
    genealogy=make.gen(output,verbose = F)
  }else{
    output=model.spatial(nx=nx,ny=ny,NG=NG,dSpace=dSpace,D=traitDimention,muH=muH,muP=muP,nH=nH,
                         nP=nP,alphaP =alphaP,alphaH =alphaH,thin=1,rP=rP,rH=rH,
                         effect = effect,verbose = 2*NG,
                         oneByOne = F,timeStep = 20,P = output$P,H=output$H)
    genealogy=make.gen(output,treeP = genealogy$P, treeH = genealogy$H,verbose = F)
  }
  
  newline=list()
  
  xP=output$P
  xH=output$H
  rangeP=sapply(1:traitDimention,function(i){range(output$P[i,])})
  rangeH=sapply(1:traitDimention,function(i){range(output$H[i,])})
  tmp=multi2di(genealogy$P)
  gammaP=gammaStat(tmp)
  betaP=maxlik.betasplit(tmp)$max_lik
  tmp=multi2di(genealogy$H)
  gammaH=gammaStat(tmp)
  betaH=maxlik.betasplit(tmp)$max_lik
  corHP=sapply(1:traitDimention,function(j){y=try(cor(genealogy$P$x.tip[j,],sapply(genealogy$P$tip.label,function(i){genealogy$H$x.tip[j,genealogy$H$tip.label==i]})))
  if(inherits(y,"try-error")){NA}else{y}})
  distanceP=compute.dist(genealogy$P,verbose = F)
  distanceH=compute.dist(genealogy$H,verbose = F)
  newline=list(xP=xP,xH=xH,betaP=betaP,betaH=betaH,rangeP=rangeP,rangeH=rangeH,gammaP=gammaP,gammaH=gammaH,cor=corHP,genealogy=genealogy,totalTime=totalTime+NG,
               rootP=totalTime+NG-genealogy$P$ini,rootH=totalTime+NG-genealogy$H$ini)
  names=names(newline)

  for(threshold in thresholds){
    cat(paste(";", threshold))
    newline=c(newline,list(compute.metrics(output,genealogy,distanceH,distanceP,threshold,rP,rH,min_interaction_strength,module=module)))
  }
  names(newline)=c(names,as.character(thresholds))
  
  return(list(newline=newline,output=output))
}

plot.time.series=function(listTime,thresholds,name,trait.id=1){
  if("0" %in% thresholds){
    if(name=="cor"){
      Y=c()
      times=c()
      traitDimention=length(listTime[[1]]$cor)
      for (i in 1:length(listTime)){
        Y=cbind(Y,as.matrix(listTime[[i]][["cor"]]))
        times=c(times,listTime[[i]][["totalTime"]])
      }
      if(traitDimention>1){
        Colors = colorRampPalette(c("steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4"))(traitDimention) 
        plot(times,Y[1,],type='l',ylim=c(-1,1),col=Colors[1])
        for(i in 2:traitDimention){
          lines(times,Y[i,],col=Colors[i])
        }}else{
          plot(times,Y)
        }
      lines(c(times[1],times[length(times)]),c(0,0),col="red")
    }else if(name=="range"){
      YP=c()
      YH=c()
      times=c()
      traitDimention=length(listTime[[1]]$cor)
      for (i in 1:length(listTime)){
        YP=cbind(YP,as.matrix(as.vector(listTime[[i]][["rangeP"]])))
        YH=cbind(YH,as.matrix(as.vector(listTime[[i]][["rangeH"]])))
        times=c(times,listTime[[i]][["totalTime"]])
      }
      if(traitDimention>1){
        Colors = colorRampPalette(c("steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4"))(traitDimention) 
        plot(times,YP[1,],type='l',ylim=range(cbind(YP,YH)),col=Colors[1])
        lines(times,YH[1,],lty=4,col=Colors[1])
        lines(times,YH[2,],lty=4,col=Colors[1])
        lines(times,YP[2,],lty=1,col=Colors[1])
        for(i in 2:traitDimention){
          lines(times,YP[2*i-1,],col=Colors[i])
          lines(times,YH[2*i-1,],lty=4,col=Colors[i])
          lines(times,YP[2*i,],col=Colors[i])
          lines(times,YH[2*i,],lty=4,col=Colors[i])
        }}else{
          plot(times,YP,ylim=range(cbind(YP,YH)))
          lines(time,YH,lty=4)
        }
    }else if(name %in% c("totalTime")){
      Y=c()
      times=c()
      for (i in 1:length(listTime)){
        y=listTime[[i]][[name]]
        if(is.null(y)){ yH=NA}
        Y=c(Y,y)
        times=c(times,listTime[[i]][["totalTime"]])
      }
      max=max(Y,na.rm = T)
      min=min(Y,na.rm = T)
      plot(times,Y,col="chocolate2",type = 'l',ylim=c(min,max))
    }else{
      YH=c()
      YP=c()
      times=c()
      for (i in 1:length(listTime)){
        yH=listTime[[i]][[paste(name,"H",sep="")]]
        if(is.null(yH)){ yH=NA}
        YH=c(YH,yH)
        yP=listTime[[i]][[paste(name,"P",sep="")]]
        if(is.null(yP)){ yP=NA}
        YP=c(YP,yP)
        times=c(times,listTime[[i]][["totalTime"]])
      }
      max=max(max(YH,na.rm = T),max(YP,na.rm = T),na.rm = T)
      min=min(min(YH,na.rm = T),min(YP,na.rm = T),na.rm = T)
      plot(times,YH,col="chocolate2",type = 'l',ylim=c(min,max))
      lines(times,YP,col="olivedrab3")
  }
  }else{
    newPlot=T
    Green=colorRampPalette(c("olivedrab1","olivedrab4"))(length(thresholds)) 
    Orange=colorRampPalette(c("chocolate1","chocolate4"))(length(thresholds)) 
    k=1
    times=c()
    if(name %in% c("modularity", "nModules","distance_PH","NODF","weightedNODF")){
      Y=list()
      for(threshold in thresholds){
        Y=c(Y,list(c()))
      for (i in 1:length(listTime)){
        y=listTime[[i]][[threshold]][[name]]
        if(is.null(y)){ y=NA}
        Y[[k]]=c(Y[[k]],y)
        if (k==1) times=c(times,listTime[[i]][["totalTime"]])
      }
        k=k+1}
      max=max(sapply(Y,function(x){max(x,na.rm=T)}),na.rm = T)
      min=min(sapply(Y,function(x){min(x,na.rm=T)}),na.rm = T)
      plot(times,Y[[1]],col=Orange[1],type = 'l',ylim=c(min,max))
      if(length(thresholds)>1){
        for(k in 2:length(thresholds)){
          lines(times,Y[[k]],col=Orange[k])
        }
      }
    }else{
      YH=list()
      YP=list()
      times=c()
      for(threshold in thresholds){
        YP=c(YP,list(c()))
        YH=c(YH,list(c()))
        for (i in 1:length(listTime)){
          yH=listTime[[i]][[threshold]][[paste(name,"H",sep="")]]
          if(is.null(yH)){ yH=NA}
          YH[[k]]=c(YH[[k]],yH)
          yP=listTime[[i]][[threshold]][[paste(name,"P",sep="")]]
          if(is.null(yP)){ yP=NA}
          YP[[k]]=c(YP[[k]],yP)
          if(k==1) times=c(times,listTime[[i]][["totalTime"]])
        }
        k=k+1}
      max=max(c(sapply(YP,function(x){max(x,na.rm=T)}),sapply(YH,function(x){max(x,na.rm=T)})),na.rm = T)
      min=min(c(sapply(YP,function(x){min(x,na.rm=T)}),sapply(YH,function(x){min(x,na.rm=T)})),na.rm = T)
      plot(times,YH[[1]],col=Orange[1],type = 'l',ylim=c(min,max))
      lines(times,YP[[1]],col=Green[1])
      if(length(thresholds)>1){
        for(k in 2:length(thresholds)){
          lines(times,YH[[k]],col=Orange[k])
          lines(times,YP[[k]],col=Green[k])
        }
      }
    }}
}