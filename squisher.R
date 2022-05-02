### Take a matrix and conversion list and squish them (old function)
# load("D:/Dropbox/projects/rlab/data/target_nbl_bcca_ensg_counts.rda")
# ensgmat<-target_nbl_bcca_ensg_counts
# ensgmat<-as.matrix(ensgmat)
# source("../shared/functions/geneids.R")
# tmp<-ens2eg(rownames(ensgmat))
# convlist<-eg2sym(tmp)
# names(convlist)<-names(tmp)
# rawcounts<-oldsquish(ensgmat,convlist=convlist,method="sum",verbose=TRUE)
oldsquish<-function(inexp,convlist,verbose=FALSE,method="sum"){
  inexp<-inexp[!is.na(convlist),]
  egs<-convlist[!is.na(convlist)]
  
  outexp<-matrix(0,nrow=length(unique(egs)),ncol=ncol(inexp))
  rownames(outexp)<-unique(egs)
  colnames(outexp)<-colnames(inexp)
  
  
  # Sum values of Rows mapping the same Gene
  if(method=="sum"){
    for(i in 1:nrow(outexp)){
      eghere<-rownames(outexp)[i]
      rowhere<-which(egs==eghere)
      if(verbose){
        message(i,"/",nrow(outexp),": ",length(rowhere))
      }
      if(length(rowhere)==1){
        outexp[eghere,]<-inexp[rowhere,]
      } else {
        mysum<-apply(inexp[rowhere,],2,sum,na.rm=TRUE)
        outexp[eghere,]<-mysum
      }
    }
  }
  
  # Average values of Rows mapping the same Gene
  if(method=="average"){
    for(i in 1:nrow(outexp)){
      eghere<-rownames(outexp)[i]
      rowhere<-names(which(egs==eghere))
      if(verbose){
        message(i,"/",nrow(outexp),": ",length(rowhere))
      }
      if(length(rowhere)==1){
        outexp[eghere,]<-inexp[rowhere,]
      } else {
        mysum<-apply(inexp[rowhere,],2,mean,na.rm=TRUE)
        outexp[eghere,]<-mysum
      }
    }
  }
  
  # Average values of Rows mapping the same Gene using median
  if(method=="median"){
    for(i in 1:nrow(outexp)){
      eghere<-rownames(outexp)[i]
      rowhere<-names(which(egs==eghere))
      if(verbose){
        message(i,"/",nrow(outexp),": ",length(rowhere))
      }
      if(length(rowhere)==1){
        outexp[eghere,]<-inexp[rowhere,]
      } else {
        mysum<-apply(inexp[rowhere,],2,median,na.rm=TRUE)
        outexp[eghere,]<-mysum
      }
    }
  }
  
  # Keep only max value of Rows mapping the same Gene (most expressed probeset, taken alternatively from every sample)
  if(method=="highest"){
    for(i in 1:nrow(outexp)){
      eghere<-rownames(outexp)[i]
      rowhere<-names(which(egs==eghere))
      if(verbose){
        message(i,"/",nrow(outexp),": ",length(rowhere))
      }
      if(length(rowhere)==1){
        outexp[eghere,]<-inexp[rowhere,]
      } else {
        mysum<-apply(inexp[rowhere,],2,max,na.rm=TRUE)
        outexp[eghere,]<-mysum
      }
    }
  }
  
  #TODO: keep the highest variance, keep the most expressed, etc.
  
  
  
  return(outexp)
}






### Take a matrix and conversion list and squish them (fast function)
# load("D:/Dropbox/projects/rlab/data/target_nbl_bcca_ensg_counts.rda")
# ensgmat<-target_nbl_bcca_ensg_counts
# ensgmat<-as.matrix(ensgmat)
# source("../shared/functions/geneids.R")
# tmp<-ens2eg(rownames(ensgmat))
# convlist<-eg2sym(tmp)
# names(convlist)<-names(tmp)
# rawcounts<-oldsquish(ensgmat,convlist=convlist,method="sum",verbose=TRUE)
library(dplyr)
squish<-function(inexp,convlist,method="sum"){
  initial_colnames<-colnames(inexp)
  obj<-data.frame(inexp)
  obj<-cbind(convlist[rownames(obj)],obj)
  colnames(obj)[1]<-"geneid"
  library(dplyr)
  if(method=="sum"){
    obj<-obj%>%group_by(geneid)%>%summarise(across(everything(),sum))
  }else if(method=="average"){
    obj<-obj%>%group_by(geneid)%>%summarise(across(everything(),mean))
  }else if(method=="median"){
    obj<-obj%>%group_by(geneid)%>%summarise(across(everything(),median))
  } else{
    message("Method ",method," not supported")
    break()
  }
  newcounts<-as.matrix(obj[,2:ncol(obj)])
  rownames(newcounts)<-pull(obj,geneid)
  newcounts<-newcounts[!is.na(rownames(newcounts)),]
  colnames(newcounts)<-initial_colnames
  return(newcounts)
}




# # Benchmark, old squish vs new squish
# ensgmat<-matrix(nrow=20000,ncol=500)
# ensgmat[]<-sample(0:10000,nrow(ensgmat)*ncol(ensgmat),replace=TRUE) # Integer (e.g. counts) example
# #ensgmat[]<-rnorm(nrow(ensgmat)*ncol(ensgmat)) # Numeric (e.g. expmat) example
# colnames(ensgmat)<-paste0("sample",1:ncol(ensgmat))
# rownames(ensgmat)<-paste0("ENSG",1:nrow(ensgmat))
# convlist<-sample(c(NA,paste0("gene",1:10000)),nrow(ensgmat),replace=TRUE)
# names(convlist)<-rownames(ensgmat)
# system.time(oldcounts<-oldsquish(ensgmat,convlist=convlist,method="sum",verbose=TRUE)) # 243.97 seconds (20000 esng, 500 samples, 10000 genes)
# oldcounts<-oldcounts[order(rownames(oldcounts)),]
# system.time(newcounts<-squish(ensgmat,convlist=convlist,method="sum")) # 9.01 seconds (20000 esng, 500 samples, 10000 genes)
# oldcounts-newcounts
# all.equal(oldcounts,newcounts)
# identical(oldcounts,newcounts) # FALSE only if input is integer







