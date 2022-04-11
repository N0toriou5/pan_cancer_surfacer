# SURFACER core pipeline to run on a server
# The PAN-Cancer Surfaceome Landscape . generate all MRAs. This code summarize steps
setwd("set_your_home_dir/")
# load all required libraries and set up the environment
library(stringr)
library(matrixStats)
library(corto)
library(ggforce)
library(edgeR)
library(DESeq2)
load("data/surfacer_2022.rda")
tissues<-c("adrenal gland","breast","brain","esophagus","liver","lung","ovary","pancreas", "prostate", "skin", 
           "stomach", "testis", "thyroid", "uterus") # NB: 'colon tissue is a valid tissue for GTEx only, while colorectal is the key word for TCGA'

### Try the TCGABiolinks mode----
library(TCGAbiolinks)
library(SummarizedExperiment)
library(recount)

###conversion of uuids to TCGA barcodes
library(TCGAutils)
library(limma)
library(biomaRt)
library(factoextra)
source("data/squisher.R")
source("data/ensembl2symbol.R")


for (tissue in tissues) {
  ########Query from Recount2 platform#######
  recount.gtex<-TCGAquery_recount2(project="GTEX", tissue=tissue)
  #to get the SE object
  name<-paste0("GTEX_",tissue)
  name<-sub(" ","_",name)
  SE.recount.gtex <- recount.gtex[[name]]
  # same for TCGA
  recount.tcga<-TCGAquery_recount2(project="TCGA", tissue=tissue)
  name<-paste0("TCGA_",tissue)
  name<-sub(" ","_",name)
  SE.recount.tcga <- recount.tcga[[name]]
  # to get the read-count matrix only for TCGA
  rse_traw <- scale_counts(SE.recount.tcga)
  traw <- assays(rse_traw)$counts
  #colData(ov.recount.tcga$TCGA_ovary)$gdc_cases.samples.portions.analytes.aliquots.submitter_id
  #colData(rse_traw)$gdc_cases.samples.portions.analytes.aliquots.submitter_id
  colnames(traw)<-colData(rse_traw)$gdc_cases.samples.portions.analytes.aliquots.submitter_id
  
  # to get the read-count matrix only for TCGA
  rse_gnorm <- scale_counts(SE.recount.gtex)
  graw<- assays(rse_gnorm)$counts
  
  
  ### Here starts the SURFACER approach on GitHub----
  
  # Intersect common genes
  traw<-traw[rowSums(traw)>100,]
  graw<-graw[rowSums(graw)>100,]
  common<-intersect(rownames(graw),rownames(traw))
  traw<-traw[common,]
  graw<-graw[common,]
  
  
  #### Here we need to subset subtypes
  # Now we need to subset the specific TCGA tumor subtype
  subtypes<-unique(recount.tcga[[name]]$gdc_cases.project.project_id)
  if (is.null(subtypes)==FALSE){
    for (subtype in subtypes) {
      subset <- which(recount.tcga[[name]]$gdc_cases.project.project_id==subtype)
      traw2<-traw[,subset]
      # What is normal, what is cancer
      # TCGA
      tcgacodes<-substr(colnames(traw2),14,16)
      table(tcgacodes) # 01A 954, 01B 14
      tumor<-traw2[,tcgacodes=="01A"]
      colnames(tumor) <- substr(colnames(tumor),1,16)
      idx<-which(duplicated(colnames(tumor)))
      ifelse(length(idx)>0,tumor<-tumor[,-(idx)],tumor<-tumor)
      
      #GTEx
      gtexcodes <- substring(colnames(graw),1,1)
      table(gtexcodes)
      gnorm <- graw[,gtexcodes=="S"]
      idx<-which(duplicated(colnames(gnorm)))
      ifelse(length(idx)>0,gnorm<-gnorm[,-(idx)],gnorm<-gnorm)
      rawcounts <- cbind(tumor,gnorm)
      group<-factor(c(rep("tumor",ncol(tumor)),rep("normal",ncol(gnorm))))
      names(group)<-colnames(rawcounts)
      group<-relevel(group,ref="normal")
      
      # convert ENsembl gene names to HUGO Symbols
      
      ensgenes<-sub('\\..*',"",rownames(rawcounts))
      idx<-which(duplicated(ensgenes))
      ifelse(length(idx)>0,rawcounts<-rawcounts[-(idx),],rawcounts<-rawcounts)
      ensgenes<-sub('\\..*',"",rownames(rawcounts))
      idx<-which(duplicated(ensgenes)) # empty
      rownames(rawcounts)<-ensgenes
      
      # create a convlist for the squishy thing
      
      symbol_list<-as.character(ens2sym(ensgenes))
      names(symbol_list)<-ensgenes
      input<-as.matrix(rawcounts[,1:ncol(rawcounts)])
      rawmat<-squish(input,symbol_list,method="sum")
      
      # VST
      #expmat<-vst(rawcounts)
      ## Use DESeq2 function
      if(!file.exists(paste0("results/000_",name,"_",subtype,"-expmat.rda"))){
        expmat <- vst(rawmat, blind = TRUE, nsub = 1000, fitType = "parametric")
        save(expmat,rawmat,gnorm,tumor,file=paste0("results/000_",name,"_",subtype,"-expmat.rda"))
      } else {load(paste0("results/000_",name,"_",subtype,"-expmat.rda"))}
      
      
      # PCA to see what is tumor what is normal
      res.pca <- prcomp(t(expmat), scale = TRUE)
      gp<-fviz_pca_ind(res.pca,
                       geom = "point",
                       col.ind = group, # color by groups
                       palette = c("#00AFBB",  "#FC4E07"),
                       addEllipses = TRUE, # Concentration ellipses
                       ellipse.type = "confidence",
                       legend.title = "Tissue",
                       title = paste0(subtype," vs. GTEx"),
                       repel = TRUE)
      png(paste0("plots/000_pca",name,"_",subtype,".png"),w=1500,h=1500,res=300)
      print(gp)
      dev.off()
      
      if(!file.exists(paste0("data/GTEx_",name,"_",subtype,"-regulon.rda"))){
        # create the surface activity network for normal breast tissue
        regulon<-corto(expmat[,colnames(gnorm)],centroids=surfacer,nbootstraps = 1000,p=1e-08,nthreads=20,
                       verbose = TRUE)
        save(regulon,file = paste0("data/GTEx_",name,"_",subtype,"-regulon.rda"))
      } else {load(paste0("data/GTEx_",name,"_",subtype,"-regulon.rda"))}
      
      # MRA Tumor vs. Normal
      if(!file.exists(paste0("results/000_tcga_",name,"_",subtype,"-mra.rda"))){
        newcounts<-expmat[rowVars(expmat)>0.01,]
        mr<-mra(newcounts[,colnames(tumor)],newcounts[,colnames(gnorm)],regulon=regulon,minsize=15,nperm=1000,nthreads = 20,verbose=TRUE)
        save(mr,file=paste0("results/000_tcga_",name,"_",subtype,"-mra.rda"))
      } else {load(paste0("results/000_tcga_",name,"_",subtype,"-mra.rda"))}
      gc()
    }
  }else {
    # What is normal, what is cancer
    # TCGA
    tcgacodes<-substr(colnames(traw),14,16)
    table(tcgacodes) # 01A 954, 01B 14
    tumor<-traw[,tcgacodes=="01A"]
    colnames(tumor) <- substr(colnames(tumor),1,16)
    #GTEx
    gtexcodes <- substring(colnames(graw),1,1)
    table(gtexcodes)
    gnorm <- graw[,gtexcodes=="S"]
    rawcounts <- cbind(tumor,gnorm)
    group<-factor(c(rep("tumor",ncol(tumor)),rep("normal",ncol(gnorm))))
    names(group)<-colnames(rawcounts)
    group<-relevel(group,ref="normal")
    
    # convert ENsembl gene names to HUGO Symbols
    
    ensgenes<-sub('\\..*',"",rownames(rawcounts))
    idx<-which(duplicated(ensgenes))
    ifelse(length(idx)>0,rawcounts<-rawcounts[-(idx),],rawcounts<-rawcounts)
    ensgenes<-sub('\\..*',"",rownames(rawcounts))
    idx<-which(duplicated(ensgenes)) # empty
    rownames(rawcounts)<-ensgenes
    
    # create a convlist for the squishy thing
    
    symbol_list<-as.character(ens2sym(ensgenes))
    names(symbol_list)<-ensgenes
    input<-as.matrix(rawcounts[,1:ncol(rawcounts)])
    rawmat<-squish(input,symbol_list,method="sum")
    
    # VST
    #expmat<-vst(rawcounts)
    ## Use DESeq2 function
    if(!file.exists(paste0("results/000_",name,"-expmat.rda"))){
      expmat <- vst(rawmat, blind = TRUE, nsub = 1000, fitType = "parametric")
      save(expmat,rawmat,gnorm,tumor,file=paste0("results/000_",name,"-expmat.rda"))
    } else {load(paste0("results/000_",name,"-expmat.rda"))}
    
    
    # PCA to see what is tumor what is normal
    res.pca <- prcomp(t(expmat), scale = TRUE)
    gp<-fviz_pca_ind(res.pca,
                     geom = "point",
                     col.ind = group, # color by groups
                     palette = c("#00AFBB",  "#FC4E07"),
                     addEllipses = TRUE, # Concentration ellipses
                     ellipse.type = "confidence",
                     legend.title = "Tissue",
                     title = paste0(name," vs. GTEx"),
                     repel = TRUE)
    png(paste0("plots/000_pca",name,".png"),w=1500,h=1500,res=300)
    print(gp)
    dev.off()
    
    if(!file.exists(paste0("data/GTEx_",name,"-regulon.rda"))){
      # create the surface activity network for normal breast tissue
      regulon<-corto(expmat[,colnames(gnorm)],centroids=surfacer,nbootstraps = 1000,p=1e-08,nthreads=20,
                     verbose = TRUE)
      save(regulon,file = paste0("data/GTEx_",name,"-regulon.rda"))
    } else {load(paste0("data/GTEx_",name,"-regulon.rda"))}
    
    # MRA Tumor vs. Normal
    if(!file.exists(paste0("results/000_tcga_",name,"-mra.rda"))){
      newcounts<-expmat[rowVars(expmat)>0.01,]
      mr<-mra(newcounts[,colnames(tumor)],newcounts[,colnames(gnorm)],regulon=regulon,minsize=15,nperm=1000,nthreads = 20,verbose=TRUE)
      save(mr,file=paste0("results/000_tcga_",name,"-mra.rda"))
    } else {load(paste0("results/000_tcga_",name,"-mra.rda"))}
    gc()
  } 
  
}
