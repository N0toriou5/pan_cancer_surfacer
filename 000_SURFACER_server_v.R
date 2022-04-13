# SURFACER core pipeline to run on a server
# The PAN-Cancer Surfaceome Landscape . generate all MRAs. 
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
  ### This code part summarizes step 1 and 2
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
  colnames(traw)<-colData(rse_traw)$gdc_cases.samples.portions.analytes.aliquots.submitter_id
  
  # to get the read-count matrix only for TCGA
  rse_gnorm <- scale_counts(SE.recount.gtex)
  graw<- assays(rse_gnorm)$counts
  
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
      
      ### This code part performs step 4
      if(!file.exists(paste0("data/GTEx_",name,"_",subtype,"-regulon.rda"))){
        # create the surface activity network for normal breast tissue
        regulon<-corto(expmat[,colnames(gnorm)],centroids=surfacer,nbootstraps = 1000,p=1e-08,nthreads=20,
                       verbose = TRUE)
        save(regulon,file = paste0("data/GTEx_",name,"_",subtype,"-regulon.rda"))
      } else {load(paste0("data/GTEx_",name,"_",subtype,"-regulon.rda"))}
      
      # MRA Tumor vs. Normal (step 5)
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

### Differential Expression Analysis (step 3)
filenames<-dir(path="results/",pattern = "-expmat.rda",full.names=TRUE) # eliminare il 3 LAML
delist<-list()

for (i in 1:length(filenames)){
  file<-filenames[i]
  load(file)
  name<-str_match(file,"-\\s*(.*?)\\s*-")[2]
 # DESeq2 block (filter out poorly expressed genes)
  group<-factor(c(rep("tumor",ncol(tumor)),rep("normal",ncol(gnorm))))
  names(group)<-colnames(rawmat)
  group<-relevel(group,ref="normal")
  design<-as.data.frame(group)
  colnames(design)<-"cond"
  dds<-DESeqDataSetFromMatrix(countData=rawmat,colData=design,design=~cond)
  dds<-dds[rowSums(counts(dds))>=5,]
  dds$cond<-relevel(dds$cond,ref="normal")
  dea<-DESeq(dds,parallel=TRUE)
  resultsNames(dea) #"cond_tumor_vs_normal"
  res<-as.data.frame(results(dea,name="cond_tumor_vs_normal"))
  res<-res[order(res$pvalue),]
  common<-intersect(rownames(res),surfacer)
  res<-res[common,]
  write.xlsx(res, file=paste0("results/",name,"_DE.xlsx"))
  delist[[i]]<-res
  names(delist)[i]<-name
}
save(delist,file="results/000_DEA.rda")

#### The following code part performs step 6 (it identifies SSHs)
filenames<-dir(path="results/",pattern = "-mra",full.names=TRUE)
mralist<-list()
for (i in 1:length(filenames)){
  name<-filenames[i]
  load(name)
  subtype <- str_match(name, "-\\s*(.*?)\\s*-")[2]
  mralist[[i]]<-mr
  names(mralist)[i]<-subtype
}
length(mralist) #20 TCGA tumor subtypes
save(mralist,file = "results/000_MRA.rda")
surfacer<-Reduce(intersect, list(names(which(abs(mralist[[1]]$nes)>=1.96)),
                                 names(which(abs(mralist[[2]]$nes)>=1.96)),
                                 names(which(abs(mralist[[3]]$nes)>=1.96)),
                                 names(which(abs(mralist[[4]]$nes)>=1.96)),
                                 names(which(abs(mralist[[5]]$nes)>=1.96)),
                                 names(which(abs(mralist[[6]]$nes)>=1.96)),
                                 names(which(abs(mralist[[7]]$nes)>=1.96)),
                                 names(which(abs(mralist[[8]]$nes)>=1.96)),
                                 names(which(abs(mralist[[9]]$nes)>=1.96)),
                                 names(which(abs(mralist[[10]]$nes)>=1.96)),
                                 names(which(abs(mralist[[11]]$nes)>=1.96)),
                                 names(which(abs(mralist[[12]]$nes)>=1.96)),
                                 names(which(abs(mralist[[13]]$nes)>=1.96)),
                                 names(which(abs(mralist[[14]]$nes)>=1.96)),
                                 names(which(abs(mralist[[15]]$nes)>=1.96)),
                                 names(which(abs(mralist[[16]]$nes)>=1.96)),
                                 names(which(abs(mralist[[17]]$nes)>=1.96)),
                                 names(which(abs(mralist[[18]]$nes)>=1.96)),
                                 names(which(abs(mralist[[19]]$nes)>=1.96)),
                                 names(which(abs(mralist[[20]]$nes)>=1.96))))

### Surface Biomarkers prioritization (step 7)
load("results/000_DEA.rda")
mytargets<-c() #### mytargets are my SSHs!
for (name in surfacer){
  ev<-TRUE
  for (j in 1:length(delist)){
    ifelse(abs(delist[[j]][name,]$log2FoldChange)>=0.5,ev<-TRUE,ev<-FALSE)
  }
  ifelse(ev==TRUE,mytargets<-c(mytargets,name),mytargets<-mytargets)  
}

### Survival Analysis (step 8)

### one gene survival
filenames<-dir(path="results/",pattern = "-expmat",full.names=TRUE)
subtypes<- str_match(filenames, "-\\s*(.*?)\\s*-")
subtypes<-gsub("-","",subtypes)[,1]
subtypes<-sort(subtypes)
survmat<-matrix(ncol=20,nrow=length(mytargets),dimnames = list(mytargets,subtypes))
for (subtype in subtypes){
  load(paste0("data/tcga_",subtype,"-survival.rda"))
  file<-filenames[grep(subtype, filenames)]
  load(file)
  common<-intersect(colnames(expmat),names(survival))
  survival<-survival[common]
  expmat<-expmat[,common]
  for (gene in mytargets){
    mygene<-gene
    oritrack<-expmat[mygene,]
    track<-oritrack
    track[]<-"Other"
    track[oritrack>median(oritrack)]<-paste0(mygene,"high")
    track[oritrack<=median(oritrack)]<-paste0(mygene,"low")
    track<-as.factor(track)
    track<-relevel(track,ref=paste0(mygene,"low"))
    
    
    # Fit model
    sdiff<-survdiff(survival~track)
    l<-sdiff$obs
    m<-(l[1]+0.1)/(l[2]+0.1)
    p<-1-pchisq(sdiff$chisq, df=length(sdiff$n) - 1)
    events_difference<-sdiff$obs-sdiff$exp
    if(events_difference[1]>events_difference[length(events_difference)]){
      sign<-"neg"
    } else {
      sign<-"pos"
    }
    if( p<=0.05){
      ifelse(sign=="pos",survmat[mygene,subtype]<- (-log10(p)),survmat[mygene,subtype]<- (-(-log10(p))))
      png(paste0("plots/003_",subtype,"_",mygene,".png"),w=1500,h=1500,res=300)
      #plotstepsurv(oritrack = oritrack,survival = survival,mygene=mygene,ngroups=4,cox_survival = TRUE,title=mygene)
      plotclassicsurv(oritrack = oritrack,survival = survival,mygene=mygene,title=paste(subtype,mygene))
      dev.off()
    }else
      {next}
    }

}

### Sample-by-sample survival
filenames<-dir(path="results/",pattern = "-expmat.rda",full.names=TRUE)
filenames<-filenames[-c(9:11)]

for (file in filenames){
  load(file)
  name<-gsub(".+TCGA_","",file)
  name<-sub("-expmat.rda","-regulon.rda",name)
  load(paste0("data/GTEx_TCGA_",name))
  # Run Surface Markers Classifier
  tum<-expmat[,colnames(tumor)]
  tum<-tum[rowVars(tum)>0.01,]
  spmra<-mra(tum,regulon=regulon,nthreads=6,verbose=TRUE) # single-patient MRA
  name<-sub("-regulon.rda","",name)
  save(spmra,file=paste0("results/000_",name,"-spmra.rda"))
}

dist2d <- function(a,b,c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
} 


filenames<-dir(path="results/",pattern = "-spmra.rda",full.names=TRUE)
for (file in filenames){
  load(file)
  name<-gsub(".+TCGA_","",file)
  name<-sub("-spmra.rda","",name)
  # Plot SP-MRA
  genes<-intersect(rownames(spmra),surfacer)
  # Create the heatmap
  ## Prepare table
  toshow<-spmra[genes,]
  # Pairwise correlation between samples (columns)
  cols.cor <- cor(toshow, use = "pairwise.complete.obs", method = "pearson")
  clustering_distance_cols = as.dist(1 - cols.cor)
  res.hc <- hclust(clustering_distance_cols, method = "ward.D2" )
  subtype <- sub(".+TCGA-","",name)
  
  gp<-fviz_nbclust(t(toshow), kmeans, diss=clustering_distance_cols,method = "wss")
  
  distances<-c()
  for (i in c(1:10)){
    p2<-c(1,gp$data$y[1])
    p3<-c(10,gp$data$y[10])
    p1<-c(i,gp$data$y[i])
    diss<-dist2d(p1,p2,p3)
    distances<-c(distances,diss)
  }
  
  k<-which.max(distances)
  #fviz_nbclust(t(toshow), kmeans, diss=clustering_distance_cols,method = "gap_stat",nboot=100)
  # fviz_nbclust(t(toshow), kmeans, diss=clustering_distance_cols,method = "silhouette")
  png(paste0("plots/004_",subtype,"_elbow.png"),w=2000,h=1500,res=450)
  plotting<-gp+geom_vline(xintercept = k, linetype = 2)+
    labs(subtitle = "Elbow method")
  print(plotting)
  dev.off()
  # Plot the heatmap
  colside<-cutree(res.hc,k=k)
  # Color function
  colfun<-colorRampPalette(c("navy","navy","blue","blue","white","red","red","red3","red3"))
  toplot<-toshow
  #toplot[toplot>10]<-10
  #toplot[toplot< -10]<-(-10)
  
  png(paste0("plots/004_",subtype,"_heatmap.png"),w=4000,h=3000,res=300)
  heatmap.3(toplot,KeyValueName = "NES",hcc=res.hc,col=colfun,ColSideColors = colside)
  dev.off()
  
  # survival
  
  load(paste0("data/tcga_",subtype,"-survival.rda"))
  #load("F:/Projects/Survivals/TCGA/BRCA/survival.rda")
  patients<-intersect(names(colside),names(survival))
  survival<-survival[patients]
  colside<-colside[patients]
  Subtype<-colside
  sfit<-survfit(survival~Subtype)
  png(paste0("plots/004_survival_surfaceR_",subtype,".png"),w=2500,h=2500,res=350)
  gp<-ggsurvplot(sfit, data=survival,
             conf.int=FALSE, pval=TRUE, risk.table=F, size = 1.2,censor.size=3,
             title=paste0("Kaplan-Meier Curve for TCGA-",subtype," Survival"))
  print(gp)
  dev.off()
}

### Druggable genome analysis (step 9)


