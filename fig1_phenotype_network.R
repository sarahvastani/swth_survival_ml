############################################
#Build a network of juvenile phenotypes
############################################

#load packages
pkgs<-c("tidyverse","scico","ggpubr","qgraph")
sapply(pkgs,library,character.only=T); rm(pkgs)
theme_set(theme_classic())

#setwd("C:/Users/Steph/OneDrive - Texas A&M University/Thrushes/survival/data/")


phenotypes<-knn1<-read_csv("C:/Users/Steph/GitHub/thrush_hybrids/migratory_traits/knn_imputed_full_dataset.csv")

colSums(!is.na(phenotypes))

phenoInput<-phenotypes%>%
  select(-t5_spring40,-set)%>%
  drop_na()


#----------------------------------
## Functions from Wilkins et al. 2015
#A. Function to get confidence intervals for bootstraps of correlation matrix

boot.cor<-function(DATA,shps,iterations)#DATA is correlation matrix, shps are modality designations, iterations is number of bootstraps
{
  t=iterations # number of iterations for bootstrapping; 
  pb <- txtProgressBar(min=1, max=t-1,style=3) #define progress bar
  cor.orig=cor(DATA,method="spearman",use="pairwise.complete.obs")
  diag(cor.orig)=0
  boot.cor=array(dim=c(nrow(cor.orig),ncol(cor.orig),t))
  cors=vector(length=t)
  #Loop to create t=iterations datasets composed of k trait observations from n individuals resampled with replacement from our original dataset
  for (i in 1:t)
  {
    setTxtProgressBar(pb,i)#update progress bar
    s=sample(1:nrow(DATA),nrow(DATA),replace=T) #sample rows with replacement
    new.dat=DATA[s,] #use the row numbers from above to construct bootstrapped sample
    m=cor(new.dat,method="spearman",use="pairwise.complete.obs") #new correlation matrix
    boot.cor[,,i]=m
    cors[i]=cor.test(m[which(upper.tri(m)==TRUE)],cor.orig[which(upper.tri(cor.orig)==TRUE)])$estimate
    #       assort.boot.orig[i]=assortment.discrete(abs(m),shps)$r #calculate assortativity. Use absolute values.
    #       assort.boot.filt[i]=assortment.discrete(abs(m.filt),shps)$r #calculate assortativity for filtered bootstrap data matrix
  }
  
  avg.cor=apply(boot.cor,c(1,2),mean,na.rm=T) #calculate average correlation matrix from 1000 resampling iterations
  #calculate bootstrap confidence intervals for each edge in the original dataset
  lower.cor=apply(boot.cor, c(1,2), quantile, probs=c(0.025))
  upper.cor=apply(boot.cor, c(1,2), quantile, probs=c(0.975))
  return(list(orig.cor=cor.orig,avg.cor=avg.cor,lower.cor=lower.cor,upper.cor=upper.cor))
}
######### End bootstrapping function


######
# B. Function to plot bootstrap correlation CIs vs empirical correlations, and extract robust edges

#Analyze output of boot.cor function, generating assorativity and redundancy metrics
#Also generates figure showing bootstrap correlations against empirical correlations

analyze.boot<-function(boot.cor.output)
{
  par(cex.lab=1.8,mar=c(5,5,3,3),cex.axis=1.4,font.lab=2,mgp=c(3.5,1,0))  
  D<-boot.cor.output
  emp.edges=D$orig.cor[which(upper.tri(D$orig.cor)==TRUE)]
  boot.avgs=D$avg.cor[which(upper.tri(D$avg.cor)==TRUE)]
  boot.lower=D$lower.cor[which(upper.tri(D$lower.cor)==TRUE)]
  boot.upper=D$upper.cor[which(upper.tri(D$upper.cor)==TRUE)]
  ci.0=c("black","red")[((boot.lower>0&boot.upper>0)|(boot.lower<0&boot.upper<0))+1]
  plot(emp.edges,boot.avgs, pch=19, ylim=c(-1,1), xlim=c(-1,1), xlab="Empirical Correlations", ylab="Bootstrap Correlations",col=ci.0,las=1)
  for (i in 1:length(emp.edges)){
    lines(c(emp.edges[i], emp.edges[i]), c(boot.lower[i],boot.upper[i]),col=ci.0[i])
  }
  min.threshold<-
    abline(v=0.3,lty=2)
  abline(v=-0.3,lty=2)
  abline(h=0, lty=1)
  
  # find threshold where bootstrap CIs don't overlap 0.
  abs(boot.lower)>0
  mini=maxi=0
  for(i in 1: length(emp.edges))
  {
    if(emp.edges[i]<mini&ci.0[i]=="black"){mini=emp.edges[i]
    }else if(emp.edges[i]>maxi&ci.0[i]=="black"){maxi=emp.edges[i]
    }else{}
  }
  
  ci.0_full=c("black","red")[((D$lower.cor>0&D$upper.cor>0)|(D$lower.cor<0&D$upper.cor<0))+1]  
  robust.indx<-which(ci.0_full=="red")
  nonrobust.indx<-which(ci.0_full=="black")
  
  return(list(thresh.min=mini,thresh.max=maxi,robust.indx=robust.indx,nonrobust.indx=nonrobust.indx))
}# End analyze.boot function

## END FUNCTIONS


# specify node shapes
# divide into groups if comparing among modalities
phenoGroupNames=rep("trait1",ncol(phenoInput))

#use 100000 for final version
pheno.boot.cor.output<-boot.cor(phenoInput,phenoGroupNames,iterations=100000)
pheno.boot.analy<-analyze.boot(pheno.boot.cor.output) #analyze and plot bootstrap output

#Get correlation matrix for unfiltered dataset
pheno.nodes<-cor(phenoInput,method="spearman")
# Remove nonrobust edges from graph
pheno.nodes[pheno.boot.analy$nonrobust.indx]<-0
diag(pheno.nodes)<-0


renameFeatures1<-
  data.frame(Feature=c("release_year","fall_detectDay1","aims_heterozygosity",
                       "bodyCondition","sex_binary","fall_bearing1","kipps",
                       "tail.length","distal","releaseDay","wing.cord","p10",
                       "aims_ancestry","tarsus.length"),
             FeatureRenamed=c("Release year","Fall timing","Heterozygosity",
                              "Body condition","Sex","Fall orientation","Kipps (wing)",
                              "Tail length","Distal (wing)","Release day","Wing cord","p10 (feather)",
                              "Ancestry","Tarsus length"),
             FeatureCategory=c("Other","Behaviour","Genetics",
                               "Morphology","Other","Behaviour","Morphology",
                               "Morphology","Morphology","Behaviour","Morphology","Morphology",
                               "Genetics","Morphology"))
Features1<-data.frame(Feature=colnames(phenoInput))%>%
  left_join(renameFeatures1)%>%
  left_join(data.frame(
    FeatureCategory=sort(unique(renameFeatures1$FeatureCategory)),
    FeatureColours=scico(4,palette='devon',categorical=F,end=0.92)))

#png(filename="../figures/phenoNetwork.png",width=600,height=400)
qgraph( pheno.nodes,
        color="white",layout="spring",
        labels=Features1$FeatureRenamed,
        minimum=0.1, #minimum edge to plot
        vsize=12, #node size
        label.norm="0000", #normalize node labels to string length
        fade=F,shape="circle",label.scale=T,
        border.color=Features1$FeatureColours,border.width=8,
        negCol="grey60",posCol="grey20",label.color=1,cut=0 )

dev.off()

write_csv(as.data.frame(pheno.nodes)%>%rownames_to_column('Feature'),
          file="C:/Users/Steph/GitHub/thrush_hybrids/migratory_traits/phenoNetworkCorr.csv")
