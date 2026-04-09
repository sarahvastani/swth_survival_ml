library(tidyverse)
library(scico)
library(sf)
library(rnaturalearth)
library(qgraph)
library(cowplot)
theme_set(theme_classic())

meta1<-read.csv("~/GitHub_data/AIMs_metaData/AIMs_metadata.20250327.csv")
lat.df<-read.csv("~/GitHub_data/survivalData/motusDetections.filtered.20250409.csv")
receivers<-read.csv("~/GitHub_data/survivalData/receiver_database_dec25.csv")
data1<-read.csv("C:/Users/Steph/GitHub/thrush_hybrids/migratory_traits/Dataset1.phiSurvival.20250707.csv")

lat.df.map<-lat.df%>%filter(name_in_vcf%in%data1$name_in_vcf)%>%
  mutate(detectDate=substr(ts,1,10))%>%
  dplyr::select(recvDeployName,recvDeployLat,recvDeployLon,detectdoy,detectyear,
                detectDate,name_in_vcf,aims_ancestry)%>%
  arrange(detectDate)%>%
  unique()
lat.df.map$aims_ancestry
meta1Release<-meta1%>%
  filter(name_in_vcf%in%data1$name_in_vcf)%>%
  mutate(recvDeployName="releaseSitePemberton")%>%
  mutate(detectdoy=data.table::yday(release_date))%>%
  mutate(recvDeployLat=as.numeric(release_gps.n),recvDeployLon=as.numeric(release_gps.w))%>%
  rename(detectyear=release_year,detectDate=release_date)%>%
  select(all_of(colnames(lat.df.map)))

lat.df.map2<-rbind(lat.df.map,meta1Release)%>%
  arrange(detectDate)

releaseSite<-meta1Release%>%
  summarise(meanLat=mean(recvDeployLat),meanLon=mean(recvDeployLon))

data1%>%
  mutate(ancestryGroup=case_when(aims_ancestry<=0.1~"inland",
                                 aims_ancestry>=0.9~"coastal",
                                 aims_ancestry>0.1&aims_ancestry<0.5~"hybrid_inland",
                                 aims_ancestry>=0.5&aims_ancestry<0.9~"hybrid_coastal"))%>%
  group_by(ancestryGroup)%>%summarise(count=n())

data1%>%
  mutate(ancestryGroup=case_when(aims_ancestry<=0.1~"inland",
                                 aims_ancestry>=0.9~"coastal",
                                 aims_ancestry>0.1&aims_ancestry<0.5~"hybrid_inland",
                                 aims_ancestry>=0.5&aims_ancestry<0.9~"hybrid_coastal"))%>%
  group_by(ancestryGroup,phi_binary)%>%summarise(count=n())%>%
  ungroup()%>%group_by(ancestryGroup)%>%mutate(percent=count/sum(count))%>%
  filter(phi_binary==1)


#load map data
world <- ne_countries(scale = "medium", returnclass = "sf")

receivers<-receivers%>%
  filter(timeStart<as.Date("2025-12-31")&timeEnd>("2019-01-01"))

map1<-ggplot() +
  geom_sf(data = world,fill="grey95",colour='grey80')+
  geom_point(data=receivers,aes(x=longitude,y=latitude),
             size=0.5,colour='grey60',shape=1)+
  geom_point(data=lat.df.map,aes(x=recvDeployLon,y=recvDeployLat,
                                 group=name_in_vcf,colour=aims_ancestry),
             size=1.2,alpha=0.5)+
  geom_path(data=lat.df.map,aes(x=recvDeployLon,y=recvDeployLat,
                                  group=name_in_vcf,colour=aims_ancestry),
            size=0.1,alpha=0.5)+
  geom_point(data=releaseSite,aes(x=meanLon,y=meanLat),colour="grey10",size=2,shape=8)+
  scale_colour_scico(palette="devon",end=0.7,direction=-1,name="Ancestry")+
  xlab("Longitude")+ylab("Latitude")+
  #theme(legend.position="none")+
  ylim(0,66)+xlim(-140,-55)+
  theme(legend.position=c(0.92,0.8),legend.background=element_rect(fill="white"))
map1
pheno.nodes<-read_csv("C:/Users/Steph/GitHub/thrush_hybrids/migratory_traits/phenoNetworkCorr.csv")%>%
  column_to_rownames("Feature")%>%
  as.matrix()

renameFeatures1<-
  data.frame(Feature=c("release_year","fall_detectDay1","aims_heterozygosity",
                       "bodyCondition","sex_binary","fall_bearing1","kipps",
                       "tail.length","distal","releaseDay","wing.cord","p10",
                       "aims_ancestry","tarsus.length"),
             FeatureRenamed=c("Release\nyear","Fall\ntiming","Hetero-\nzygosity",
                              "Body\ncondition","Sex","Fall\norien-\ntation","Kipps\n(wing)",
                              "Tail\nlength","Distal\n(wing)","Release\nday","Wing\ncord","p10\n(feather)",
                              "Ancestry","Tarsus\nlength"),
             FeatureCategory=c("Other","Behaviour","Genetics",
                               "Morphology","Other","Behaviour","Morphology",
                               "Morphology","Morphology","Behaviour","Morphology","Morphology",
                               "Genetics","Morphology"))
ColDf1<-data.frame(
  FeatureCategory=sort(unique(renameFeatures1$FeatureCategory)),
  FeatureColours=scico(4,palette='devon',categorical=F,end=0.92))
Features1<-data.frame(Feature=colnames(pheno.nodes))%>%
  left_join(renameFeatures1)%>%
  left_join(ColDf1)
network1<-qgraph( pheno.nodes,
                  color="white",layout="spring",
                  labels=Features1$FeatureRenamed,
                  minimum=0.1, #minimum edge to plot
                  vsize=12, #node size
                  label.norm="0000", #normalize node labels to string length
                  fade=F,shape="circle",label.scale=T,
                  border.color=Features1$FeatureColours,border.width=8,
                  negCol="grey60",posCol="grey20",label.color=1,cut=0 )
# png("C:/Users/Steph/GitHub/thrush_hybrids/migratory_traits/network1.png",res=3000)
# plot(network1)
# dev.off()

Features1<-Features1%>%mutate(FeatureRenamed2=gsub("\n"," ",gsub("-\n","",FeatureRenamed)))

legend1=get_legend(ggplot(Features1,aes(x=FeatureCategory,fill=FeatureCategory))+
  geom_bar()+
  scale_fill_manual(values=ColDf1$FeatureColours,name=""))

network2<-ggdraw() +
  draw_image("C:/Users/Steph/GitHub/thrush_hybrids/migratory_traits/network1.png")
gg1<-ggpubr::ggarrange(map1,network2,legend1,labels=c("A","B",""),
                  nrow=1,ncol=3,widths=c(2,2,0.7))
gg1

ggsave("C:/Users/Steph/GitHub/thrush_hybrids/migratory_traits/Fig1.pdf",
       plot=gg1,
       height = 4,width=10,bg='white')
ggsave("C:/Users/Steph/GitHub/thrush_hybrids/migratory_traits/Fig1.png",
       plot=gg1,dpi=1000,
       height = 4,width=10,bg='white')
