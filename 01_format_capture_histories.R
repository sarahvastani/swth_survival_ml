###################################################################
#Clean up raw radio tag data to estimate survival
###################################################################

#load packages
pkgs<-c("tidyverse","viridis","ggpubr","sf","rnaturalearth","readxl")
sapply(pkgs,library,character.only=T); rm(pkgs)
theme_set(theme_classic())

setwd("C:/Users/Steph/TAMU_OneDrive/Thrushes/survival/data/raw/")

#This is the output from pre-cleaning scripts 1 and 2
#see 04_01_download_motus.R and 04_02_add_motus_tag_ids.R for additional filtering
#lat.df.Jul23<-read.csv("2.filtered_data_Jul_1.csv")
lat.df.raw<-read_csv("C:/Users/Steph/GitHub_data/survivalData/1.filtered_motus_data_dec_21_24.csv")
lat.df<-lat.df.raw

#read in meta data - output from 01_phenotypes_cleaning.R
meta_thrush<-read_csv("C:/Users/Steph/GitHub_data/AIMs_metaData/AIMs_metadata.20250327.csv")
meta_thrush<-meta_thrush%>%select(!starts_with("aims_ancestry_")&!starts_with("aims_heterozygosity_"))

lat.df<-lat.df%>%dplyr::rename(motustagid=motusTagID)
#only verlapping column name should be motustagid
names(lat.df)[names(lat.df)%in%names(meta_thrush)]

#load map data
world <- ne_countries(scale = "medium", returnclass = "sf")

poly.all<-data.frame(X.co=rep(c(0,0.5,1)), #make triangle polygon for plotting
                     Y.co=rep(c(0,1,0)))

dupTag<-meta_thrush%>%filter(!is.na(motustagid))%>%
  mutate(motustagid2=as.character(motustagid))%>%
  filter(duplicated(motustagid2))%>%pull(motustagid)
meta_thrush%>%filter(motustagid==dupTag)

#CI07H04 was tagged in 2021, promptly flew into a window, and DH27H01 was given the same tag in 2022
#rename motustagid for CI07H04 from 58075 to 58075a

meta_thrush<-meta_thrush%>%
  mutate(motustagid=as.character(motustagid))%>%
  mutate(motustagid=if_else(name_in_vcf=="CI07H04_S118_L001",paste(motustagid,"a",sep=""),motustagid))

lat.df<-lat.df%>%
  mutate(motustagid=as.character(motustagid))%>%
  mutate(motustagid=if_else(motustagid=="58075"&detectyear==2021,"58075a",motustagid))

#merge detection with relevant metadata
lat.df<-lat.df%>%
  left_join(meta_thrush%>%dplyr::select(reference,name_in_vcf,capture_site,release_site,
                                               release_year,motustagid,aims_ancestry,aims_heterozygosity,
                                        release_gps.n,release_gps.w))

ggplot()+
  geom_polygon(data=poly.all,aes(x=X.co,y=Y.co),fill=NA,colour="black")+
  geom_point(data=meta_thrush%>%filter(reference%in%lat.df$reference),
             aes(x=aims_ancestry,y=aims_heterozygosity),size=1,shape=1,stroke=1)+
  coord_equal()+xlab("ancestry")+ylab("heterozygosity")



#clean up detections
#based on manual visual surveys of routes
lat.df<-lat.df%>%
  
  filter(!is.na(tagDeployStart))%>%
  
  #remove towers with many faulty detections
  filter(!recvDeployName%in%c("LacEdouard-Champs","Lambs Gap")&recvDeployLon<(-70))%>%
  
  #unlikely detections based on time of year and location (north in winter / south in summer)
  filter(!(name_in_vcf=="BI03H07"&(detectdoy>292|detectdoy<291)))%>%
  filter(!name_in_vcf%in%c("AI01H17","AI02H10"))%>%
  filter(!(recvDeployName=="Estero Punta Banda"&detectdoy%in%334:335))%>%
  
  #unlikely detections based on other detections for the same bird
  #ex. inland migration route with one coastal tower - remove coastal tower
  filter(!(recvDeployName=="Golfo de Santa Clara - RV parl"&
             name_in_vcf%in%c("AI01H14","AI02H06","AI02H08","BH31H12","BI03H07","CI03H02_S84_L001","BI02H10",
                              "CH31H04_S78_L001")))%>%
  filter(!(name_in_vcf%in%c("AH28H01","AH28H13","AI01H14","AI02H06","AI02H08")&
             recvDeployName=="Drasher"))%>%
  filter(!(name_in_vcf=="CH30H05_S68_L001"&detectdoy>230&detectdoy<234))%>%
  filter(!(name_in_vcf%in%c("BI04H14","BI07H03")&recvDeployName=="Phoenix"))%>%
  filter(!(name_in_vcf%in%c("BH30H01","CH31H03_S77_L001","CH31H05_S79_L001")&
             recvDeployName=="Mackay Island NWR, NC"))%>%
  filter(!(name_in_vcf%in%c("CI03H05_S87_L001")&recvDeployName=="McGill_Bird_Observatory"))%>%
  filter(!(name_in_vcf=="DH31H08"&recvDeployName=="Hopkins Forest"))%>%
  filter(!(name_in_vcf=="AI01H14"&recvDeployName=="GA_SSI_DOCK"&detectyear==2021))%>%
  filter(!(name_in_vcf=="CI03H02_S84_L001"&recvDeployName%in%c("Florida Panther NWR, FL")))%>%
  filter(!(name_in_vcf=="DH29H07"&recvDeployName%in%c("Florida Panther NWR, FL")&substr(ts,1,7)=="2022-09"))%>%
  filter(!(name_in_vcf=="DH29H07"&recvDeployLon>-100&substr(ts,6,7)%in%c("07")))%>%
  filter(!(name_in_vcf=="DH29H14"&recvDeployName%in%c("Florida Panther NWR, FL")&substr(ts,1,7)=="2022-09"))%>%
  filter(!(name_in_vcf=="DH29H14"&recvDeployLon>-100&substr(ts,6,7)%in%c("06","07")))%>%
  filter(!(name_in_vcf=="DH30H05"&recvDeployName%in%c("Florida Panther NWR, FL")))%>%
  filter(!(name_in_vcf=="DH30H05"&recvDeployLon>-100&substr(ts,6,7)=="07"))%>%
  filter(!(name_in_vcf=="EH28H14"&recvDeployName=="Quilchena"))%>%
  filter(!(name_in_vcf=="EH25H07"&recvDeployName=="Buckner"&detectyear==2024))%>%
  filter(!(name_in_vcf=="EH27H08"&recvDeployName=="Kennekuk 6"&detectyear==2024))%>%
  filter(!(name_in_vcf=="EH27H13"&recvDeployName=="Ruthven"&detectyear==2024))%>%
  filter(!(name_in_vcf=="EH27H16"&recvDeployName=="Ruthven"&detectyear==2024))%>%
  filter(!(name_in_vcf=="EH29H02"&recvDeployName%in%c("Kennekuk Cove County Park","GFAFB - Lagoons")&detectyear==2024))%>%
  filter(!(name_in_vcf=="EH29H25"&recvDeployName%in%c("Ruthven","Buckner")&detectyear==2024))%>%
  filter(!(name_in_vcf=="EH31H01"&recvDeployName=="Buckner"&detectyear==2024))%>%
  filter(!(name_in_vcf=="FH30H07"&recvDeployName=="Buckner"&detectyear==2024))%>%
  filter(!(name_in_vcf=="FH30H16"&recvDeployName=="Buckner"&detectyear==2024))%>%
  filter(!(name_in_vcf=="FH31H06"&recvDeployName=="Buckner"&detectyear==2024))%>%
  filter(!(name_in_vcf=="FH31H08"&recvDeployName=="Buckner"&detectyear==2024))%>%
  filter(!(name_in_vcf=="FI01H01"&recvDeployName%in%c("Ruthven","Buckner")&detectyear==2024))%>%
  filter(!(name_in_vcf=="FI01H08"&recvDeployName=="GFAFB - Lagoons"&detectyear==2024))%>%
  filter(!(name_in_vcf=="FI02H06"&recvDeployName=="Ruthven"&detectyear==2024))%>%
  filter(!(name_in_vcf=="FI03H02"&recvDeployName=="Ruthven"&detectyear==2024))%>%
  filter(!(name_in_vcf=="FI03H06"&recvDeployName=="Buckner"&detectyear==2024))%>%
  filter(!(name_in_vcf=="FI03H15"&recvDeployName=="Buckner"&detectyear==2024))%>%
  filter(!(name_in_vcf=="FI04H02"&recvDeployName%in%c("Buckner","Ruthven")&detectyear==2024))%>%
  filter(!(name_in_vcf=="FH30H09"&recvDeployName=="Ruthven"&detectyear==2024))%>%
  
  #no detections prior to year 2
  filter(!(name_in_vcf%in%c("BI03H05_S181_L001","BI05H03_S185_L001")))%>%
  
  #manually set tag deploy start for tag that was activated in Costa Rica in March but deployed in BC in August
  mutate(tagDeployStart=if_else(name_in_vcf=="FH30H07",as_datetime("2024-08-30 12:00:00 UTC"),tagDeployStart))%>%
  
  
  #apply general day by location cutoffs
  filter(!(recvDeployLat<30&detectdoy>136&detectdoy<243))%>% #below 30N, June 15 to August 30
  filter(!(recvDeployLat>35&(detectdoy>334|detectdoy<59)))%>% #above 35N, Dec to Feb
  
  mutate(deployDay=substr(tagDeployStart,1,10),tsDay=substr(ts,1,10))%>%
  filter(!is.na(name_in_vcf))%>%
  mutate(tagDays=as.integer(gsub(" days"," ",difftime(tsDay,deployDay,units="days"))))%>%
  
  #remove detections for tag turned on before deployment (issue for FH30H07)
  filter(tagDays>0)


#note - discrepancies between banding sheets and tag latitudes
#banding sheet latitudes match listed release sites - replace tag deploy longitudes
lat.df%>%
  filter(tagDepLat!=release_gps.n)%>%
  dplyr::select(name_in_vcf,motustagid,tagDepLat,release_gps.n,release_site,capture_site)%>%
  distinct()

lat.df<-lat.df%>%
  mutate(tagDepLat=release_gps.n,
         tagDepLon=release_gps.w)

today=gsub("-","",Sys.Date())

#output csv of filtered detections for estimating phenotypes & mapping

# write.csv(lat.df,
#           paste("C:/Users/Steph/GitHub_data/survivalData/motusDetections.filtered",today,"csv",sep="."),
#           row.names=F)

###################################################################
#format capture histories for input to CJS models
###################################################################


ch.df<-data.frame()
ch_multi.df<-data.frame()

#for each individual in juvenile survival dataset
#released in Pemberton, not translocated
for(i in meta_thrush%>%filter(tag_type=="radio"&#age_release=="HY"&
                              ((release_site=="Pemberton"&capture_site=="Pemberton")|
                              (release_site=="Whistler"&capture_site=="Whistler")))%>%
    filter(!name_in_vcf%in%c("BH30H01","BH29H01","BH29H03","BH29H04","BH29H05","BH29H07"))%>%
    pull(name_in_vcf)){
  x1<-rep(0,300)
  x2<-rep(0,300)
  
  y<-lat.df%>%filter(name_in_vcf==i)%>%
    mutate(ch_multi=case_when(recvDeployLat>40&tagDays<150~"A",
                              recvDeployLat<40~"B",
                              recvDeployLat>40&tagDays>150~"C"))%>%
    select(name_in_vcf,tagDays,ch_multi)%>%distinct()
  
  x2[y%>%filter(tagDays<301)%>%pull(tagDays)]<-
    y%>%filter(tagDays<301)%>%pull(ch_multi)
  x2[1]<-"A"
  
  y1<-y%>%pull(tagDays)
  x1[c(1,y1[y1<301])]<-1 #add 1 to detections because bird was caught the day it was tagged
  if(length(y1)>0){if(max(y1)>300){
    x1[300]<-1
    x2[300]<-"C"    }} #fill in last day if bird caught later
  
  ch.df<-rbind(ch.df,c(i,x1))
  ch_multi.df<-rbind(ch_multi.df,c(i,x2))
  
}

colnames(ch.df)[1]<-"name_in_vcf"
colnames(ch_multi.df)[1]<-"name_in_vcf"

ch.df<-ch.df%>%unite(ch,-name_in_vcf,sep="")
ch_multi.df<-ch_multi.df%>%unite(ch_multi,-name_in_vcf,sep="")
ch.df<-ch.df%>%left_join(ch_multi.df)%>%left_join(meta_thrush)%>%
  mutate(ch=paste("days",ch,sep="_"))%>% #stops capture history from being saved as a very large number
  select(name_in_vcf,aims_ancestry,aims_heterozygosity,
         release_site,release_year,tag_type,sex_binary,age_release,
         release_gps.n,release_gps.w,motustagid,ch,ch_multi)

#group histories into sets of 10 days
#this is might be one approach to deal with clusters of detections

ch.df<-ch.df%>%mutate(ch2=gsub("days_","",ch))

ch2<-data.frame(name_in_vcf=ch.df$name_in_vcf)

days<-10

for(i in 1:round((nchar(ch.df$ch2)[1]/days))){
  colID=paste("t",i,sep="")
  ch2[colID]<-rep(0,nrow(ch.df))
  j=i*days-(days-1) #get starting position
  ch2[grep("1",substr(ch.df$ch2,j,j+(days-1))),colID]<-1
}

ch2<-ch2%>%unite(ch,-name_in_vcf,sep="")

ch2<-ch2%>%left_join(meta_thrush)%>%
  mutate(ch_10days=paste("days",ch,sep="_"))%>% #stops capture history from being saved as a very large number
  select(reference,name_in_vcf,motustagid,aims_ancestry,release_site,release_year,sex_binary,ch_10days)

today
#write.csv(ch2,paste("C:/Users/Steph/GitHub/thrush_hybrids/updateSurvival/thrush_survival_ch",today,"csv",sep="."))



###################################################################
##Look at individual towers and birds for filtering
#results applied above in lines 62-116
###################################################################


#print data for particular towers
tower.names<-c("Estero Punta Banda","Golfo de Santa Clara - RV parl")
tower.check<-lat.df%>%filter(recvDeployName%in%tower.names)%>%
   mutate(ts=substr(ts,1,10),tagDeployStart=substr(tagDeployStart,1,10))%>%
   group_by(name_in_vcf,tagDeployStart,ts,recvDeployName,recvDeployLat,recvDeployLon,motustagid)%>%
     summarise(nDetect=n(),max_run=max(runLen),mean_run=mean(runLen))

#check data for a particular bird
bird.names<-c("CI03H02_S84_L001")
bird.check<-lat.df%>%filter(name_in_vcf%in%bird.names)%>%
    group_by(name_in_vcf,ts,tagLifespan,tagDeployStart,recvDeployID,recvDeployLat,
             recvDeployLon,recvDeployName,detectyear,detectdoy)%>%
    mutate(ts=substr(ts,1,10))%>%mutate(tagDeployStart=substr(tagDeployStart,1,10))%>%
    summarise(nDetect=n(),max_run=max(runLen),mean_run=mean(runLen))

ggplot(data = world) +
  geom_sf(fill="gray98")+
  geom_point(data=bird.check,
             aes(x=recvDeployLon,y=recvDeployLat,colour=detectdoy,
                             shape=as.factor(detectyear)),size=3)+
  coord_sf(xlim=c(-150,-40),ylim=c(10,60),expand=FALSE)


options(dplyr.print_min = Inf) 

#adjust so that start day is Aug 24 (day first bird tagged) not Jan 1
lat.df<-lat.df%>%
  mutate(tag.days.yr=ifelse(tagDays>365,tagDays-365,tagDays),
         yr=case_when(tagDays>365&tagDays<365*2~"year2",tagDays<366~"year1",tagDays>365*2~"year3"))%>%
  mutate(doy.aug24=case_when(ts>paste(detectyear,"08-24",sep="-")~difftime(ts,paste(detectyear,"08-24",sep="-"),units="days"),
                             ts<=paste(detectyear,"08-24",sep="-")~difftime(ts,paste((detectyear-1),"08-24",sep="-"),units="days")))%>%
  mutate(doy.aug24=round(as.numeric(gsub("days ","",doy.aug24))))%>%
  mutate(tag.days.yr=ifelse(tagDays>(365*2),tagDays-365*2,tag.days.yr))

#verified that we did collect data this year
lat.df%>%filter(substr(ts,1,4)=="2024"&recvProjID==280)%>%
  select(recvDeployName,reference,detectdoy)%>%unique()

#only check for 2024 detections
#lat.df<-lat.df%>%filter(substr(ts,1,4)=="2024")


print_maps=F

map_location="C:/Users/Steph/GitHub_data/survivalData/check_20250407/"

#print out a map of every included detection for each bird to manually check
if(print_maps==T){
  for(bird1 in unique(lat.df$name_in_vcf)){
    bird.lat<-lat.df%>%filter(name_in_vcf==bird1)%>%
      select(name_in_vcf,ts,tagDeployStart,recvDeployID,recvDeployLat,recvDeployLon,
             recvDeployName,detectyear,tagDays,detectdoy,tag.days.yr,doy.aug24,yr)%>%
      mutate(ts=substr(ts,1,10))%>%mutate(tagDeployStart=substr(tagDeployStart,1,10))%>%
      distinct()%>%arrange(tagDays)%>%
      filter(!is.na(tagDays))
    
    bird.title<-paste(bird1,"| towers:",length(unique(bird.lat$recvDeployName)),
                      "| days survived:",max(bird.lat$tagDays))
    vir.start<-min(bird.lat$doy.aug24)/365
    vir.end<-max(bird.lat$doy.aug24)/365
    
    png(paste(map_location,bird1,".png",sep=""),width=20,height=12,units="cm",res=300)
    print(ggplot(data = world) +
            geom_sf(fill="gray98")+
            geom_path(data=bird.lat,aes(x=recvDeployLon,y=recvDeployLat),size=0.8)+
            scale_colour_viridis(begin=vir.start,end=vir.end,
                                 name="day of year\nday one = Aug 24")+
            scale_shape_manual(values=unique(as.numeric(substr(bird.lat$yr,5,5))-1),
                               name="year")+
            geom_point(data=bird.lat,aes(x=recvDeployLon,y=recvDeployLat,
                                         colour=tag.days.yr,shape=yr),size=3,stroke=2)+
            theme_minimal()+ggtitle(bird.title)+
            coord_sf(xlim=c(-150,-40),ylim=c(2,60),expand=FALSE))
    
    dev.off() }}


