################################################################################
# extract migratory phenotypes from filtered motus data
################################################################################


#load packages
pkgs<-c("tidyverse","viridis","ggpubr","geosphere","sf","rnaturalearth")
sapply(pkgs,library,character.only=T); rm(pkgs)
theme_set(theme_classic())

#filter motus detections - output from survival_cleaning.R
lat.df<-read_csv("C:/Users/Steph/GitHub_data/survivalData/motusDetections.filtered.20250409.csv")

#metadata - estimated in collect_AIMs_metadata.R
meta_thrush<-read_csv("C:/Users/Steph/GitHub_data/AIMs_metaData/AIMs_metadata.20250327.csv")

#survival data - estimated in survival_cleaning.R
#only used to format for Sarah's analysis
# thrush_surv<-read.csv("C:/Users/Steph/GitHub/thrush_hybrids/updateSurvival/motusSurvivalAll.20250327.csv")


#edit to quantNorm function - remove NAs from estimation
quantNorm =function(x){
  x1<-x[!is.na(x)]
  x[!is.na(x)]<-qnorm(rank(x1,ties.method = "average")/(length(x1)+1))
  x
}

#check that quantnorm output makes sense
v1<-c(10,9,NA,1,2,NA,-7,-10)
quantNorm(v1)

meta_thrush<-meta_thrush%>%
  
  #only keep juveniles that were radio tagged
  filter(age_release=="HY"&tag_type=="radio")%>%
  
  select(reference,name_in_vcf,age_release,tag_type,
         release_site,release_gps.n,release_gps.w,
         capture_site,capture_gps.n,capture_gps.w,
         sex_binary,
         release_year,release_date)%>%
  
  #get day of year for bird release date
  mutate(releaseDay=as.POSIXlt(release_date, format = "%Y-%m-%d")$yday)

#define a projection - sticking with a commonly used one here
projection<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#keep detections from first year tagged, Aug 1 onwards
#summarise as detection events - one detection per motus station per bird per day
lat.df.fall1<-lat.df%>%filter(detectyear==release_year)%>%
  filter(detectdoy>212)%>%
  select(name_in_vcf,motustagid,tsDay,tagDeployStart,tagDepLat,tagDepLon,
         recvDeployID,recvDeployLat,recvDeployLon,recvDeployName,
         detectyear,tagDays,detectdoy,deployDay,
         release_site,capture_site,release_year)%>%
  mutate(tagDeployStart=substr(tagDeployStart,1,10))%>%
  distinct%>%arrange(tagDays)%>%
  filter(!is.na(tagDays))



#get mean coords for birds tagged in Pemberton
#falls near Shadow Lake - should be correct
pemberton<-lat.df.fall1%>%filter(release_site=="Pemberton")%>%
  select(name_in_vcf,tagDepLat,tagDepLon)%>%
  distinct()%>%summarise(lon=mean(tagDepLon),lat=mean(tagDepLat))

#estimate the distance between release site and detection site
#estimation needs to be pairwise, otherwise function will generate a matrix of all possible combos
lat.df.fall1$releaseSiteDistance<-NA
lat.df.fall1$pembertonDistance<-NA
for(i in 1:nrow(lat.df.fall1)){
  lat.df.fall1$releaseSiteDistance[i]<-distm(lat.df.fall1[i,c("tagDepLon","tagDepLat")],
                                     lat.df.fall1[i,c("recvDeployLon","recvDeployLat")],
                                     fun = distHaversine)
  lat.df.fall1$pembertonDistance[i]<-distm(pemberton,
                                        lat.df.fall1[i,c("recvDeployLon","recvDeployLat")],
                                        fun = distHaversine)
  }



# #retain points within 300km of the release site


#determine 'tiers' of detections based on distance from deployment site or latitude band
lat.df.fall1.plotTiers<-lat.df.fall1%>%
  mutate(distTier=case_when(pembertonDistance<300000~"tier1",
                            pembertonDistance>300000&pembertonDistance<1500000~"tier2",
                            pembertonDistance>1500000&pembertonDistance<4600000~"tier3",
                            pembertonDistance>4600000~"tier4"))%>%
  mutate(distTier2=case_when(pembertonDistance<300000~"tier1",
                            recvDeployLat>45~"tier2",
                            recvDeployLat>35~"tier3",
                            recvDeployLat>25~"tier4",
                            recvDeployLat<25~"tier5"))%>%
  arrange(-tagDays)


world <- ne_countries(scale = "medium", returnclass = "sf")
ggpubr::ggarrange(ggplot(data = world) +
                    geom_sf(fill="gray98")+
                    geom_point(data=lat.df.fall1.plotTiers,size=2,
                               aes(x=recvDeployLon,y=recvDeployLat,colour=distTier))+
                    coord_sf(xlim=c(-150,-40),ylim=c(0,60),expand=FALSE)+
                    scale_colour_viridis(discrete=T,end=0.9)+
                    ggtitle("distance"),
                  ggplot(data = world) +
                    geom_sf(fill="gray98")+
                    geom_point(data=lat.df.fall1.plotTiers,size=2,
                               aes(x=recvDeployLon,y=recvDeployLat,colour=distTier2))+
                    coord_sf(xlim=c(-150,-40),ylim=c(0,60),expand=FALSE)+
                    scale_colour_viridis(discrete=T,end=0.9)+
                    ggtitle("latitude"))

#Work with distTier - makes more sense for estimating bearings
#estimate as distance from Pemberton and as distance from release site
#difference matters for translocated birds and birds released in Tatlayoko, Bella Coola etc


lat.df.fall1<-lat.df.fall1%>%
  mutate(pembertonDistTier=case_when(pembertonDistance<300000~"tier1",
                                     pembertonDistance>300000&pembertonDistance<1500000~"tier2",
                                     pembertonDistance>1500000&pembertonDistance<4600000~"tier3",
                                     pembertonDistance>4600000~"tier4"))%>%
  
  mutate(releaseSiteDistTier=case_when(releaseSiteDistance<300000~"tier1",
                                       releaseSiteDistance>300000&releaseSiteDistance<1500000~"tier2",
                                       releaseSiteDistance>1500000&releaseSiteDistance<4600000~"tier3",
                                       releaseSiteDistance>4600000~"tier4"))%>%
  arrange(-tagDays)

ggplot(data = world) +
  geom_sf(fill="gray98")+
  geom_point(data=lat.df.fall1,size=2,
             aes(x=recvDeployLon,y=recvDeployLat,colour=pembertonDistTier))+
  coord_sf(xlim=c(-150,-40),ylim=c(0,60),expand=FALSE)+
  scale_colour_viridis(discrete=T,end=0.9)

ggplot(data = world) +
  geom_sf(fill="gray98")+
  geom_point(data=lat.df.fall1,size=2,
             aes(x=recvDeployLon,y=recvDeployLat,colour=releaseSiteDistTier))+
  coord_sf(xlim=c(-150,-40),ylim=c(0,60),expand=FALSE)+
  scale_colour_viridis(discrete=T,end=0.9)



################################################################################
#estimate bearing and day at each time point
################################################################################

for(tier in sort(unique(lat.df.fall1$pembertonDistTier))){
  
  tier.df<-lat.df.fall1%>%
    filter(pembertonDistTier==tier)%>%
    #keep the first record within the tier for each bird
    arrange(tagDays)%>%filter(!duplicated(motustagid))
  
  day1<-paste("fall_pemberton_detectDay",tier,sep="_")
  tier.df[day1]<-tier.df["detectdoy"]
  
  bear1<-paste("fall_pemberton_bearing",tier,sep="_")
  tier.df[bear1]<-NA
  
  for(i in 1:nrow(tier.df)){

    tier.df[i,bear1]<-
      bearingRhumb(p1=cbind(pemberton$lon,pemberton$lat),
                   p2=cbind(tier.df$recvDeployLon[i],
                            tier.df$recvDeployLat[i]))
  }
  
  meta_thrush<-meta_thrush%>%
    left_join(tier.df%>%select(name_in_vcf,all_of(c(day1,bear1))))
  
  rm(tier.df,day1,bear1)
  
}

for(tier in sort(unique(lat.df.fall1$releaseSiteDistTier))){
  
  tier.df<-lat.df.fall1%>%
    filter(releaseSiteDistTier==tier)%>%
    #keep the first record within the tier for each bird
    arrange(tagDays)%>%filter(!duplicated(name_in_vcf))
  
  day1<-paste("fall_releaseSite_detectDay",tier,sep="_")
  tier.df[day1]<-tier.df["detectdoy"]
  
  bear1<-paste("fall_releaseSite_bearing",tier,sep="_")
  tier.df[bear1]<-NA
  
  for(i in 1:nrow(tier.df)){
    
    tier.df[i,bear1]<-
      bearingRhumb(p1=cbind(tier.df$tagDepLon[i],
                            tier.df$tagDepLat[i]),
                   p2=cbind(tier.df$recvDeployLon[i],
                            tier.df$recvDeployLat[i]))
  }
  
  meta_thrush<-meta_thrush%>%
    left_join(tier.df%>%select(name_in_vcf,all_of(c(day1,bear1))))
  
  rm(tier.df,day1,bear1)
  
}

################################################################################
#calculate winter longitude / bearing
################################################################################

meta_thrush$fall_lon_at_lat35<-NA
for(bird1 in unique(lat.df.fall1$name_in_vcf)){
  
  #select focal bird
  bird.fall1<-lat.df.fall1%>%filter(name_in_vcf==bird1)
  
  #filter for points above 35 degrees latitude and below 35 degrees
  above35<-bird.fall1%>%filter(recvDeployLat>35)
  below35<-bird.fall1%>%filter(recvDeployLat<35)
  if(nrow(above35)>0&nrow(below35)>0){
    
    #keep detections closest to 35
    recv1=above35%>%arrange(recvDeployLat)%>%slice(n())
    recv2=below35%>%arrange(recvDeployLat)%>%slice(1)
    lon1<-gcLon(recv1[,c("recvDeployLon","recvDeployLat")],
          recv2[,c("recvDeployLon","recvDeployLat")],
          35)
    meta_thrush$fall_lon_at_lat35[meta_thrush$name_in_vcf==bird1]<-lon1[lon1<0]
    
  }
}



################################################################################
#calculate winter longitude / bearing
################################################################################


#keep detections from first winter tagged
#summarise as detection events - one detection per motus station per bird per day
#first winter - Nov 1 onward in tagged year or before March 1 in next year
lat.df.winter1<-lat.df%>%filter((detectyear==release_year&detectdoy>305)|
                  (detectyear==(release_year+1))&detectdoy<60)%>%
  
  select(name_in_vcf,motustagid,tsDay,tagDeployStart,tagDepLat,tagDepLon,
         recvDeployID,recvDeployLat,recvDeployLon,recvDeployName,
         detectyear,tagDays,detectdoy,deployDay,
         release_site,capture_site,release_year)%>%
  mutate(tagDeployStart=substr(tagDeployStart,1,10))%>%
  distinct()%>%
  filter(!is.na(tagDays))

#estimate the distance between release site and detection site
#estimation needs to be pairwise, otherwise function will generate a matrix of all possible combos
lat.df.winter1$pembertonDistance<-NA
for(i in 1:nrow(lat.df.winter1)){
  lat.df.winter1$pembertonDistance[i]<-distm(pemberton,
                                        lat.df.winter1[i,c("recvDeployLon","recvDeployLat")],
                                        fun = distHaversine)
}
  
  
lat.df.winter1<-lat.df.winter1%>%
  
  #retain points past 4000km of the release site
  #roughly relates to Florida/southern Mexico and further
  filter(pembertonDistance>4000000)%>%
  
  #keep the longest distance associated with each bird
  arrange(-pembertonDistance)%>%filter(!duplicated(name_in_vcf))%>%
  
  #get receiver location and normalize
  mutate(winter_longitude=recvDeployLon,winter_latitude=recvDeployLat)



ggplot(data = world) +
  geom_sf(fill="gray98")+
  geom_point(data=lat.df.winter1,size=2,colour='orangered3',
             aes(x=recvDeployLon,y=recvDeployLat))+
  coord_sf(xlim=c(-150,-40),ylim=c(0,60),expand=FALSE)

meta_thrush<-meta_thrush%>%
  left_join(lat.df.winter1%>%select(name_in_vcf,winter_longitude,winter_latitude))



#get the proportion of individuals with data for each trait

colSums(!is.na(meta_thrush%>%select(starts_with("fall_")|starts_with("winter_"))))/
  nrow(meta_thrush)

#limit to birds caught & released in Pemberton / Whistler

#proportion
colSums(!is.na(meta_thrush%>%filter(release_site%in%c("Pemberton","Whistler")&
                                      capture_site%in%c("Pemberton","Whistler"))%>%
                 select(starts_with("fall_")|starts_with("winter_"))))/
  nrow(meta_thrush%>%filter(release_site%in%c("Pemberton","Whistler")&
                              capture_site%in%c("Pemberton","Whistler")))

#total count
colSums(!is.na(meta_thrush%>%filter(release_site%in%c("Pemberton","Whistler")&
                                      capture_site%in%c("Pemberton","Whistler"))%>%
                 select(starts_with("fall_")|starts_with("winter_"))))


today=gsub("-","",Sys.Date())
outFile=paste("C:/Users/Steph/GitHub_data/survivalData/motusMigratoryPhenotypes",today,"csv",sep=".")
write.csv(meta_thrush,outFile)

