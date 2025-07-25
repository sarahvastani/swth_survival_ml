library(tidyverse)
#library(mgcv)
theme_set(theme_classic())


library(marked)

dataset1<-read_csv("C:/Users/Steph/GitHub/thrush_hybrids/migratory_traits/Dataset1.csv")

ch1<-read_csv("C:/Users/Steph/GitHub/thrush_hybrids/updateSurvival/thrush_survival_ch.20250707.csv")

dataset1<-dataset1%>%
  left_join(ch1%>%select(reference,ch_10days))%>%
  mutate(ch=gsub("days_","",ch_10days))

model1="hmmCJS"

model1.parameters=list(Phi=list(formula=~time),
                       p=list(formula=~time))
dataset1$phi<-NA
for(i in 1:nrow(dataset1)){
  if(!is.na(dataset1[i,'ch'])){
    #use try() so that for loop doesn't fail and exit if parameters are wrong
    cjsfit1=crm(dataset1[i,],model=model1,
                model.parameters=model1.parameters)
    dataset1$phi[i]<-prod(cjsfit1$results$reals$Phi$estimate)
    rm(cjsfit1)
    print(i)
  }
}

ggplot(dataset1,
       aes(x=as.character(ch),y=phi))+
  geom_point()
dataset1<-dataset1%>%
  mutate(phi_binary=if_else(phi>0.5,1,0))
dataset1%>%filter(phi>0.9)%>%pull(phi)%>%unique()

cor(dataset1$t5_spring40,dataset1$phi,method='pearson')

print(
dataset1%>%filter(phi_binary!=t5_spring40)%>%
  select(reference,ch,t5_spring40,phi),
n=50)
dataset1%>%filter(is.na(phi))
dataset1%>%filter(phi>0.5&phi<0.9)%>%pull(ch)

#write_csv(dataset1,"C:/Users/Steph/GitHub/thrush_hybrids/migratory_traits/Dataset1.phiSurvival.20250707.csv")
#write_csv(dataset1%>%
#            mutate(ch=paste('days_',ch,sep=''))%>%
#            select(reference,ch,phi_binary),
#          "C:/Users/Steph/GitHub/thrush_hybrids/updateSurvival/thrush_survival_phiBinary.10days.20250708.csv")


dataset2<-read_csv("C:/Users/Steph/GitHub/thrush_hybrids/migratory_traits/Dataset1.phiSurvival.20250707.csv")
dataset2$phi
dataset2%>%filter(phi>0.5&phi<0.9)%>%pull(phi)
dataset2%>%filter(phi<0.5&phi>0.1)%>%pull(phi)
dataset2%>%filter(phi>0.9)%>%pull(phi)
dataset2%>%filter(phi<0.1)%>%pull(phi)%>%max()


dataset1$t5_spring40
rm(cjsfit2)
