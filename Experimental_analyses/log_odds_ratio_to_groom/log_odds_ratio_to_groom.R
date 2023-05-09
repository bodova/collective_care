# Log-odds-ratio to groom as a function of delta (current spore load difference between treated-ants)

## In Fig. 6 (and Suppl. Fig. 9) we determine, how much more likely it is for a nestmate to groom the treated-ant with a current higher-load than the other-treated ant, given a difference in spore loads. In Supp. Fig. 7, we determine the log odds ratio based on a non-updated initial spore load. 

### Note that original labels in the data files “POST_windowed_450_outPairwise_EventStartCounts.csv” and “POST_windowed_450_groomTot_spore_MM_paper.csv” contain "H" for high fungal load and "L" for low fungal load, which are termed F (for H) and f (for L) in the publication. The control group (C in the publication) is referred to as "Tx" (Triton X). Replicates are labelled according to their original dish number (note that an overall of 9 dishes across groups could not be analysed due to technical errors, so that only 99 replicates could be included in the study). The file 'labeling.csv' shows replicate number alongside with the original dish number. 

#### Note that only groups with spores were included in the analysis (HH, HL, LL, HTx, LTx). 

# TABLE
# delta_bin_interval, delta_bin_id, delta_bin_center,
# treatment, group (used in Suppl. Fig.9 )
# num_choose_higher:number of times that the ant with higher current load was groomed 
# num_choose_lower :number of times that the ant with lower current load was groomed 
# lor_groom:  log{(num_choose_higher)/(num_choose_lower)}

#------------------#
rm(list = ls())
basedir='~/collective-care/Experimental_analyses/'
setwd(basedir)
source(paste0(basedir,'aux/load_with_scripts.R'))

#LIBRARIES
library(tidyverse)
library(data.table)
library(stringr)

#COLOR PALETTE
treatfill<-c('LL'="#EDB116",'HH'= "#D95319",  'HL'="#EDB116", 'LTx'="#999999",'HTx'= "#999999") #
treatcols<-c('LL'="#EDB116",'HH'= "#D95319",  'HL'="#D95319", 'LTx'="#EDB116",'HTx'= "#D95319") #

#READ DATA
tau_groom<-read.table ('./data/POST_windowed_450_outPairwise_EventStartCounts.csv',sep=",",header=F)
#rows: ants, columns: 2 dish+ 3 performing ant identifier+4 spore data variables + 3 receiving ant identifier+ 4 spore data variables 
# time series  of events of  grooming (performed) every 30s -450frames

spore_sp<-read.table ('./data/POST_windowed_450_groomTot_spore_MM_paper.csv',sep=",",header=F) 
#rows:ants, columns:2 dish+ 3 ant identifier variables + 4 spore data variables 
# time series of spore load estimated every 30s -450frames (see Modelling_simulations/MM_backcompute in this repository)

#Add headers for 9 column variables + 180 windows
partheader<-c('treat','rep','color','strain','dose','headRFP','headGFP','bodyRFP','bodyGFP','color_R','strain_R','dose_R','headRFP_R','headGFP_R','bodyRFP_R','bodyGFP_R') #_R for recipient

# where ‘treat’ refers to the treatment of the group 'HH','HL','HTx','LL','LTx'; ‘rep’ to the dish per treatment group; ’color’ the colorID of the performing individual; ‘strain’ the application of either the GFP(G) or the RFP(R) strain of the fungus, and N to no application for the nestmates, ‘dose’ high(H), low(L) or no (X) fungal dose applied; ‘headRFP’ and ‘headGFP’ to the RFP- resp. GFP-labelled spores quantified by PCR in the ant’s head; ‘bodyRFP’ and ‘bodyGFP’ to the RFP- resp. GFP-labelled spores quantified from the ant’s body; ‘colorR’ the color of the receiving ant, as well as strain, dose, head- and body spores of the respective labels quantified on this receiving ant. The first column with 'ReplicateNr' gives the replicate number.

colnames(tau_groom)[1:16]<-c(partheader);
partheader<-c('treat','rep','color_R','strain','dose','headRFP','headGFP','bodyRFP','bodyGFP')
colnames(spore_sp)[1:9]<-c(partheader)

#exclude self behaviors (treated to treated)
tau_groom<-subset(tau_groom,tau_groom$color!=tau_groom$color_R)
#exclude TxTx
tau_groom<-subset(tau_groom,treat!="TxTx");spore_sp<-subset(spore_sp,treat!="TxTx")

#grooming from N towards treated
tau_groom<-subset(tau_groom,strain=="N")# from nestmates
tau_groom<-subset(tau_groom,strain_R%in%c("R","G","T"))# towards treated
spore_sp<-subset(spore_sp,strain!="N")
#check tau_groom (n=656 rows) spore_sp(n=164 rows) 

#convert to long with tidyr::gather #last window should be num_windows+9
tau_groom<-tidyr::gather(tau_groom,window,groom, V17:V196)
spore_sp<-tidyr::gather(spore_sp,window,spores, V10:V189)
#rename windows, and correct window number
tau_groom$window<-as.numeric(str_sub(tau_groom$window,start=2))-16
spore_sp$window<-as.numeric(str_sub(spore_sp$window,start=2))-9 
#for Tx treated there should be no initial spore counts, NA-->0
spore_sp[spore_sp$dose=='X',]$spores=0
treatlist<-c('HTx','LTx','HL','HH','LL')

# given that the current (instant) spore loads estimated at the start of 
# a window can be largely influenced by grooming activity happening in
# the same window (i.e. since we backcompute from final loads) we examine 
# the grooming choices in later windows (shifted by wshift) 

cortable<-NULL #values for correlation at different shifts
#for (wshift in seq(10,1,-1)) {

wshift=1
tdata<-tau_groom%>%filter(treat%in%(treatlist))
sdata<-spore_sp%>%filter(treat%in%(treatlist))

tdata<-subset(tdata,window<=180-wshift)
sdata<-subset(sdata,window>wshift)
sdata$window<-sdata$window-wshift

#assign rank according to spores rankinstant: 1 = higher spore load 2 = lower spore load 
sdata<-sdata%>%group_by(treat,rep,window)%>% dplyr::mutate(rankinstant=rank(-spores))
sdata$rankinstant<-as.factor(sdata$rankinstant)
sdata<-as.data.frame(sdata)
sdata<-sdata%>%dplyr::select(treat,rep,window,color_R, spores,rankinstant)
xx<-merge(tdata,sdata)

x<-xx%>%dplyr::select(treat,rep,color,window,groom,spores,rankinstant) 
rankinstant1<-subset(x,rankinstant==1)
rankinstant2<-subset(x,rankinstant==2)
colnames(rankinstant2)<-c('treat','rep','color','window','groom2','spores2','rankinstant2')
ri<-merge(rankinstant1,rankinstant2)

data<-ri
data$sp_diff<-data$spores-data$spores2
data<-data[!(data$groom==0 & data$groom2==0),] #exclude windows of no decisions towards treated 

data<-data[ order( data[,1], data[,2], data[,3],data[,4] ),]
data$groom_diff<-data$groom-data$groom2 

data$groom_high<-rep(NA);data$groom_low<-rep(NA)
data$groom_high[which(data$groom_diff>0)] <- data$groom[which(data$groom_diff>0)]
data$groom_low [which(data$groom_diff>0)] <- data$groom2[which(data$groom_diff>0)]
data$groom_high[which(data$groom_diff<=0)] <- data$groom[which(data$groom_diff<=0)]
data$groom_low [which(data$groom_diff<=0)] <- data$groom2[which(data$groom_diff<=0)] 
 
data$dish<-paste0(data$treat,data$rep)
#write.csv(data,paste0("~/Documents/collcare_package/data/per_window_N2treated_and_delta.csv"))
#--------------------------------------------------------------------------------------------------------------
slices=10 # slice data fairly into spore load difference (transformed) categories
data$bins<-rep(0);data$bins_interval<-rep(0)
for (tr in treatlist){ #per-treatment slicing
  v<-round(c(unname(unique(quantile(data[data$treat==tr,]$sp_diff,seq(0,1,1/slices))))),10) 
  v[slices+1]<-v[slices+1]+.1 #to make sure to include all data, upper bound
  v[1]<-v[1]-.1 #to make sure to include all data, lower bound
  #quantile, choose breaks roughly of equal content
  data[data$treat==tr ,]$bins<-as.factor(cut(data[data$treat==tr,]$sp_diff,breaks=v, include.lowest = T))
  #bins from 1 to slices
  data[data$treat==tr ,]$bins_interval<-as.character(cut(data[data$treat==tr ,]$sp_diff,breaks=v,include.lowest = T))
}

data$bins_int<-gsub("[()]","",data$bins_interval) #remove ()
data$bins_int<-gsub("\\[|\\]","",data$bins_int) #remove []
data<-data %>% tidyr::separate(bins_int, c("x1", "x2"),sep=",")
data$bins_x<-(as.numeric(as.character(data$x1))+as.numeric(as.character(data$x2)))/2

#write.csv(data,paste0("~/Documents/collcare_package/data/log_odds_ratio_vs_spore_load_diff.csv"))

table(data$bins,data$treat)

TABLE<-data%>%dplyr::select(treat, bins,bins_x,bins_interval,groom_high,groom_low)%>%group_by(treat,bins_interval,bins,bins_x)%>% dplyr::summarise(num_choose_higher= sum(groom_high),num_choose_lower= sum(groom_low))%>%ungroup()

TABLE$lor_groom<-log((TABLE$num_choose_higher)/(TABLE$num_choose_lower))

TABLE<-TABLE %>% mutate(group = fct_collapse(treat,
                                 two_treated = c("HH", "HL", "LL"),
                                 one_treated = c("HTx", "LTx")))
# looks at groups with 2 or only one spore-treated ants separately 

cor<-cor.test(TABLE$lor_groom,TABLE$bins_x,method="spearman")
cortable<-rbind(cortable,c(cor$estimate,cor$p.value,wshift))
cortable<-as.data.frame(cortable)
library(ggpubr)
A<-ggplot(TABLE, aes(x= bins_x,y = lor_groom))+
  theme_pubr()+
  theme(legend.title = element_blank())+
  scale_color_manual(values=treatcols)+
  scale_fill_manual(values=treatfill)+
  geom_hline(yintercept=0,color="gray")+
  scale_x_continuous(limits=c(0,3.6e5))+scale_y_continuous(limits=c(-2.1,2.1))+
  geom_point(size=1.8, shape=21,stroke=2,aes(color=treat,fill=treat), position=position_dodge(width=0.01))
 
A<-A+annotate("text",x = c(100000,100000), y = c(-2,2), label = c("Groom lower spore load ", "Groom higher spore load"),size=5)+ 
    annotate("text",x = c(250000), y = c(-2), label =paste0("cor=", c(round(cor$estimate,2))),size=5)+
 labs( y= paste0("Log odds ratio "), x= "Spore load difference")


signcols= c("TRUE"=19,"FALSE"=1)

plot(A)


#write.csv(TABLE,paste0("~/Documents/collcare_package/data/log_odds_ratio_vs_spore_load_diff_TAB.csv")) 

#--------------------------------------------------------------------------------------------
cor<-cor.test(TABLE[TABLE$group=="one_treated",]$lor_groom,TABLE[TABLE$group=="one_treated",]$bins_x,method="spearman")
Suppl_Fig_9a<-ggplot(subset(TABLE,TABLE$group=="one_treated"), aes(x= bins_x,y = lor_groom))+
  theme_pubr()+
  theme(legend.title = element_blank())+
  scale_color_manual(values=treatcols)+
  scale_fill_manual(values=treatfill)+
  geom_hline(yintercept=0,color="gray")+
  scale_x_continuous(limits=c(0,3.6e5))+scale_y_continuous(limits=c(-2.1,2.1))+
  geom_point(size=1.8, shape=21,stroke=2,aes(color=treat,fill=treat), position=position_dodge(width=0.01))+
  annotate("text",x = c(100000,100000), y = c(-2,2), label = c("Groom lower spore load ", "Groom higher spore load"),size=5)+ 
  annotate("text",x = c(250000), y = c(-2), label =paste0("cor=", c(round(cor$estimate,2))),size=5)+
  labs( y= paste0("Log odds ratio "), x= "Spore load difference")
plot(Suppl_Fig_9a)
cor<-cor.test(TABLE[TABLE$group=="two_treated",]$lor_groom,TABLE[TABLE$group=="two_treated",]$bins_x,method="spearman",exact=F)
Suppl_Fig_9b<-ggplot(subset(TABLE,TABLE$group=="two_treated"), aes(x= bins_x,y = lor_groom))+
  theme_pubr()+
  theme(legend.title = element_blank())+
  scale_color_manual(values=treatcols)+
  scale_fill_manual(values=treatfill)+
  geom_hline(yintercept=0,color="gray")+
  scale_x_continuous(limits=c(0,3.6e5))+scale_y_continuous(limits=c(-2.1,2.1))+
  geom_point(size=1.8, shape=21,stroke=2,aes(color=treat,fill=treat), position=position_dodge(width=0.01))+
  annotate("text",x = c(100000,100000), y = c(-2,2), label = c("Groom lower spore load ", "Groom higher spore load"),size=5)+ 
  annotate("text",x = c(250000), y = c(-2), label =paste0("cor=", c(round(cor$estimate,2))),size=5)+
  labs( y= paste0("Log odds ratio "), x= "Spore load difference")
plot(Suppl_Fig_9b)

#--------------------------------------------------------------------------------------------
# divide time so that each part contains roughly the same number of grooming events
# log odds ratio to groom based on APPLIED INITIAL SPORE DOSE
ALLTABLES<-NULL
treatlist<-c('HTx','LTx','HL')
#split observation period into three parts containing roughly a third of the observed events each 
partitions<-list(c(1,36),c(37,90),c(91,180)) # equivalent to minutes 0-18, 18-45, and 45-90 of the experiment (see Suppl. Fig. 7)

wshift=1
for (p in 1:3){
  firstw<-partitions[[p]][1]
  lastw<-partitions[[p]][2]
  tdata<-tau_groom%>%filter(treat%in%(treatlist))
  tdata<-subset(tdata,window<=180-wshift)
  tdata<-tdata[tdata$window%in%c(firstw:lastw),]
  
  #assign rank based on applied dose
  tdata$rank=case_when(
    tdata$treat=="LTx" & tdata$dose_R=="L"~1,
    tdata$treat=="HTx" & tdata$dose_R=="H"~1,
    tdata$treat=="HL"  & tdata$dose_R=="H"~1,
    TRUE~2)
  
  x<-tdata%>%dplyr::select(treat,rep,color,window,groom,rank) 
  rank1<-subset(x,rank==1)
  rank2<-subset(x,rank==2)
  colnames(rank2)<-c('treat','rep','color','window','groom2','rank2')
  ri<-merge(rank1,rank2)
  
  data<-ri
  data<-data[!(data$groom==0 & data$groom2==0),] 
  data<-data[ order( data[,1], data[,2], data[,3],data[,4] ),]
  data$groom_diff<-data$groom-data$groom2 
  
  data$groom_high<-rep(NA);data$groom_low<-rep(NA)
  data$groom_high[which(data$groom_diff>0)] <- data$groom[which(data$groom_diff>0)]
  data$groom_low [which(data$groom_diff>0)] <- data$groom2[which(data$groom_diff>0)]
  data$groom_high[which(data$groom_diff<=0)] <- data$groom[which(data$groom_diff<=0)]
  data$groom_low [which(data$groom_diff<=0)] <- data$groom2[which(data$groom_diff<=0)] 
  
  data$dish<-paste0(data$treat,data$rep)
  
  TABLE<-data%>%dplyr::select(treat,dish,groom_high,groom_low)%>%group_by(treat,dish)%>%
    dplyr::summarise(num_choose_higher= sum(groom_high),num_choose_lower= sum(groom_low))%>%
                       ungroup() #per dish #delted n = n() inside summarise
  
  #Haldanes correction, when one category is zero you get LOR=0 or INF, add one to these cells
  TABLE[TABLE$num_choose_higher==0|TABLE$num_choose_lower==0,]$num_choose_higher<-TABLE[TABLE$num_choose_higher==0|TABLE$num_choose_lower==0,]$num_choose_higher+1
  TABLE[TABLE$num_choose_higher==0|TABLE$num_choose_lower==0,]$num_choose_lower<-TABLE[TABLE$num_choose_higher==0|TABLE$num_choose_lower==0,]$num_choose_lower+1
  
  TABLE$lor_groom<-log((TABLE$num_choose_higher)/(TABLE$num_choose_lower))
  TABLE<-TABLE%>%group_by(treat)#%>% dplyr::mutate(se=sd(lor_groom,na.rm=T)/sqrt(n))
  TABLE$part<-rep(p)
  ALLTABLES<-rbind(ALLTABLES,TABLE)
}

#write.csv(ALLTABLES,paste0("~/Documents/collcare_package/data/log_odds_ratio_in_time.csv"))

TABLE_mean<-ALLTABLES%>%dplyr::select(lor_groom,treat,part)%>%group_by(treat,part)%>%
  dplyr::summarise(y=mean(lor_groom,na.rm=T),sd=sd(lor_groom,na.rm=T),n=n(),se=sd/sqrt(n)) 
TABLE_mean$lowerse<-TABLE_mean$y-TABLE_mean$se 
TABLE_mean$upperse<-TABLE_mean$y+TABLE_mean$se

TABLE_mean$treat<-factor(TABLE_mean$treat,levels=c("HL","HTx","LTx"))

TABLE_mean$x<-c(1,2,3,5,6,7,9,10,11)
#write.csv(TABLE_mean,paste0("~/Documents/collcare_package/data/log_odds_ratio_in_time_TAB.csv")) 

P<-ggplot(TABLE_mean, aes(x=x,y = y,ymin=lowerse,ymax=upperse))+
  geom_vline(xintercept=c(4,8),color="gray90",size=.3)+
  geom_hline(yintercept=0,color="gray50",linetype="FFFF",size=.2)+
  geom_errorbar(width=0.2)+
  geom_point(size=2, shape=21,stroke=1.8,aes(color=treat,fill=treat))+
  theme_pubr(base_size=14)+theme(legend.position="none")+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, 'cm'),legend.spacing.x = unit(.7, 'cm'),
        axis.ticks = element_line(colour = "black", size = .2),# ticks black
        axis.ticks.length=unit(-0.15, "cm"), #ticks inward
        axis.text.x= element_text(size=11,margin=unit(c(0.5,0.5,0.5,0.5), "cm")),#ticks inward
        axis.text.y= element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")))+ #ticks inward
  scale_color_manual(values=treatcols)+
  scale_fill_manual(values=treatfill)+
  scale_y_continuous(limits=c(-2,2),breaks=c(-2,-1,0,1,2))+
  scale_x_continuous(breaks=c(1,2,3,5,6,7,9,10,11),
                     labels=c("0-18","18-45","45-90","0-18",
                              "18-45","45-90","0-18","18-45","45-90"))

Suppl_Fig_7<-P+annotate("text",x = c(3.8,3.8), y = c(-1.9,1.9), 
              label = c("Groom lower applied spore load ",
               "Groom higher applied spore load"),size=5)+ 
  labs(y= "Log odds ratio", x= "Time after treatment (min)")

print(Suppl_Fig_7)

# 
# > sessionInfo()
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux 10 (buster)
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
# LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.3.5.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggpubr_0.4.0      data.table_1.14.2 forcats_0.5.1     dplyr_1.0.9       purrr_0.3.4       readr_2.1.2      
# [7] tidyr_1.2.0       tibble_3.1.7      ggplot2_3.3.6     tidyverse_1.3.1   stringr_1.4.0    


