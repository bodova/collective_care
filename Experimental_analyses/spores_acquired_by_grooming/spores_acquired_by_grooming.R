# We here look at the spores found in the infrabuccal pockets of the ants (in their head) and relate it to their grooming activity by modelling the probability of finding spores in the head after having groomed an ant carrying those spores (Supp. Figure 4). 

## Note that original labels in the dataset 'POST_windowed_450_outPairwise' contain "H" for high fungal load and "L" for low fungal load, which are termed F (for H) and f (for L) in the publication. Control treatment (C) is labeled by Tx (Triton X). The replicates are labelled according to their original dish number (note that an overall of 9 dishes across groups could not be analysed due to technical errors, so that only 99 replicates could be included in the study). The file 'labeling.csv' shows replicate number alongside with the original dish number. 

#================================================================================#

rm(list = ls())
basedir='~/collective-care/Experimental_analyses/'
setwd(basedir)
source(paste0(basedir,'aux/load_with_scripts.R'))

library(glmmTMB)
library(DHARMa)
library(stringr)
library(tidyverse)

# Data: windowed tables where each row corresponds to the time series
#of each ant and one of the two treated ants they can groom (i.e. two rows per untreated ant)

tau_groom<-read.table ('./data/POST_windowed_450_outPairwise.csv',sep=",",header=F,na.strings = c("NA","None"))

partheader<-c('treat','rep','color','strain','dose','headRFP','headGFP','bodyRFP','bodyGFP','color2','strain2','dose2','headRFP2','headGFP2','bodyRFP2','bodyGFP2')
colnames(tau_groom)[1:16]<-c(partheader)

# where ‘treat’ refers to the treatment of the group 'HH','HL','HTx','LL','LTx','TxTx'; ‘rep’ to the dish number per treatment group; ’color’ the colorID of the performing individual; ‘strain’ the application of either the GFP(G) or the RFP(R) strain of the fungus, and N to no application for the nestmates, ‘dose’ high(H), low(L) or no (X) fungal dose applied; ‘headRFP’ and ‘headGFP’ to the RFP- resp. GFP-labelled spores quantified by PCR in the ant’s head; ‘bodyRFP’ and ‘bodyGFP’ to the RFP- resp. GFP-labelled spores quantified from the ant’s body; ‘color2’ the color of the receiving ant, as well as strain, dose, head- and body spores of the respective labels quantified on this second ant. 

totg=tau_groom[,17:196]%>%colSums()%>%sum()
tau_groom[,17:196]%>%colSums()%>%cumsum()%>%as_tibble()%>%tidyr::pivot_longer(cols = starts_with("V"))%>%mutate(name=1:180,value=value/totg)->temp
temp[temp$value<=0.34,]%>%tail() #33 percent grooming happens at window 36, within 18min
temp[temp$value<=0.67,]%>%tail() #66 percent grooming happens at window 90, within 45min

#exclude self grooming 
tau_groom<-tau_groom[tau_groom$color!=tau_groom$color2,]
#only nestmates actors
tau_groom<-tau_groom[tau_groom$strain=="N",]
#only spore-treated receivers
tau_groom<-tau_groom[tau_groom$strain2%in%c("G","R"),]

library(tidyverse)
totg=tau_groom[,17:196]%>%colSums()%>%sum()
tau_groom[,17:196]%>%colSums()%>%cumsum()%>%as_tibble()%>%tidyr::pivot_longer(cols = starts_with("V"))%>%mutate(name=1:180,value=value/totg)->temp
temp[temp$value<=0.33,]%>%tail() #33 percent grooming happens at window 32, within 16min
temp[temp$value<=0.67,]%>%tail() #66 percent grooming happens at window 82, within 41min

#partheader<-c('Treatment','Replicate',
#              'Color_P','Label_P','IndT_P','head_RFP_P','head_GFP_P','body_RFP_P','body_GFP_P',
#              'Color_R','Label_R','Dose_R','head_RFP_R','head_GFP_R','body_RFP_R','body_GFP_R')
#colnames(tau_groom)[1:16]<-c(partheader)
#
#write.table (tau_groom,'~/Documents/collcare_package/data/spores_acquired_by_grooming.csv')

#convert to long with tidyr::gather #last window should be num_windows+9
tau_groom<-tidyr::gather(tau_groom,window,groom, V17:V196)
#rename windows, and correct window number
tau_groom$window<-as.numeric(str_sub(tau_groom$window,start=2))-16

#spores in head
tau_groom$head<-ifelse(tau_groom$strain2=="G",tau_groom$headGFP,
                       ifelse(tau_groom$strain2=="R",tau_groom$headRFP,0))
#spores in body
tau_groom$body<-ifelse(tau_groom$strain2=="G",tau_groom$bodyGFP,
                       ifelse(tau_groom$strain2=="R",tau_groom$bodyRFP,0))

#select given time period
initial_w=1;final_w=180 #ok model assumptions ***
#initial_w=1;final_w=36 # ok model assumptions *** which corresponds roughly the time where half of the events happen
d<-tau_groom%>%subset(window%in%c(initial_w:final_w)) 

#sum grooming within interval and compare to head spores collected 
dh<-d%>%dplyr::group_by(treat,rep,color,color2,strain2)%>% #dyad
  #subset(head>0)%>%
  dplyr::summarise(sumgroom=sum(groom),head=max(head))
ants_tot<-dim(dh)[1]

print(paste0(dim(dh)[1]," ants out of ", ants_tot, " groomed within the interval ",initial_w,":", final_w ))

dh$ant<-as.factor(paste0(dh$treat,dh$rep,dh$color))
dh$dish<-as.factor(paste0(dh$treat,dh$rep))
dh$sumgroom_min<-dh$sumgroom/(15*60) 
dh$sqrthead<-sqrt(dh$head)

# Create the model: probability of finding each of the two spores given grooming duration towards rfp or gfp, separately 
dh$binhead<-factor(dh$sqrthead>0)
rfp<-subset(dh,strain2=="R")
full <- glmmTMB(binhead ~ sumgroom_min+(1|dish/ant),family=binomial, data = rfp)
null <- glmmTMB(binhead ~ 1+(1|dish/ant),family=binomial, data = rfp)
anova(null,full,test="Chisq")
simres <- simulateResiduals(full);plot(simres)
gfp<-subset(dh,strain2=="G")
full <- glmmTMB(binhead ~ sumgroom_min+(1|dish/ant),family=binomial, data = gfp)
null <- glmmTMB(binhead ~ 1+(1|dish/ant),family=binomial, data = gfp)
anova(null,full,test="Chisq")
simres <- simulateResiduals(full);plot(simres)

# Plot
# plot Suppl. Figure 4a
pA<-ggplot(dh,aes(y=sumgroom_min,x=binhead))+
  theme_pubr(base_size=14)+
  theme(legend.title = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.25),# ticks black
        axis.ticks.length=unit(-0.15, "cm"), #ticks inward
        axis.text.x= element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),#ticks inward
        axis.text.y= element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")))+ #ticks inward
  geom_violin(trim=T)+
  geom_boxplot(width=0.25,aes(color=strain2))+
  annotate("text", label=c("***"), x=1.5,y=40, size=8)+ 
  ylab("Total grooming time (min)")+xlab("Presence of spores in head")+
  scale_color_manual("Spores",values=c("R"="#ff5555","G"="#2BC20E"),labels = c("GFP","RFP"))+
  coord_flip()+
  scale_y_continuous(limits=c(0,50),breaks=seq(0,50,10))+
  scale_x_discrete(labels=c("FALSE" = "No", "TRUE" = "Yes"))

library(cowplot)
p=plot_grid(pA, labels = c("(A)"),ncol=1)
print(p)
summary(full)
#get odds
c<-fixef(full)
OR<-round(exp(c$cond[[2]]),3)
confint(full)[2,]
exp(confint(full)[2,])

library(ggeffects)
ggpredict(full)
p<-ggpredict(full, "sumgroom_min")
plot(p)

# plot regression between spores in their head (for those with >0) vs time grooming 
ggplot(dh[dh$head>0,],aes(x=sumgroom_min,y=head))+theme_pubr()+
  theme(legend.title = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.25),# ticks black
        axis.ticks.length=unit(-0.15, "cm"), #ticks inward
        axis.text.x= element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),#ticks inward
        axis.text.y= element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")))+ #ticks inward
  geom_point(aes(color=strain2),alpha=0.3,size=1)+
  geom_smooth(method="lm", color="black" )+
  xlab("Time grooming (min)")+ylab("Spores in head")+
  scale_color_manual(values=c("R"="red","G"="forestgreen"))

#---for such a nonlinear relationship we can report the strength of Spearman correlation
# for  increasingly large intervals, from the end  
corlist<-NULL 
intsize=180/9 
for (i in seq(0,180-intsize,intsize)){
  initial_w=i
  final_w=180
  d<-tau_groom%>%subset(window%in%c(initial_w:final_w))
  dh<-d%>%dplyr::group_by(treat,rep,color,color2,strain2)%>% #dyad
    subset(head>0)%>%
    dplyr::summarise(sumgroom=sum(groom),head=max(head))
  ants_tot<-dim(dh)[1]
  dh$sumgroom_min<-dh$sumgroom/(15*60) 
  nonzero<-dh[dh$head>0,]
  c<-cor.test(nonzero$head,nonzero$sumgroom_min,method="spearman", exact=F)
  corlist<-rbind(corlist,c(initial_w,final_w,c$statistic, c$estimate, c$p.value))
}
colnames(corlist)<-c("initial_w","final_w","statistic","rho","pvalue")
corlist<-as.data.frame(corlist)
signcols= c("FALSE"=19,"TRUE"=1)
library(ggpubr)
corlist$center<-(corlist$initial_w+corlist$final_w)/2


# plot Suppl. Figure 4b

pB<-ggplot(corlist, aes(initial_w/2,rho,shape=pvalue>=0.05))+
    scale_shape_manual("Significance",values=signcols,labels=c("p-value ≤ 0.05","non-significant"))+
    geom_point(size=2)+theme_pubr()+labs(y="Spearman rho of head spores and grooming time", x="Time after treatment (min)")+
    theme_pubr(base_size=14)+
    theme(legend.title = element_blank(),
          axis.ticks = element_line(colour = "black", size = 0.25),# ticks black
          axis.ticks.length=unit(-0.15, "cm"), #ticks inward
          axis.text.x= element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),#ticks inward
          axis.text.y= element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")))+ #ticks inward
    scale_x_continuous(breaks=seq(0,180,10))+ 
    scale_y_continuous(limits=c(0.09,0.25),breaks=c(0.1,0.15,0.2,0.25))
pB<-pB+geom_text(aes(x=initial_w/2+5,label = paste0(initial_w/2,"-",final_w/2))) 

#Suppl Figure 4
p=plot_grid(pA,pB, labels = c("A","B"),ncol=2,label_size = 16)
print(p)

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
#   [1] ggeffects_1.1.2    cowplot_1.1.1      DHARMa_0.4.5       glmmTMB_1.1.3      openxlsx_4.2.5     RColorBrewer_1.1-3
# [7] ggpubr_0.4.0       plyr_1.8.7         multcomp_1.4-19    TH.data_1.1-1      MASS_7.3-57        survival_3.3-1    
# [13] mvtnorm_1.1-3      lme4_1.1-29        Matrix_1.4-1       emmeans_1.7.5      forcats_0.5.1      stringr_1.4.0     
# [19] dplyr_1.0.9        purrr_0.3.4        readr_2.1.2        tidyr_1.2.0        tibble_3.1.7       ggplot2_3.3.6     
# [25] tidyverse_1.3.1   

