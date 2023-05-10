# Hygiene behavior of the ants (Figure 1 and Supplementary Figure 2)
# Plots and statistical analysis of the hygiene behavior of the ants, including selfgrooming (Fig. 1b), performed grooming (Fig. 1c), received grooming (Fig. 1d) and poison uptake behavior (Suppl. Fig. 2) 

## Note that original labels in the dataset 'behavior_aggregated' contain "H" for high fungal load and "L" for low fungal load, which are termed F (for H) and f (for L) in the publication. Control treatment (C) is labeled by "C". 
## Also, the replicates are numbered according to their original dish number, but see 'labeling.csv' containing the original alongside the updated replicate labels given in the  where 99 replicates were included in the study. 

### Note the following behavior abbreviations used in the raw data files: selfgrooming: "selfG", grooming_performed: "groomOut", grooming_received: groom_In". Poison uptake behavior from the ants abdominal (gaster) tip is labelled "gaster". Also the difference between post-treatment and pre-treatment values, is abbreviated with two letters  SD: selfgrooming delta : "SD", delta in grooming performed: "PD", delta in grooming received: "RD", and delta in poison uptake behavior from the gaster: "GD".


#================================================================================#

rm(list = ls())
basedir='~/collective_care/Experimental_analyses/'
setwd(basedir)
source(paste0(basedir,'aux/load_with_scripts.R'))

# load data
xdata<-read.table ("./data/behaviour_aggregated.csv",sep=",",header=TRUE)
xdata$type <- factor(xdata$type, levels = c("H","L","C","N")) #relevel
levels(xdata$type)[levels(xdata$type)=="H"] <- "F" # termed F (high fungal load)
levels(xdata$type)[levels(xdata$type)=="L"] <- "f" # termed f (low fungal load)

# unified color palette 
typecols= c("F"="#D95319","f"="#EDB116","C"="#999999","N"="#0072BD")

## ------STATS  
#*(a)account for pseudoreplication by averaging same individual treatment within a dish
#*(b)account for multiple testing 
library(rcompanion)
library(effectsize)
# With effect sizes. The matched-pairs rank biserial correlation coefficient (rc) is a recommended effect size statistic.  It is included in King, Rosopa, and Minimum (2000).
# NOTE When the data in the second group are greater than in the first group, rc is negative


TAB<-NULL
i=0
for (var in c("selfG","groomOut","groomIn","gaster","head")){

  for (typ in c("F","f","C","N")){
    i=i+1
    xdata$VARPRE<-xdata[[paste0(var,"PRE")]]
    xdata$VAR<-xdata[[var]]
    dat<-xdata%>%dplyr::select(dish, type, VAR,VARPRE)%>%
      subset(type==typ)%>%dplyr::group_by(dish)%>%
      dplyr::summarise(mean=mean(VAR),meanpre=mean(VARPRE) ) #*(a)
    W<-wilcox.test(dat$meanpre,dat$mean,paired=TRUE,conf.int = TRUE,exact = F)
    Y=c(dat$meanpre,dat$mean)
    G=factor(c(rep("pre",length(dat$meanpre)),
               rep("post",length(dat$mean))))
    rc<-round(effectsize(W)$r_rank_biserial,3)
    #rc<-wilcoxonPairedRC(x = Y ,g=G)
    TAB<-rbind(TAB,c(var,typ,W$p.value, W$statistic[[1]],round(W$estimate[[1]],5),round(W$conf.int[[1]],5),round(W$conf.int[[2]],5),rc))
  }}
TAB<-as.data.frame(TAB)
colnames(TAB)<-c("behaviour","ind.treat","Wilc.p.value","Vstatistic", "Pseudomedian","95CI[1]","95CI[2]", "rc")

TAB$adj.Wilc.p.value<-0
#*(b) adjust p-values within individual treatment
for (tr in c("F","f","C","N")){
  TAB[TAB$ind.treat==tr,]$adj.Wilc.p.value<-p.adjust(as.vector(TAB[TAB$ind.treat==tr,]$Wilc.p.value),method="BH")
}
TAB$sign<-ifelse(TAB$adj.Wilc.p.value<0.001,"***",ifelse(TAB$adj.Wilc.p.value<0.01,"**",ifelse(TAB$adj.Wilc.p.value<0.05,"*","n.s")))

TAB$adj.Wilc.p.value<-round(TAB$adj.Wilc.p.value,5)

STAB<-dplyr::select(TAB,c("behaviour","ind.treat","sign"))
table.sign<-tidyr::spread(STAB,ind.treat, sign)

table.sign<-table.sign[c(4,3,2,1),c(1,4,3,2,5)]
table.sign<-as.data.frame(table.sign)
table.sign$behaviour<-c("SD","PD","RD","GD")
table.sign<-as.matrix(table.sign)


## ------PLOTS 
scaleFUN <- function(x) sprintf("%.2f", x) # format, two decimal points on y-axis
y_labs<-c("Change in selfgrooming", "Change in performed grooming", 
          "Change in received grooming","Change in poison uptake behaviour")
ylimstable<-rbind(c(-0.02,0.055),c(-0.2,0.55),c(-0.2,0.55),c(-0.02,0.055),c(-0.02,0.055))

dat<-xdata%>%dplyr::select(dish,type,SD,PD,RD,GD,HD)%>%group_by(dish,type)%>%
  dplyr::summarise(SD=mean(SD),PD=mean(PD),RD=mean(RD),GD=mean(GD),drop=F) %>%as.data.frame() 

i=0
plot_list = list()
for(var in c("SD","PD","RD","GD")){
  i=i+1
  print(var)
  #get summary statistics
  dat_summ<-data_summary(dat,varname=var,groupnames=c("type")) 
  colnames(dat_summ)[5]<-"lo" # fix non-unique name, lower bound CI around mean
  colnames(dat_summ)[6]<-"hi" # fix non-unique name, upper uper CI around mean
  dat_summ$type<-factor(dat_summ$type,levels=c("F","f","C","N"))
  dat_summ$x<-c(1,2,3,4)
  dat_summ$sem <- dat_summ$sd/sqrt(dat_summ$n) #SEmean
  dat_summ$xM<-dat_summ$mean # mean
  this_sign_label<-table.sign[i,2:5]
  this_lbl_ypos<-max(dat_summ$hi)+2*max(dat_summ$sem)
  this_ylim<-ylimstable[i,]
  g<-ggplot(aes(y=xM,x=type,fill=type,color=type),data=dat_summ)+
    geom_point(size=4)+
    geom_rect(aes(ymin=lo,ymax=hi,xmin=x-0.3, xmax=x+0.3,alpha=0.5),size=0.1)+ #95%CI 
    geom_errorbar(aes(ymin=xM-sem, ymax=xM+sem,width=0.2),size=1)+  #sem error bars
    annotate("text", label=this_sign_label,x=c(1:4),y=rep(this_lbl_ypos,4),size=6)+ #test
    geom_abline(intercept=0,slope=0,color="black") + #line at zero
    scale_y_continuous(labels=scaleFUN)+theme(axis.ticks = element_line(color = "black"))+ #two decimals
    scale_fill_manual(values=typecols)+
    scale_color_manual(values=typecols)+
    theme_pubr(base_size=14)+
    theme(legend.position = "none",legend.title = element_blank(),
          axis.ticks = element_line(colour = "black", size = 0.25),# ticks black
          axis.ticks.length=unit(-0.15, "cm"), #tick size
          axis.text.x= element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),#ticks inward
          axis.text.y= element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")))+ #ticks inward
    ylab(paste0(y_labs[i]))+coord_cartesian(ylim=this_ylim)
  g<-g+xlab("Individual treatment") #x-label all figs
  plot_list[[i]] = g
}

library(cowplot)
plot_grid(plot_list[[1]], labels = c("b"),ncol=1)
plot_grid(plot_list[[2]], labels = c("c"),ncol=1)
plot_grid(plot_list[[3]], labels = c("d"),ncol=1)
plot_grid(plot_list[[4]], labels = c("S2"),ncol=1)



#sesionInfo()
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux 10 (buster)
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
# LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.3.5.so
# 
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] cowplot_1.1.1      DescTools_0.99.45  gmodels_2.18.1.1  
# [4] effectsize_0.7.0   rcompanion_2.4.15  openxlsx_4.2.5    
# [7] RColorBrewer_1.1-3 ggpubr_0.4.0       ggplot2_3.3.6     
# [10] plyr_1.8.7         dplyr_1.0.9        forcats_0.5.1     
# [13] multcomp_1.4-19    TH.data_1.1-1      MASS_7.3-57       
# [16] survival_3.3-1     mvtnorm_1.1-3      lme4_1.1-29       
# [19] Matrix_1.4-1       emmeans_1.7.5     