# Nestmate selfgrooming following performed grooming towards (Suppl. Fig. 3a), and received grooming from (Suppl. Fig. 3b)  others.  

## This R script uses the datafile created by 'selfG_triggers.R' named '3m_selfG_triggers.csv' (which can also be directly downloaded from the data repository) to both plot and statistically compare the self grooming after interactions with other spore-, control- or untreated ants. In particular, we test if selfgrooming of a nestmate depends on whether it previously performed or received allogrooming by a spore-treated (F or f, raw labels H or L), control-treated (C, raw label T for individual and Tx for group treatment) or another nestmate (N). For the control individuals, we were further interested, whether it plays a role if they were together in a nest with a spore-treated (C with F or f) or with another control individual (C with C), i.e. in a nest with or without pathogen threat.

rm(list = ls())
basedir='~/collective-care/Experimental_analyses/'
setwd(basedir)

# LIBRARIES---------------
library(tidyverse)
library(stringr)
library(ggpubr)
library(rstatix)
library(cowplot)

# COLOR PALETTE-----------
typefill<-c("Ff"="#FF5E05","CwFf"="#999999","CwC"="#999999","N"="#0072BD")
typecols= c("Ff"="#FF5E05","CwFf"="#FF5E05","CwC"="#999999","N"="#0072BD")


#PLOTS AND STATS

i=0
plot_list = list()
for (this_case in c('case_1', 'case_2')){
  i=i+1
  df<-read.csv('./data/3m_selfG_triggers_body_only.csv')
  df<-read.csv('./data/3m_selfG_triggers_head_only.csv')
  df<-read.csv('./data/3m_selfG_triggers.csv')
  df<-df%>%filter(case==this_case)
  #df%>%View()
  
  xlab<-if(this_case=='case_1'){'After performing grooming to'}else{'After receiving grooming from'}
  
  case<-this_case
  print(case)
  
  df$base <- factor(df$base, levels = c("HL","T","N"), labels=c("Ff","C","N"))
  
  # datapoint: one mean value per dish, per type
  df%>%group_by(treat,rep,base)%>%
    dplyr::summarise(mean_per_dish=mean(prop_selfgrooming))%>%
    data.frame()->summ.dish
  ggplot(aes(x=treat,y=mean_per_dish,fill=base),data=summ.dish)+ geom_boxplot()+ 
    labs(x=xlab,y="Self grooming by N")+
    theme_pubr()

  # datapoint: one mean value_per dish, per type, C divided by context
  df$var<-case_when(df$base %in%c("Ff")~"Ff",
                    df$base %in%c("C") & df$treat %in%c("HTx","LTx")~"CwFf",
                    df$base %in%c("C") & df$treat %in%c("TxTx")~"CwC",
                    df$base %in%c("N")~"N" )
  
  df$var <- factor(df$var, levels = c("Ff", "CwFf", "CwC","N"))
  
  df%>%group_by(treat,rep,var)%>%dplyr::summarise(mean_per_dish=mean(prop_selfgrooming))%>%data.frame()->summ.dish

  if (case=="case_1"){
    g<-ggplot(summ.dish, aes(x=var,y=mean_per_dish,col=var))+ geom_boxplot(aes(fill=var),alpha=0.6)+ #geom_jitter()+
      labs(x=xlab,y="Self grooming by N")+
      scale_x_discrete(labels=c("F,f", expression("C"[" wF,f"]), expression("C"[" wC"]),"N"))+
      scale_fill_manual(values=typefill)+
      scale_color_manual(values=typecols)+
      scale_y_continuous(limits=c(0,0.6))+
      annotate("text", label=c("a","a","ab","b"),x=c(1:4),y=rep(0.5,4),size=6)+ #paired Wilcoxon pre/post 
      theme_pubr(base_size=14)+theme(legend.position = "none",legend.title = element_blank(),
                                     axis.ticks = element_line(colour = "black", size = 0.25),# ticks black
                                     axis.ticks.length=unit(-0.15, "cm"), #ticks inward
                                     axis.text.x= element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),#ticks inward
                                     axis.text.y= element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")))#ticks inward
    plot_list[[i]] = g} 
  
  if (case=="case_2"){
    g<-ggplot(summ.dish, aes(x=var,y=mean_per_dish,col=var))+ geom_boxplot(aes(fill=var),alpha=0.6)+ #geom_jitter()+
      labs(x=xlab,y="Self grooming by N")+
      scale_x_discrete(labels=c("F,f", expression("C"[" wF,f"]), expression("C"[" wC"]),"N"))+
      scale_fill_manual(values=typefill)+
      scale_color_manual(values=typecols)+
      scale_y_continuous(limits=c(0,0.6))+
     
      annotate("text", label=c("ns"),x=c(2.5),y=rep(0.5,1),size=6)+ #paired Wilcoxon pre/post 
      theme_pubr(base_size=14)+theme(legend.position = "none",legend.title = element_blank(),
                                     axis.ticks = element_line(colour = "black", size = 0.25),# ticks black
                                     axis.ticks.length = unit(-0.15, "cm"), #ticks inward
                                     axis.title.y = element_blank(),
                                     axis.text.x= element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),#ticks inward
                                     axis.text.y= element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")))#ticks inward
    plot_list[[i]] = g} 
  
  
  #STATS
  print( kruskal.test(mean_per_dish~var, data = summ.dish))
  print(kruskal_effsize(mean_per_dish~var, data = summ.dish))
  print(dunn_test(mean_per_dish~var, data = summ.dish, p.adjust.method = "BH"))
  
}
p=plot_grid(plot_list[[1]],plot_list[[2]], labels = c("A","B"),ncol=2)
print(p)
# 
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux 10 (buster)
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
# LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.3.5.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] cowplot_1.1.1   rstatix_0.7.0   ggpubr_0.4.0    forcats_0.5.1  
# [5] stringr_1.4.0   dplyr_1.0.9     purrr_0.3.4     readr_2.1.2    
# [9] tidyr_1.2.0     tibble_3.1.7    ggplot2_3.3.6   tidyverse_1.3.1
