# We here analyse the pellets produced by the groups of ants in dependence of their treatment (Suppl. Fig. 5a), as well as the timeline for observed pellet expulsions for individual ants

## Note that only the 5 treatment groups with spore treatments are included, as none of the control (CC) groups produced any pellets 


### Note that original labels in the datasets contain "H" for high fungal load and "L" for low fungal load, which are termed F (for H) and f (for L) in the publication. Controls are labeled as 'C' in the datafile 'pellets_retrieved.csv' (which gives the numbers of pellets collected per experimental dish at the end of the experiment), and as 'T' in the data files 'final_spore_load.csv'(giving the measured spore load at the end of the experiment) and 'pellet-observations.csv' (giving the timeline of pellet expulsion events that could be observed by individual ants). Replicates are labelled according to their original dish number (note that an overall of 9 dishes across groups could not be analysed due to technical errors, so that only 99 replicates could be included in the study). The file 'labeling.csv' shows replicate number alongside with the original dish number. 

#================================================================================#

rm(list = ls())
basedir='~/collective-care/Experimental_analyses/'
setwd(basedir)
source(paste0(basedir,'aux/load_with_scripts.R'))


#-------------------------------------
#Number of pellets recovered from each group (Supp. Fig. 5a) 

pelletnum = read.table("./data/pellets_retrieved.csv",header=T,sep=",") 
pelletnum$Treatment <- factor(pelletnum$Treatment, levels = c("HH","HL","LL","HC","LC"))

levels(pelletnum$Treatment)[levels(pelletnum$Treatment)=="HH"] <- "FF"
levels(pelletnum$Treatment)[levels(pelletnum$Treatment)=="HL"] <- "Ff"
levels(pelletnum$Treatment)[levels(pelletnum$Treatment)=="LL"] <- "ff"
levels(pelletnum$Treatment)[levels(pelletnum$Treatment)=="HC"] <- "FC"
levels(pelletnum$Treatment)[levels(pelletnum$Treatment)=="LC"] <- "fC"

#install.packages("effectsize", repos = "https://easystats.r-universe.dev/") #for effect size
library("effectsize")
res=lm(Totpellets~Treatment,pelletnum)
n=lm(Totpellets~1,pelletnum)

eta_squared(res, partial = FALSE)
omega_squared(res)

anova(n,res,test="Chisq")
glht(res, linfct=mcp(Treatment="Tukey"))
cld(glht(res, linfct=mcp(Treatment="Tukey")))
summary(glht(res,linfct=mcp(Treatment="Tukey")),test=adjusted(type="BH"))
cld<-c("a","ab","ab","b","c")

pA<-ggplot(pelletnum,aes(x=Treatment,y=Totpellets))+
  geom_dotplot(binaxis="y",stackdir="center",binwidth=1,dotsize=.3,aes(fill=Treatment))+
  labs(x = "Treatment group",y="Number of pellets collected per group") +theme(legend.position="none")+
  scale_fill_manual(values = c("black", "gray40", "gray75","gray75","gray95"))+
  scale_y_continuous(breaks=seq(0,14,2),lim=c(0,13))+
  annotate("text", label=cld,x=c(1:5),y=rep(13,5), size=6)+
  theme_pubr(legend="none",base_size=14)+
  theme(legend.position = "none",legend.title = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.25),# ticks black
        axis.ticks.length=unit(-0.15, "cm"), #ticks inward
        axis.text.x= element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),#ticks inward
        axis.text.y= element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) #ticks inward


p=plot_grid(pA, labels = c("(A)"),ncol=1)
print(pA)

#-------------------------------------
#Mean number of spores in pellets is not different among treatment groups

## no figure associated, only statistics

### We get the number of spores from the final measured spore load of the pellet pool per dish, which is included in the datafile measured_spore_load, including the 'ReplicateNr'; the group 'treatment' (HH, HL, HT, LL, LT, TT); the 'Replicate' dish number; the 'Color' of the individual (labelled 'pellet' for pellet); the 'Part' of the ant (head, body, or pellet pool); the individual treatment 'IndTreat' with high (H), low (L) fungal load, the control treatment (T) or no treatment (N), with NA for pellets; the 'Label' of the spore variant (GFP, RFP or NA for nestmates and pellets); and the measured spores with the respective label 'RFP' or 'GFP'

pellspores = read.table("./data/final_spore_load.csv",header=T,sep=",") 
pellspores<-pellspores[pellspores$Color=="pellets",] #Colour
pellspores$Treatment<-as.factor(pellspores$Treatment)
pellspores$Totspores<-pellspores$RFP+pellspores$GFP
pellspores<-pellspores%>% dplyr::select(Treatment, Replicate,Totspores)


levels(pellspores$Treatment)[levels(pellspores$Treatment)=="HH"] <- "FF"
levels(pellspores$Treatment)[levels(pellspores$Treatment)=="HL"] <- "Ff"
levels(pellspores$Treatment)[levels(pellspores$Treatment)=="LL"] <- "ff"
levels(pellspores$Treatment)[levels(pellspores$Treatment)=="HT"] <- "FC"
levels(pellspores$Treatment)[levels(pellspores$Treatment)=="LT"] <- "fC"

pell<-merge(pelletnum,pellspores,by=c("Treatment","Replicate"))
pell<-subset(pell,Totspores>0 & Totpellets>0)
pell$Meanperpell<-pell$Totspores/pell$Totpellets

k<-kruskal.test(data=pell,Meanperpell~Treatment)
epsilonSquared(x = pell$Meanperpell, g = pell$Treatment)
pell$all<-rep(1)
dat_summ<-data_summary(pell,varname="Meanperpell",groupnames=c("all"))
dat_summ$sem <- dat_summ$sd/sqrt(dat_summ$n)
round(t(dat_summ),2)
data_summary(pell,varname="Meanperpell",groupnames=c("Treatment")) 

#-----------------------------------------------------------------------------------
#first pellets needs less time than subsequent (Suppl. Fig. 5b)

#tau is a measure of how much they groom before they drop a pellet
#M is a measure of when they have groomed halfway through before they drop a pellet

#annotated pellets being disgorged
pell<-read.csv("./data/pellet_observations.csv",header=T) 
pell<-pell%>%filter(!(Treatment=="HC"& Replicate==18))####These were not removed from submission
fps<-15
pell<-pell[complete.cases(pell), ]
pell<-pell%>%select("Treatment","Replicate","Color","Type","Number_correct","Relative_frame","Groom_between")
colnames(pell)<-c("treat","rep","color","type","numpellet","relframe","tau")

## note that these refer to: "Treatment" to treatment of the group (HH, HL, LL, HC, LC); "Replicate" to original dish number; "Color" to the color of the performing ant; "Type" to the individual ant treatment (high H, low L fungal load, T for TritonX control and N for no treatment in the nestmates); "Number_correct" as the number of the spit-out pellet, corrected for the fact that the ants sometimes picked up an already expelled pellet once more without intermediate grooming, hence these pellets were not produced-and-expelled, but only picked-and-expelled, and were thus removed from the analysis;"Relative_frame": frame of observation,"Groom_between": grooming period in-between consecutive expulsions.

typecols= c("H"="#D95319","L"="#EDB116","T"="#999999","N"="#0072BD")

# We here focus on the pellet production by nestmates (who produce most pellets and have themselves not been treated); we also exclude the first pellet, as the infrabuccal pocket filling state at the start of the experiment is unknown and hence the time grooming to first pellet not as defined as for the consecutive pellets
pell<-subset(pell,type=="N") #only nestmates
pell<-subset(pell,numpellet!=1) #exclude first pellet


# we look at the grooming nestmates do towards exposed individuals before pellet production
# we contrast how much grooming it takes for a pellet in first 18 min vs later
library(forcats)
pell$early<-rep("3")
pell$early[pell$relframe/(fps*60)<=60]<-"2" 
pell$early[pell$relframe/(fps*60)<=30]<-"1"
pell$min<-pell$relframe/(fps*60)
pell$early<-as.factor(pell$early)

# filter outliers 
pell<-pell%>%filter(pell$min<90)


median<-median(pell[pell$type=="N",]$relframe/(fps*60))
sum<-pell%>%dplyr::group_by(treat)%>%dplyr::summarise(Number=n())

pell%>%dplyr::group_by(early)%>%dplyr::summarise(Number=n())

sd<-subset(pell,type%in%c("N"))
k<-kruskal.test(tau/(fps*60)~early,sd)

library(rcompanion)
epsilonSquared(x = sd$tau/fps*60, g = sd$early)
pairwise.wilcox.test(sd$tau/fps*60,sd$early)
cld<-c("a","b","c")

pB<-ggplot(sd,aes(x=early,y=tau/(fps*60)))+
  geom_dotplot(binaxis="y",stackdir="center",binwidth=0.5,dotsize=1.2,fill="#0072BD")+
  labs(x = "Time after treatment (min)",y="Grooming between pellet expulsion events (min)")+
  scale_y_continuous(limits=c(0,32))+ #prev 25
  scale_x_discrete(labels=c("1" = "0-30", "2" = "30-60","3"="60-90"))+
  annotate("text", label=cld,x=c(1:3),y=rep(32,3), size=6)+theme_pubr(base_size=14)+
  theme(legend.position = "none",legend.title = element_blank(),
      axis.ticks = element_line(colour = "black", size = 0.25),# ticks black
      axis.ticks.length=unit(-0.15, "cm"), #ticks inward
      axis.text.x= element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),#ticks inward
      axis.text.y= element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm"))) #ticks inward
print(pB)

library(cowplot)
p=plot_grid(pB, labels = c("(B)"),ncol=1)

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
#   [1] rcompanion_2.4.15  DescTools_0.99.45  gmodels_2.18.1.1   cowplot_1.1.1      effectsize_0.7.0   openxlsx_4.2.5    
# [7] RColorBrewer_1.1-3 ggpubr_0.4.0       ggplot2_3.3.6      plyr_1.8.7         dplyr_1.0.9        forcats_0.5.1     
# [13] multcomp_1.4-19    TH.data_1.1-1      MASS_7.3-57        survival_3.3-1     mvtnorm_1.1-3      lme4_1.1-29       
# [19] Matrix_1.4-1       emmeans_1.7.5     

