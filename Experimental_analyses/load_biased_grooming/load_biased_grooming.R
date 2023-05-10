# Grooming bias by the nestmates (Figure 2), in particular grooming choice (Fig 2c) and duration (Fig. 2d) as a function of current spore load difference. 

## Note that original labels in the data files 'ratios_real' and 'ratios_null' contain "H" for high fungal load and "L" for low fungal load, which are termed F (for H) and f (for L) in the publication. The control treatment (C) is labeled by T (Triton X). The 180 30sec windows of the experiment are here labeled from 0-179 (1-180 in the source data accompanying the publication). Further, the replicates are labelled according to their original dish number (note that an overall of 9 dishes across groups could not be analysed due to technical errors, so that only 99 replicates could be included in the study). The file 'labeling.csv' shows replicate number alongside with the original dish number. 

#For each grooming choice there is a given ratio between the current spore loads of the two treated ants(L_target/(L_target + L_other))

#We simulated random choices between target and non-target to get the expected distribution of grooming events wrt.  spore load ratios, and then compare expected (null) and observed (real) distribution, via Kolmogorov-Smirnov test

#The expected and observed ratios were computed using python scripts (see Preparation, in this repository).

rm(list = ls())
basedir='~/collective_care/Experimental_analyses/'
setwd(basedir)
source(paste0(basedir,'aux/load_with_scripts.R'))

# load observed and expected ratios (Fig. 2c) 
observed<-read.table ('./data/ratios_real.csv',sep=",",header=T) 
observed<-observed%>%subset(treatment%in%c("HH","HL","LL"))%>%droplevels()
expected<-read.table ('./data/ratios_null.csv',sep=",",header=T)
expected<-expected%>%subset(treatment%in%c("HH","HL","LL"))%>%droplevels()

observed<-observed%>%dplyr::select(treatment, replicate, load_ratio)
expected<-expected%>%dplyr::select(treatment, replicate, load_ratio)
observed$Simulated<-rep(0)
expected$Simulated<-rep(1)

df <- data.frame(rbind(observed,expected))
df$Simulated<-as.factor(df$Simulated)

binnum=10
bbreaks=seq(range(df$load)[1],range(df$load)[2],by= diff(range(df$load_ratio))/binnum )

# get difference between histograms 
plot<-ggplot(df,aes(x=load_ratio)) + 
  geom_histogram(data=subset(df,Simulated == 0),breaks=bbreaks,color="gray30",fill = "gray40", alpha = 1,aes(y = stat(count / sum(count))))+ 
  geom_histogram(data=subset(df,Simulated == 1),breaks=bbreaks,color="gray30",fill = "gray80", alpha = 0.3,aes(y = stat(count / sum(count))))+
  labs(y="Probability",x="Spore proportion on groomed ant")+mytheme3
plot_build <- ggplot_build(plot)
difference<-plot_build[["data"]][[1]][["y"]]-plot_build[["data"]][[2]][["y"]]
df2<-data.frame(cbind(difference,plot_build[["data"]][[1]][["x"]]))

# plot
plot<-ggplot(df,aes(x=load_ratio)) + 
  geom_histogram(data=subset(df,Simulated == 0),
                 breaks=bbreaks,color="gray50",fill = "gray70", alpha = 1,aes(y = stat(count / sum(count)))) 
binwidth = layer_data(plot) %>% mutate(w=xmax-xmin) %>% pull(w) %>% median
plot<-plot + stat_bin(data=subset(df,Simulated == 1),breaks=bbreaks[1:11], aes(y=stat(count/sum(count))),pad=TRUE,geom="step",position=position_nudge(x=-0.5*binwidth),size=1)+scale_x_continuous(limit=c(0.07,0.93))+
  labs(y="Probability",x="Spore proportion on groomed ant")+mytheme3
print(plot)
p2<-plot+geom_smooth(data=df2,aes(difference,x=V2),color="#00BB87",size=2,span=1,fill=NA) +
  theme(panel.background=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(color = "black"),axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5)) +
  annotate("text",x=0.5,y=.15,label="***",size=8)+
  annotate("text",x=.25,y=-0.005,label="groom lower",size=5,fontface = 'italic')+
  annotate("text",x=.75,y=-0.005,label="groom higher",size=5,fontface = 'italic')+
  geom_segment(aes(x=bbreaks[1], xend = bbreaks[1], y = 0, yend = 0.0362), color = "black",size=1)+
  coord_cartesian(ylim=c(-0.017, 0.17))
print(p2)

#duration heatmap (Fig. 2d) 
#reload data
observed<-read.table ('./data/ratios_real.csv',sep=",",header=T)
observed<-subset(observed,treatment%in%c("HL","HH","LL"))
cap<-120
fps<-15
sum(observed$duration<(fps*cap))/length(observed$duration) 
observed<-observed[observed$duration<=fps*cap,]
observed$duration<-observed$duration/fps #fps=15 convert from frames to seconds 60sec

# create plot
rx=range(observed$load_ratio)
binsizex=(rx[2]-rx[1])/(binnum)
ry=range(observed$duration)
ry[2]<-cap #upper limit is one minute
binsizey=(ry[2]-ry[1])/(binnum)

g1 <- ggplot(observed, aes(load_ratio, duration)) +
  geom_bin2d(binwidth = c(.1,cap*.1))+ # modify bin width
  labs(y="Grooming event duration (s)",x="Current spore load ratio")

pg<-ggplot_build(g1)
pg<-as.data.frame(pg$data)
dpg<-pg%>%group_by(ybin)%>%dplyr::summarise(sum=sum(value))
pg<-join(pg,dpg)
pg$nvalue<-pg$value/pg$sum

p3<-ggplot(pg, aes(x,y,fill=value))+geom_tile()+
  scale_fill_distiller(palette = "Greys", direction=1) + 
  labs(y="Grooming event duration (s)",x="Current spore load ratio",
       fill="Number of grooming events")+theme_pubr(legend="top")+ 
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(limits=c(0,121),breaks=seq(0,120,20),expand =c(0,0))

library(cowplot)
plot_grid(p2, labels = c("G"), ncol = 1, rel_heights=c(5, 5), align = 'v', axis = 'l')
plot_grid(p3, labels = c("H"), ncol = 1, rel_heights=c(5, 5), align = 'v', axis = 'l')

#STATISTICS
observed<-read.table ('./data/ratios_real.csv',sep=",",header=T)
observed<-subset(observed,treatment%in%c("HL","HH","LL"))
expected<-read.table ('./data/ratios_null.csv',sep=",",header=T)
expected<-subset(expected,treatment%in%c("HL","HH","LL"))

#ECDF PLOT
ecdf_plot=ggplot() +
  geom_step(data = observed, aes(x = load_ratio, y = ecdf(load_ratio)(load_ratio)),color = "gray50",size=1.5) +
  geom_step(data = expected, aes(x = load_ratio, y = ecdf(load_ratio)(load_ratio)),color = "black") +
  xlim(range(c(observed$load_ratio, expected$load_ratio))) +
  coord_fixed(ratio = 1, xlim=c(0,1), ylim = c(0,1)) +
  xlab("Load Ratio") +
  ylab("Cumulative Probability") +
  theme_pubr()

# QQPLOT 
qq.out <- qqplot(x=expected$load_ratio, y=observed$load_ratio, plot.it=FALSE)
qq.out <- as.data.frame(qq.out)
xylim <- range( c(qq.out$x, qq.out$y) )

qq_plot=ggplot(qq.out, aes( x= x, y = y)) + 
  geom_point(size=.5,color="#00BB87") + 
  geom_abline( intercept=0, slope=1) +
  coord_fixed(ratio = 1, xlim=c(0,1), ylim = c(0,1)) +
  ylab("Expected Quantiles") + xlab("Observed Quantiles")+theme_pubr()

library(cowplot)
plot_grid(ecdf_plot,qq_plot, labels = c("a","b"), ncol = 2, rel_heights=c(5, 5))

#library("Matching") #ks.boot 
#Note: each execution has a slightly different outcome, ran the test several times to ensure robustness
#ks.boot(observed$load_ratio,expected$load_ratio,nboots=1000,alternative="two.sided")
#D = 0.053304, p-value = 1.632e-06

###does duration of grooming events depend on current load ratio? 
observed<-read.table ('./data/ratios_real.csv',sep=",",header=T)
observed<-subset(observed,treatment%in%c("HL","HH","LL"))
# ks tests, pairwise between columns in the heatmap
observed$bin<-cut(observed$load_ratio, breaks = seq(0,1, by=.1),include.lowest = TRUE)
observed$binnum<-cut(observed$load_ratio, breaks = seq(0,1, by=.1), labels=FALSE,include.lowest = TRUE)
breaks = seq(0,1, by=.1)

resultKS<-NULL
for (i in c(1:10)){
for (j in c(1:10)){
  if(i==j){break}
  x<-observed[observed$binnum==i,]$duration
  y<-observed[observed$binnum==j,]$duration
  K<-ks.test(x,y, exact=NULL)
  if (K$p.value<0.05){
    print(paste0(i,"-th column is potentially different from column ", j, "p-value:",round(K$p.value,2)))
  }
  resultKS<-rbind(resultKS, K$p.value)
}
}

# Inspect potentially different columns
#ggplot(data=subset(observed,binnum%in%c(6,3,5,8)),aes(duration))+geom_density()+facet_grid(bin~.)+theme_pubr()
#resultsKSadjusted<-p.adjust(resultKS,method="bonferroni")
#print(resultsKSadjusted<0.05)# 

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
#   [1] cowplot_1.1.1      openxlsx_4.2.5     RColorBrewer_1.1-3 ggpubr_0.4.0       ggplot2_3.3.6      plyr_1.8.7        
# [7] dplyr_1.0.9        forcats_0.5.1      multcomp_1.4-19    TH.data_1.1-1      MASS_7.3-57        survival_3.3-1    
# [13] mvtnorm_1.1-3      lme4_1.1-29        Matrix_1.4-1       emmeans_1.7.5   
