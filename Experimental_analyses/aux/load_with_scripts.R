# Instead of loading everything every time, save this file somewhere and copy and  paste the below code into the beginning of new scripts: Usage: source('path_to_file/load_with_scripts.R')

#================================================================================#
########## functions ###
source('./aux/data_summary.R')
median_cl_boot <- function(x, conf = 0.95) {
  lconf <- (1 - conf)/2
  uconf <- 1 - lconf
  require(boot)
  bmedian <- function(x, ind) median(x[ind])
  bt <- boot(x, bmedian, 1000)
  bb <- boot.ci(bt, type = "perc")
  data.frame(y = median(x), ymin = quantile(bt$t, lconf), 
             ymax = quantile(bt$t,uconf))
}


########## libraries ###
##stats, math##
library(emmeans)
library(lme4)
library(multcomp)
library(rcompanion)
##data wrangling##
library(forcats)
library(dplyr)
library(plyr)
##plotting##
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
##other##
library(openxlsx)

######## themes for graphical purposes#####
mytheme1<- theme(legend.text = element_text(face = "italic",family = "Helvetica", size = rel(1)), 
                 axis.title = element_text(family = "Helvetica",size = rel(1)), 
                 axis.text = element_text(family = "Helvetica",size = rel(.9)), 
                 axis.line = element_line(size = 1,colour = "black"), 
                 axis.ticks = element_line(colour="black",size = rel(1)),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(), 
                 panel.background = element_rect(fill = "white"), 
                 legend.position="none",
                 legend.key = element_rect(fill = "whitesmoke"), 
                 legend.title = element_text(size = rel(1),family = "Helvetica"), 
                 plot.title = element_text(face = "bold",size = rel(1),family = "Helvetica"))

mytheme2<- theme(legend.text = element_text(face = "italic",family = "Helvetica", size = rel(1)), 
                 axis.title = element_text(family = "Helvetica",size = rel(1.5)), 
                 axis.text = element_text(family = "Helvetica",size = rel(1.5)), 
                 axis.line = element_line(size = 1,colour = "black"), 
                 axis.ticks = element_line(colour="grey",size = rel(1.4)),
                 panel.grid.major = element_line(colour="grey",size = rel(0.5)), 
                 panel.grid.minor = element_blank(), 
                 panel.background = element_rect(fill = "white"), 
                 legend.key = element_rect(fill = "whitesmoke"), 
                 legend.title = element_text(size = rel(1.5),family = "Helvetica"), 
                 plot.title = element_text(face = "bold",size = rel(1.7),family = "Helvetica"))

mytheme3<- theme(legend.text = element_text(face = "italic",family = "Helvetica", size = rel(1)), 
                 axis.title = element_text(family = "Helvetica",size = rel(1.5)), 
                 axis.text = element_text(family = "Helvetica",size = rel(1.5)), 
                 axis.line = element_line(size = 1,colour = "black"), 
                 axis.ticks = element_line(colour="grey",size = rel(1.4)),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), 
                 panel.background = element_rect(fill = "white"), 
                 legend.position="none",
                 legend.key = element_rect(fill = "whitesmoke"), 
                 legend.title = element_text(size = rel(1.5),family = "Helvetica"), 
                 plot.title = element_text(face = "bold",size = rel(1.7),family = "Helvetica"))

mytheme4<- theme(legend.text = element_text(face = "italic",family = "Helvetica", size = rel(1)), 
                 axis.title = element_text(family = "Helvetica",size = rel(1.5)), 
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank(),
                 axis.text = element_text(family = "Helvetica",size = rel(1)), 
                 axis.line = element_line(size = 1,colour = "black"), 
                 axis.ticks = element_line(colour="grey",size = rel(1.4)),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), 
                 panel.background = element_rect(fill = "white"), 
                 legend.position="none",
                 legend.key = element_rect(fill = "whitesmoke"), 
                 legend.title = element_text(size = rel(1.5),family = "Helvetica"), 
                 plot.title = element_text(face = "bold",size = rel(1.7),family = "Helvetica"))



