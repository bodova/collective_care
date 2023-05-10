# Scripts to wrangle G_all.csv (see Preparation/sampleAnalysis.py) into input for behaviour_aggregated.R

rm(list = ls())
basedir='~/collective_care/Experimental_analyses/'
setwd(basedir)
source(paste0(basedir,'aux/load_with_scripts.R'))

# data
xdata =read.table ("./data/G_ALL.csv",sep=",",header=TRUE)
# cleaning data
#replicate as factor
str(xdata);xdata$replicate<-as.factor(xdata$replicate) 
# replace N,X for N,M 
levels(xdata$level)<-c(levels(xdata$level),"M")
xdata$level[which(xdata$treatment == "N" & xdata$level == "X" )] <- "M"
# add dish, exposed
xdata$dish=interaction(xdata$treatmentgroup,xdata$replicate,drop=TRUE)
xdata$exposed=as.factor(paste0(xdata$treatment,xdata$level))
# relevel some factors
xdata$treatment <- factor(xdata$treatment, levels = c("G","R","T","N"))
xdata$exposed <- factor(xdata$exposed, levels = c("GH","RH","GL","RL","TX","NM"))
# add type
tmpexposed<-xdata$exposed
xdata$type<-revalue (tmpexposed,c("GH"="H","RH"="H","GL"="L","RL"="L","TX"="C","NM"="N"))
xdata$type <- factor(xdata$type, levels = c("H","L","C","N"))

# add treattype
xdata$treattype<-as.factor(paste0(xdata$treatmentgroup,'.',xdata$type))
xdata$treattype<-factor(xdata$treattype,levels=c("HH.H","HH.N","HL.H","HL.L","HL.N","HTx.H","HTx.T","HTx.N","LL.L","LL.N","LTx.L","LTx.T","LTx.N","TxTx.T","TxTx.N"))

xdata<-xdata%>%select(treatmentgroup, replicate,dish, period,
                      color, treatment,level,type,head,
                      gaster,selfG,groomIn,groomOut)
# split by period
predata<-xdata[xdata$period=="PRE",]
names(predata)[names(predata) == "gaster"] <- "gasterPRE"
names(predata)[names(predata) == "head" ] <-"headPRE"
names(predata)[names(predata) == "selfG"] <- "selfGPRE"
names(predata)[names(predata) == "groomIn"] <- "groomInPRE"
names(predata)[names(predata) == "groomOut"] <- "groomOutPRE"
predata<-subset(predata,select = -c(period))
posdata=xdata[xdata$period=="POST",]
posdata<-subset(posdata,select = -c(period))
xdata_wide<-join(predata,posdata)

# baselevel change  #PRE-->POST   POST-PRE
xdata_wide$HD=xdata_wide$head-xdata_wide$headPRE
xdata_wide$SD=xdata_wide$selfG-xdata_wide$selfGPRE
xdata_wide$PD=xdata_wide$groomOut-xdata_wide$groomOutPRE
xdata_wide$RD=xdata_wide$groomIn-xdata_wide$groomInPRE
xdata_wide$GD=xdata_wide$gaster-xdata_wide$gasterPRE

if (TRUE){ #-Inf to NA removes 225/594 
  xdata_wide$HD[is.infinite(xdata_wide$HD)]=NA
  xdata_wide$SD[is.infinite(xdata_wide$SD)]=NA
  xdata_wide$PD[is.infinite(xdata_wide$PD)]=NA
  xdata_wide$RD[is.infinite(xdata_wide$RD)]=NA
  xdata_wide$GD[is.infinite(xdata_wide$GD)]=NA
  
}
#xdata_wide<-na.omit(xdata_wide)
#write
#write.table (xdata_wide,"./data/behaviour_aggregated.csv",sep=",")

