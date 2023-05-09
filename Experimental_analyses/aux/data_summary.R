data_summary <- function(data, varname, groupnames){
  require(plyr)
  require(gmodels) #ci for mean, and confidence interval around the mean (CI lower and CI upper)
  require(DescTools) #ci for median, and confidence intervals around the median (lwr.ci and upr.ci)
  summary_func <- function(x, varname){
    c(mean = mean(x[[ varname]], na.rm=TRUE),
      n =  sum(!is.na(x[[ varname]])),
      sd = sd(x[[ varname]], na.rm=TRUE),
      #ci(x[[ varname]],na.rm=TRUE)[2],
      #ci(x[[ varname]],na.rm=TRUE)[3],
      #MeanCI(x[[ varname]],na.rm=TRUE,method="boot")[1], 
      MeanCI(x[[ varname]],na.rm=TRUE,method="boot")[2], 
      MeanCI(x[[ varname]],na.rm=TRUE,method="boot")[3],
      MedianCI(x[[ varname]],na.rm=TRUE,method="boot")[1], 
      MedianCI(x[[ varname]],na.rm=TRUE,method="boot")[2], 
      MedianCI(x[[ varname]],na.rm=TRUE,method="boot")[3],
      max= max(x[[ varname]]))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  return(data_sum)
}

#To test
#dat<-iris
#dat$Factor<-as.factor(dat$Petal.Width>1)
#dat_summ<-data_summary(dat,varname="Petal.Width",groupnames=c("Species","Factor"),which="mean") 
#test NA
#dat$Sepal.Length[1:10]<-NA #introduce NA
#dat_summ<-data_summary(dat,varname="Sepal.Length",groupnames=c("Species"),which="median")
#head(dat_summ)

