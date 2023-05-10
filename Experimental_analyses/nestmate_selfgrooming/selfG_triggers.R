# Nestmate selfgrooming following performed or received grooming (Suppl. Fig. 3) 
 # (running time: ca. 20min)

## We here define a function to find certain types of events, involving certain individuals,
# one event preceding another (current) by a given lag. Start by filtering the events of interest,
# the individuals involved, and for every preceding event, look forward n windows (max_lag), and count how often there was a preceding event and how often there was none. 

### Note that original labels in the basic annotation data files (one file per video) contain "H" for high fungal load and "L" for low fungal load, which are termed "F" (for H) and "f" (for L) in the publication. The control group ("C" in the publication) is labelled by "Tx" when referring to the treatment group and by "T" in the individual ant treatments, given that it consisted of a sham treatment of Triton X. Data files are named according to their original dish number (note that an overall of 9 dishes across groups could not be analysed due to technical errors, so that only 99 replicates could be included in the study). The file 'labeling.csv' shows replicate number alongside the original dish number. 

#### Each file (given in tables_events) contains the detailed behavioural annotations. Allogrooming behaviour is labeled by 'G' and selfgrooming contains both grooming one owns body (S) or head (H). 


# SETUP------------------
rm(list = ls())
basedir='~/collective_care/Experimental_analyses/'
setwd(basedir)
source(paste0(basedir,'aux/load_with_scripts.R'))


# LIBRARIES---------------
library(tidyverse)
library(data.table)
library(stringr)
library(tictoc)
library(ggpubr)

# FUNCTION------------
find_fraction_curr_after_prec <- function(df,
                                        case, 
                                        max_lag,                        #in frames 
                                        events_prec,                     
                                        type_prec, type_R_prec,         
                                        events_curr,                     
                                        type_curr, type_R_curr           
                                        ) {
  # Initialize vars
  result<-data.frame(color_a=NULL,prop_selfgrooming=NULL,tot_duration=NULL)
 
  for (a in unique(df$color)) {
  
    # Subset the data frame to only include events involving a given individual
    idf <- df[df$color == a | df$color_R == a ,]
    
    if (nrow(idf) == 0) {
      print(paste0("no events involving ",a, " skipping..."))
      next}  
    
    # Get the rows where allogrooming was performed/received by given types (A+/A-)
    prec_rows <-NULL; 
    prec_rows<-idf%>%filter(event_type %in% events_prec)
    
    if (nrow(prec_rows)==0) {
      #print(paste0("no allogrooming by treated ",a, " skipping..."))
      next}
    
    #case 1 allogrooming performed towards group
    if (case == 'case_1' ) {
      prec_rows <- prec_rows %>% filter(color == a, type %in% type_prec & type_R %in% type_R_prec)}
    
    #case 2 allogrooming received from group
    if (case == 'case_2' ) {
      prec_rows <- prec_rows %>% filter(type %in% type_prec & type_R %in% type_R_prec & color_R == a)}
    
    if (nrow(prec_rows) == 0) {
      print(paste0("no events of interest involving ",a, " skipping..."))
      next}  
    
    # Get the end times of such events to iterate over
    prec_timestamps <- prec_rows$end 
    proportion_a<-NULL
    num_self_events<-0

    for (i in 1:length(prec_timestamps)) {

      # Get window bounds and length
      window_end<-prec_timestamps[i]+max_lag
      next_event<-prec_rows%>%filter(start>prec_timestamps[i] & start<window_end)
      if (nrow(next_event)>0){window_end=min(next_event$start)}
      window_length<-window_end-prec_timestamps[i]
      
      #Get selfgrooming events in window
      Sa<-idf%>%filter(event_type %in% events_curr & type%in%type_curr &
                                start>prec_timestamps[i] & end<window_end)
      
      num_self_events<-num_self_events+nrow(Sa)
      proportion_a[i]    <-sum(Sa$duration)/window_length
                       
    }
  
    result<-rbind(result,data.frame(color_a=a,
                                    prop_selfgrooming=proportion_a,
                                    num_self_events=num_self_events))
  
    }
  
  return (result)
}

# USE CASES: ----------------------------------------------------------------------------------
# # 
# # N selfgrooming (current event)
# type_curr=c('N');type_R_curr=c('N') 
# # 1) after performing allogrooming towards treated: case_1
# type_prec=c('N');type_R_prec=c('H','L','T') # 1) where N performed towards treated
# type_prec=c('N');type_R_prec=c('H','L') # 1.1) where N performed towards spore treated 

# # 2) after receiving allogrooming from treated : case_2
# type_prec=c('H','L','T');type_R_prec=c('N') # 2) where N received from treated
# type_prec=c('H','L');type_R_prec=c('N') # 2.1) where N received from spore treated


#---FIXED PARAMS----------------------------------
fps=15

#---PARAMETERS TO CHANGE--------------------------
max_lag=fps*60*3 #3min
treatments=c('HH','HL','HTx','LL','LTx','TxTx')
replicates=1:18


path=paste0(basedir,'data/tables_events/')


# Initialize result dataframes---------
case_1_data<-NULL
case_2_data<-NULL


tic()
for (per in c('POST')){
  for (treat in treatments ){
    for (rep in replicates){
      filename<- paste0(path,'EVENTS_',per,'_',treat,'_',rep,'.csv')
      # CHECK IF FILENAME EXISTS 
      if (!file.exists(filename)){print(paste0(per,'_',treat,'_',rep,' does not exist, skipping...'))}
      else{
        print(filename)
        # READ AND PREPARE DATA ALL TREATMENTS AND REPS----------
        events_df<-read.table(filename,sep=',')
        header<-c('treat','rep','start','duration','color','strain','dose','event_type','color_R','strain_R','dose_R')
        colnames(events_df)<-header
        events_df$type<-case_when(events_df$strain=='N'~'N',
                                  events_df$strain=='T'~'T',
                                  events_df$dose=='H'~'H',
                                  events_df$dose=='L'~'L')
        events_df$type_R<-case_when(events_df$strain_R=='N'~'N',
                                    events_df$strain_R=='T'~'T',
                                    events_df$dose_R=='H'~'H',
                                    events_df$dose_R=='L'~'L')
        events_df$end<-events_df$start+events_df$duration
        events_df<-events_df%>%filter (event_type%in%c('G','S','H'))%>%
          dplyr::select('treat','rep','start','duration', 'end','color','type','color_R','type_R','event_type')
        events_df$event_type<-dplyr::recode(events_df$event_type,'G'='allogrooming','H'='selfgrooming','S'='selfgrooming')
        events_df$per<-rep(per)
        
        A_HL_case_1<-NULL;A_T_case_1<-NULL;A_N_case_1<-NULL;
        A_HL_case_2<-NULL;A_T_case_2<-NULL;A_N_case_2<-NULL;
        
        # APPLY FUNCION
        A_HL_case_1 <- find_fraction_curr_after_prec (                                                                
                                                                events_df,
                                                                case='case_1',
                                                                max_lag =max_lag,
                                                                events_prec="allogrooming",
                                                                type_prec=c('N'),
                                                                type_R_prec=c('H','L'),
                                                                events_curr="selfgrooming",
                                                                type_curr=c('N'),
                                                                type_R_curr=c('N'))
    
        A_T_case_1 <- find_fraction_curr_after_prec (                                                                
                                                                  events_df,
                                                                  case='case_1',
                                                                  max_lag =max_lag,
                                                                  events_prec="allogrooming",
                                                                  type_prec=c('N'),
                                                                  type_R_prec=c('T'),
                                                                  events_curr="selfgrooming",
                                                                  type_curr=c('N'),
                                                                  type_R_curr=c('N'))
        
        A_N_case_1 <- find_fraction_curr_after_prec (                                                                
                                                                  events_df,
                                                                  case='case_1',
                                                                  max_lag =max_lag,
                                                                  events_prec="allogrooming",
                                                                  type_prec=c('N'),
                                                                  type_R_prec=c('N'),
                                                                  events_curr="selfgrooming",
                                                                  type_curr=c('N'),
                                                                  type_R_curr=c('N'))
        
        if(nrow(A_HL_case_1)>0){
          case_1_data  = rbind(case_1_data,
                               data.frame(per=rep(per),treat=rep(treat),rep=rep(rep),case=rep('case_1'),base='HL',
                                          color=A_HL_case_1$color_a, prop_selfgrooming=A_HL_case_1$prop_selfgrooming))}  
        if(nrow(A_T_case_1)>0){
          case_1_data  = rbind(case_1_data,
                               data.frame(per=rep(per),treat=rep(treat),rep=rep(rep),case=rep('case_1'),base='T',
                                          color=A_T_case_1$color_a, prop_selfgrooming=A_T_case_1$prop_selfgrooming))}   
        if(nrow(A_N_case_1)>0){
          case_1_data  = rbind(case_1_data,
                               data.frame(per=rep(per),treat=rep(treat),rep=rep(rep),case=rep('case_1'),base='N',
                                          color=A_N_case_1$color_a, prop_selfgrooming=A_N_case_1$prop_selfgrooming))} 
        A_HL_case_2 <- find_fraction_curr_after_prec (                                                                
                                                        events_df,
                                                        case='case_2',
                                                        max_lag =max_lag,
                                                        events_prec="allogrooming",
                                                        type_prec=c('H','L'),
                                                        type_R_prec=c('N'),
                                                        events_curr="selfgrooming",
                                                        type_curr=c('N'),
                                                        type_R_curr=c('N'))
        
        A_T_case_2 <- find_fraction_curr_after_prec (                                                                
                                                        events_df,
                                                        case='case_2',
                                                        max_lag =max_lag,
                                                        events_prec="allogrooming",
                                                        type_prec=c('T'),
                                                        type_R_prec=c('N'),
                                                        events_curr="selfgrooming",
                                                        type_curr=c('N'),
                                                        type_R_curr=c('N'))
        
        A_N_case_2 <- find_fraction_curr_after_prec (                                                                
                                                        events_df,
                                                        case='case_2',
                                                        max_lag =max_lag,
                                                        events_prec="allogrooming",
                                                        type_prec=c('N'),
                                                        type_R_prec=c('N'),
                                                        events_curr="selfgrooming",
                                                        type_curr=c('N'),
                                                        type_R_curr=c('N'))
        if(nrow(A_HL_case_2)>0){
          case_2_data  = rbind(case_2_data,
                             data.frame(per=rep(per),treat=rep(treat),rep=rep(rep),case=rep('case_2'),base='HL',
                                         color=A_HL_case_2$color_a, prop_selfgrooming=A_HL_case_2$prop_selfgrooming))}  
        if(nrow(A_T_case_2)>0){
          case_2_data  = rbind(case_2_data,
                               data.frame(per=rep(per),treat=rep(treat),rep=rep(rep),case=rep('case_2'),base='T',
                                          color=A_T_case_2$color_a, prop_selfgrooming=A_T_case_2$prop_selfgrooming))}   
        if(nrow(A_N_case_2)>0){
          case_2_data  = rbind(case_2_data,
                               data.frame(per=rep(per),treat=rep(treat),rep=rep(rep),case=rep('case_2'),base='N',
                                          color=A_N_case_2$color_a, prop_selfgrooming=A_N_case_2$prop_selfgrooming))}   
        
          }                              
      }#rep
    }#treat
  }#period
toc()

final_data<-rbind(case_1_data,case_2_data)

write.csv(final_data,file='./3m_selfG_triggers_head_only.csv', row.names=F)

# #--------------
# sessionInfo()
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux 10 (buster)
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
# LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.3.5.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] tictoc_1.1        data.table_1.14.2 cowplot_1.1.1     rstatix_0.7.0     ggpubr_0.4.0      forcats_0.5.1     stringr_1.4.0    
# [8] dplyr_1.0.9       purrr_0.3.4       readr_2.1.2       tidyr_1.2.0       tibble_3.1.7      ggplot2_3.3.6     tidyverse_1.3.1  
