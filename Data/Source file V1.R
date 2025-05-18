

# Libraries ---------------------------------------------------------------

#ensure trainees install packages prior to loading these files

library(dplyr)
library(tidyr)
library(flexfit) # this needs to change
library(ggplot2)
library(epitools)
library(MASS)
#add more libraries here as needed

# Lab 1 function ---------------------------------------------------------

#read_and_tidy function

read_and_tidy<- function(file_name, plate_number, num_wells, antigen_names, control_samples, background_samples, standard_curve_values, bead_threshold){
  
  data<- read.csv(file_name)
  
  MFI_row<-which(data=="Median", arr.ind=TRUE)[1]
  
  CountVal_matrix<-which(data=="Count", arr.ind=TRUE)
  
  BeadCount_row<-CountVal_matrix[which(CountVal_matrix[,1]>MFI_row),1]
  
  MFI_matrix<- data[(MFI_row+2):(MFI_row+1+num_wells),1:(length(antigen_names)+3)]
  colnames(MFI_matrix)<- data[MFI_row+1,1:(length(antigen_names)+3)]
  
  MFI_long<- MFI_matrix|>
    dplyr::select(-c("Total Events"))|>
    tidyr::gather(key="Antigen", value="MFI", c(antigen_names))
  
  BeadCount_matrix<- data[(BeadCount_row+2):(BeadCount_row+1+num_wells),1:(length(antigen_names)+3)]
  colnames(BeadCount_matrix)<- data[BeadCount_row+1,1:(length(antigen_names)+3)]
  
  BeadCount_long<- BeadCount_matrix|>
    dplyr::select(-c("Total Events"))|>
    tidyr::gather(key="Antigen", value="BeadCount", c(antigen_names))
  
  output_df<- MFI_long|>
    left_join(BeadCount_long)|>
    mutate(Plate= plate_number, 
           Sample_Type= ifelse(Sample%in% background_samples, "BG", 
                               ifelse(Sample %in% control_samples, "Ctrl", 
                                      ifelse(Sample%in% standard_curve_values$Sample, "StdCurve", "TestSample"))))
  
  
  output_df$MFI<- as.numeric(output_df$MFI)
  output_df$BeadCount<- as.numeric(output_df$BeadCount)
  
  output_df<- output_df|>
    mutate(Low_Beads= ifelse(BeadCount<bead_threshold,1,0))
  return(output_df)
}


# Lab 2 function ---------------------------------------------------------


#filter low bead count
filter_low_beads <- function(bg_df){
  
  filtered_bg_df<- bg_df|>
    filter((Sample_Type!="StdCurve"&Low_Beads==0)|Sample_Type=="StdCurve")
  
  return(filtered_bg_df)
}

#remove background (subtraction or division)
rm_background <- function(clean_df, method){
  
  bg_df<- clean_df%>%
    filter(Sample_Type=="BG")|>
    group_by(Antigen)|>
    summarize(Median_BG= median(MFI))
  
  
  if(method=="division"){
    bg_rm_df<- clean_df|>
      left_join(bg_df)|>
      mutate(MFI_BG= MFI/Median_BG)
  }
  
  if(method=="subtraction"){
    bg_rm_df<- clean_df|>
      left_join(bg_df)|>
      mutate(MFI_BG= MFI-Median_BG)
  }
  
  return(bg_rm_df)
}

#mfi transformations
get_mfi_transforms <- function(filt_df){
  
  transformed_df<- data.frame()
  for(i in 1:length(antigen_names)){
    ag_df<- filt_df|>
      filter(Antigen==antigen_names[i]&Sample_Type=="TestSample"&MFI_BG>=0)
    
    boxcox_fit<- boxcox(MFI_BG~1, plotit=F, data=ag_df)
    lambda<- boxcox_fit$x[which(boxcox_fit$y==max(boxcox_fit$y))]
    
    
    ag_df<- ag_df|>
      mutate(Log_BG= log(MFI_BG), BoxCox_BG= (MFI_BG^lambda - 1)/lambda, 
             Lambda=lambda)
    
    transformed_df<- rbind(transformed_df, ag_df)
    
  }
  
  transformed_df_long<- transformed_df|>
    gather(Transform_Method, BG_Trans_Value, 'MFI_BG':'BoxCox_BG')|>
    mutate(Transform_Method_Clean= ifelse(Transform_Method=="MFI_BG", "None", 
                                          ifelse(Transform_Method=="Log_BG", "Natural Log", 
                                                 "Box-Cox")))
  annotate_df<- transformed_df_long|>
    group_by(Antigen)|>
    summarize(Lambda= paste0("Lambda=",round(Lambda[1], 3)), Transform_Method_Clean=Transform_Method_Clean[1])
  
  dist_plot<- ggplot(transformed_df_long, aes(x=BG_Trans_Value, color=Transform_Method_Clean))+
    geom_density()+theme_bw()+facet_wrap(~Antigen, scales='free')+xlab("Transformed MFI")+
    scale_color_brewer(palette="Dark2")+geom_text(data=annotate_df, mapping = aes(
      x=-Inf, y=-Inf, label=Lambda),hjust=-0.1, vjust=-1, color="black")+
    ggtitle(paste0("Plate ", transformed_df_long$Plate))+ylab("Density")
  
  return(dist_plot)
}

#Standardization

  #1. standard curves
plot_std_curves<- function(plate_norm_df, std_curve_values, input){
  
  
  if(input== "bgMFI"){
    std_curve_df<- plate_norm_df|>
      filter(Sample_Type=="StdCurve")|>
      left_join(std_curve_values)
    std_curve_df$Replicate<- as.factor(std_curve_df$Replicate)
    
    sample_summary<- plate_norm_df|>
      filter(Sample_Type=="TestSample")|>
      group_by(Antigen)|>
      summarize(MedianMFI= log(median(MFI_BG, na.rm=T)), 
              Upper95= log(quantile(MFI_BG, 0.95,na.rm=T)), 
              Lower95= log(quantile(MFI_BG, 0.05,na.rm=T)))
    
    std_curve_plot<- ggplot(data= std_curve_df, aes(x=log(Dilution), y= log(MFI_BG), group=Replicate, color=Replicate))+
      geom_point()+geom_line()+facet_wrap(~Antigen)+xlab("Log Dilution")+ylab("Log Background Corrected MFI")+
      ggtitle(paste0("Standard Curves for Plate ", std_curve_df$Plate))+theme_bw()+scale_color_brewer(palette="Paired")+
      geom_hline(data=sample_summary, aes(yintercept=MedianMFI), lty=2, color="red")+
      geom_hline(data=sample_summary, aes(yintercept=Upper95),  color="red")+
      geom_hline(data=sample_summary, aes(yintercept=Lower95),  color="red")
  }
  
  if(input=="MFI"){
    std_curve_df<- plate_norm_df|>
      filter(Sample_Type=="StdCurve")|>
      left_join(std_curve_values)
    
    std_curve_df$Replicate<- as.factor(std_curve_df$Replicate)
    
    sample_summary<- plate_norm_df|>
      filter(Sample_Type=="TestSample")|>
      group_by(Antigen)|>
      summarize(MedianMFI= log(median(MFI, na.rm=T)), 
              Upper95= log(quantile(MFI, 0.95,na.rm=T)), 
              Lower95= log(quantile(MFI, 0.05,na.rm=T)))
    
    std_curve_plot<- ggplot(data= std_curve_df, aes(x=log(Dilution), y= log(MFI), group=Replicate, color=Replicate))+
      geom_point()+geom_line()+facet_wrap(~Antigen)+xlab("Log Dilution")+ylab("Log Raw MFI")+
      ggtitle(paste0("Standard Curves for Plate ", std_curve_df$Plate))+theme_bw()+scale_color_brewer(palette="Paired")+
      geom_hline(data=sample_summary, aes(yintercept=MedianMFI), lty=2, color="red")+
      geom_hline(data=sample_summary, aes(yintercept=Upper95),  color="red")+
      geom_hline(data=sample_summary, aes(yintercept=Lower95),  color="red")
    
  }
  return(list(std_curve_plot, std_curve_df))
  
}


############################
FUNinv <- function(y, par) {
  -par["Scale"]*log(((par["Aup"] - par["Alow"])/
                       (y - par["Alow"]))^(1/par["a"]) - 1) + par["Xmid"]
}


get_concentration_FlexFit<- function(plate_df_norm, std_curve_values, input){
  antigen_list<- unique(plate_df_norm$Antigen)
  
  std_curve_df<- plate_df_norm|>
    filter(Sample_Type=="StdCurve")|>
    left_join(std_curve_values)
  
  std_curve_df<- as.data.frame(std_curve_df)
  
  sample_log_conc_df<- data.frame()
  curve_fit_df<-data.frame()
  for(i in 1:length(antigen_list)){
    if(input=="MFI"){
      
      antigen_curve_complete<- std_curve_df|>
        filter(Antigen==antigen_list[i])|>
        mutate(assay_id= Plate, analyte=Antigen, expected_conc= Dilution, fi= MFI,
               well_role="Standard", sample_id=Sample, replicate= Replicate,
               Log_Dilution= log(Dilution), Log_MFI_val= log(MFI))|>
        dplyr::select(assay_id, analyte, expected_conc, Dilution, fi, well_role, sample_id, 
                      replicate,Log_Dilution,Log_MFI_val) 
      
      
      
      std_curveFit_antigen<- flexfit::fitStd( std= antigen_curve_complete, 
                                              xvar="Log_Dilution", yvar="Log_MFI_val", interactive = F)
      
      
      
      if(is.null(std_curveFit_antigen$par)==F){
        sample_norm_df<- plate_df_norm|>
          filter(Sample_Type%in%c("Ctrl", "TestSample"), Antigen==antigen_list[i])|>
          mutate(log_MFI= log(MFI),
                 Log_Conc_bg=FUNinv(log(MFI),std_curveFit_antigen$par))
        
       
        sample_log_conc_df<- rbind(sample_log_conc_df, sample_norm_df)
      }
      else{
        sample_norm_df<- plate_df_norm|>
          filter(Sample_Type%in%c("Ctrl", "TestSample"), Antigen==antigen_list[i])|>
          mutate(Log_Conc_bg=NA, 
                 log_MFI= log(MFI))
        
    
        
        sample_log_conc_df<- rbind(sample_log_conc_df, sample_norm_df)
      }
      
    }
    
    
    
  
    
    if(input=="subMFI"){
      
      antigen_curve_complete<- std_curve_df|>
        filter(Antigen==antigen_list[i])|>
        mutate(assay_id= Plate, analyte=Antigen, expected_conc= Dilution, fi= MFI_BG_corrected,
               well_role="Standard", sample_id=Sample, replicate= Replicate,
               Log_Dilution= log(Dilution), Log_MFI_BG_val= log(MFI_BG_corrected))|>
        dplyr::select(assay_id, analyte, expected_conc, Dilution, fi, well_role, sample_id, 
                      replicate,Log_Dilution,Log_MFI_BG_val) 
      
      
      
      std_curveFit_antigen<- flexfit::fitStd( std= antigen_curve_complete, 
                                              xvar="Log_Dilution", yvar="Log_MFI_BG_val", interactive = F)
      
     
      
      
      if(is.null(std_curveFit_antigen$par)==F){
        sample_norm_df<- plate_df_norm|>
          filter(Sample_Type%in%c("Ctrl", "TestSample"), Antigen==antigen_list[i])|>
          mutate(Log_MFI_BG_val= log(MFI_BG_corrected),
                 Log_Conc_bg=FUNinv(log(MFI_BG_corrected),std_curveFit_antigen$par))
        
        sample_log_conc_df<- rbind(sample_log_conc_df, sample_norm_df)
        
      
      }
      else{
        sample_norm_df<- plate_df_norm|>
          filter(Sample_Type%in%c("Ctrl", "TestSample"), Antigen==antigen_list[i])|>
          mutate(Log_MFI_BG_val= log(MFI_BG_corrected),
                 Log_Conc_bg=NA)
        
        sample_log_conc_df<- rbind(sample_log_conc_df, sample_norm_df)
      
      }
      
    }

  }
  
  
  
  
  
  return(sample_log_conc_df)
  
}



#Normalization



get_norm_df<- function(long_df, method){
  if(method=="MFI"){
    norm_df_full=data.frame()
    norm_methods= c("log_MFI", "Log_Conc_bg")
    for(i in 1:length(norm_methods)){
      std_curve_fit= norm_methods[i]
      wide_df<- long_df|>
        dplyr::select(Location, Sample, Antigen, std_curve_fit, Plate, Sample_Type)|>
        pivot_wider(id_cols=c('Plate', 'Sample', 'Location', 'Sample_Type'), names_from= Antigen, values_from = std_curve_fit)|>
        filter_all(all_vars(!is.infinite(.)))
      
   
      ctrl_df<- long_df|>
        dplyr::select(Antigen, Plate, Sample, Sample_Type, std_curve_fit)|>
        filter(Sample_Type%in%c("StdCurve", "Ctrl"))|>
        filter_all(all_vars(!is.infinite(.)))
      
      n_obs_per_plate<- long_df|>
        group_by(Plate)|>
        summarize(n_obs=n())
      #simple LM with ctrls
      ctrl_df$Plate<- as.factor(ctrl_df$Plate)
      
      lm_ctrls<- lm(ctrl_df[,std_curve_fit]~Antigen+Plate+Sample, data= ctrl_df)
      plate_coefs<- lm_ctrls$coefficients[startsWith(names(lm_ctrls$coefficients),"Plate")]
      
      long_df<- long_df|>
        mutate(plate_effect= rep(c(0, plate_coefs), c(n_obs_per_plate$n_obs)))|>
        mutate(LM_norm_MFI= long_df[,std_curve_fit]-plate_effect)
      
      norm_df<-long_df|>
        dplyr::select(Plate, Sample, Location, Sample_Type, Antigen,std_curve_fit, LM_norm_MFI)|>
        pivot_longer(cols= !c('Plate', 'Sample', 'Location', 'Sample_Type', 'Antigen'), names_to="Norm_Method", values_to="Norm_MFI")|>
        mutate(Input_Value = std_curve_fit)
      
      norm_df_full<- rbind(norm_df_full, norm_df)
    }
    
  }
  
 
    
    
  
  return(norm_df_full)
}



