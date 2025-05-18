

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
    
    
    
    std_curve_plot<- ggplot(data= std_curve_df, aes(x=log(Dilution), y= log(MFI_BG), group=Replicate, color=Replicate))+
      geom_point()+geom_line()+facet_wrap(~Antigen)+xlab("Log Dilution")+ylab("Log Background Corrected MFI")+
      ggtitle(paste0("Standard Curves for Plate ", std_curve_df$Plate))+theme_bw()+scale_color_brewer(palette="Paired")
  }
  
  if(input=="MFI"){
    std_curve_df<- plate_norm_df|>
      filter(Sample_Type=="StdCurve")|>
      left_join(std_curve_values)
    
    std_curve_df$Replicate<- as.factor(std_curve_df$Replicate)
    
    std_curve_plot<- ggplot(data= std_curve_df, aes(x=log(Dilution), y= log(MFI), group=Replicate, color=Replicate))+
      geom_point()+geom_line()+facet_wrap(~Antigen)+xlab("Log Dilution")+ylab("Log Raw MFI")+
      ggtitle(paste0("Standard Curves for Plate ", std_curve_df$Plate))+theme_bw()+scale_color_brewer(palette="Paired")
    
  }
  return(list(std_curve_plot, std_curve_df))
  
}

#do we need functions for LJ plots? 

  #2. concentrations

#uses flex fit
# which is the correct function? get_concentration or get_conc_combined_curve
get_concentration<- function(plate_df_norm, std_curve_values, input, replicate_type){
  antigen_list<- unique(plate_df_norm$Antigen)
  
  std_curve_df<- plate_df_norm|>
    filter(Sample_Type=="StdCurve")|>
    left_join(std_curve_values)
  
  
  std_curve_df<- as.data.frame(std_curve_df)
  
  
  curve_fit_df<-data.frame()
  for(i in 1:length(antigen_list)){
    for(j in 1:length(replicate_type)){
      if(replicate_type[j]=="HylE"&antigen_list[i]!="Typhoid HlyE"){
        curve_fit<- data.frame(sample= replicate_type[j],
                               method=c("flexFit", "nCal Regular", "ncal Bayes", "nplr"),
                               xmid= NA,
                               upper_lim= NA,
                               lower_lim= NA,
                               Antigen= antigen_list[i])
        
        curve_fit_df<- rbind(curve_fit_df, curve_fit)
      }
      else{
        if(input=="MFI"){
          
          antigen_curve_complete<- std_curve_df|>
            filter(Antigen==antigen_list[i]&startsWith(Replicate, replicate_type[j]))|>
            mutate(assay_id= Plate, analyte=Antigen, expected_conc= Dilution, fi= MFI,
                   well_role="Standard", sample_id=Sample, replicate= str_sub(Replicate, -1),
                   Log_Dilution= log(Dilution), Log_MFI_val= log(MFI),prop_mfi= convertToProp(MFI))|>
            dplyr::select(assay_id, analyte, expected_conc, Dilution, fi, well_role, sample_id, 
                          replicate,Log_Dilution,Log_MFI_val, prop_mfi, replicate) 
          
          
          
          std_curveFit_antigen<- flexfit::fitStd(std= antigen_curve_complete, 
                                                 xvar="Log_Dilution", yvar="Log_MFI_val", interactive = F)
          
          combined_ncal_fit<- nCal::ncal(log(fi)~ expected_conc, antigen_curve_complete, 
                                         return.fits=T)
          
          joint_ncal_fit<- nCal::ncal(log(fi)~ expected_conc, antigen_curve_complete, return.fits=T,
                                      bcrm.fit = T, plot=F)
          
          
          
          combined_ncal_params <- attr(combined_ncal_fit, "fits")[[1]]
          combined_ncal_coefs<- combined_ncal_params$coefficients
          
          joint_ncal_params <- attr(joint_ncal_fit, "fits")
          
          direct_compare_combined<- cla2gh(coef(combined_ncal_params))
          
          
          npfit<- nplr(x= antigen_curve_complete$Dilution, y=antigen_curve_complete$prop_mfi)
          np_params<- getPar(npfit)
          
          
          
          
          if(is.null(std_curveFit_antigen$par)==F){
            
            
            curve_fit<- data.frame(sample= replicate_type[j],
                                   method=c("flexFit", "nCal Regular", "ncal Bayes", "nplr"),
                                   xmid= unlist(c(unname(std_curveFit_antigen$par[3]), unname(direct_compare_combined[3]),
                                                  unname(joint_ncal_params$coefficients[4]),unname(log(10^np_params$params[3])))),
                                   upper_lim= unlist(c(unname(std_curveFit_antigen$par[1]),unname(direct_compare_combined[2]),
                                                       unname(joint_ncal_params$coefficients[2]),unname(np_params$params[2]))),
                                   lower_lim= unlist(c(unname(std_curveFit_antigen$par[2]),unname(direct_compare_combined[1]),
                                                       unname(joint_ncal_params$coefficients[1]),unname(np_params$params[1]))),
                                   Antigen= antigen_list[i])
            
            curve_fit_df<- rbind(curve_fit_df, curve_fit)
            
          }
          else{
            
            curve_fit<- data.frame(sample= replicate_type[j],
                                   method=c("flexFit", "nCal Regular", "ncal Bayes", "nplr"),
                                   xmid= unlist(c(NA, unname(direct_compare_combined[3]),
                                                  unname(joint_ncal_params$coefficients[4]),unname(log(10^np_params$params[3])))),
                                   upper_lim= unlist(c(NA,unname(direct_compare_combined[2]),
                                                       unname(joint_ncal_params$coefficients[2]),unname(np_params$params[2]))),
                                   lower_lim= unlist(c(NA,unname(direct_compare_combined[1]),
                                                       unname(joint_ncal_params$coefficients[1]),unname(np_params$params[1]))),
                                   Antigen= antigen_list[i])
            
            curve_fit_df<- rbind(curve_fit_df, curve_fit)
            
            
          }
          
        }
        
        
        
        if(input=="bgMFI"){
          
          antigen_curve_complete<- std_curve_df|>
            filter(Antigen==antigen_list[i]&startsWith(Replicate, replicate_type[j]))|>
            mutate(assay_id= Plate, analyte=Antigen, expected_conc= Dilution, fi= MFI_BG,
                   well_role="Standard", sample_id=Sample, replicate= str_sub(Replicate, -1),
                   Log_Dilution= log(Dilution), Log_MFI_bg_val= log(MFI_BG),prop_mfi= convertToProp(MFI))|>
            dplyr::select(assay_id, analyte, expected_conc, Dilution, fi, well_role, sample_id, 
                          replicate,Log_Dilution,Log_MFI_bg_val,prop_mfi) 
          
          
          
          std_curveFit_antigen<- flexfit::fitStd( std= antigen_curve_complete, 
                                                  xvar="Log_Dilution", yvar="Log_MFI_bg_val", interactive = F)
          
          combined_ncal_fit<- nCal::ncal(log(fi)~ expected_conc, antigen_curve_complete, 
                                         return.fits=T)
          
          joint_ncal_fit<- nCal::ncal(log(fi)~ expected_conc, antigen_curve_complete, return.fits=T,
                                      bcrm.fit = T, plot=F)
          
          
          
          combined_ncal_params <- attr(combined_ncal_fit, "fits")[[1]]
          combined_ncal_coefs<- combined_ncal_params$coefficients
          
          joint_ncal_params <- attr(joint_ncal_fit, "fits")
          
          
          direct_compare_combined<- cla2gh(coef(combined_ncal_params))
          
          
          npfit<- nplr(x= antigen_curve_complete$Dilution, y=antigen_curve_complete$prop_mfi)
          np_params<- getPar(npfit)
          
          
          
          if(is.null(std_curveFit_antigen$par)==F){
            
            
            
            
            curve_fit<- data.frame(sample=replicate_type[j],
                                   method=c("flexFit", "nCal Regular", "ncal Bayes", "nplr"),
                                   xmid= unlist(c(unname(std_curveFit_antigen$par[3]), unname(direct_compare_combined[3]),
                                                  unname(joint_ncal_params$coefficients[4]),unname(log(10^np_params$params[3])))),
                                   upper_lim= unlist(c(unname(std_curveFit_antigen$par[1]),unname(direct_compare_combined[2]),
                                                       unname(joint_ncal_params$coefficients[2]),unname(np_params$params[2]))),
                                   lower_lim= unlist(c(unname(std_curveFit_antigen$par[2]),unname(direct_compare_combined[1]),
                                                       unname(joint_ncal_params$coefficients[1]),unname(np_params$params[1]))),
                                   Antigen= antigen_list[i])
            
            curve_fit_df<- rbind(curve_fit_df,curve_fit)
          }
          else{
            
            
            curve_fit<- data.frame(sample= replicate_type[j],
                                   method=c("flexFit", "nCal Regular", "ncal Bayes", "nplr"),
                                   xmid= unlist(c(NA, unname(direct_compare_combined[3]),
                                                  unname(joint_ncal_params$coefficients[4]),unname(log(10^np_params$params[3])))),
                                   upper_lim= unlist(c(NA,unname(direct_compare_combined[2]),
                                                       unname(joint_ncal_params$coefficients[2]),unname(np_params$params[2]))),
                                   lower_lim= unlist(c(NA,unname(direct_compare_combined[1]),
                                                       unname(joint_ncal_params$coefficients[1]),unname(np_params$params[1]))),
                                   Antigen= antigen_list[i])
            
            curve_fit_df<- rbind(curve_fit_df, curve_fit)
          }
          
        }
      }
      print(replicate_type[j])
    }
    print(antigen_list[i])
  }
  
  
  
  
  
  return(list(curve_fit_df))
  
}


get_conc_combined_curve<- function(complete_plates, std_curve_complete, input){
  
  std_curve_df<- complete_plates|>
    filter(Sample_Type=="StdCurve")|>
    left_join(std_curve_complete)
  
  std_curve_df<- as.data.frame(std_curve_df)
  
  
  antigen_list<- unique(complete_plates$Antigen)
  sample_log_conc_df<- data.frame()
  for(i in 1:length(antigen_list)){
    if(input=="MFI"){
      
      antigen_curve_complete<- std_curve_df|>
        filter(Antigen==antigen_list[i])|>
        mutate(assay_id= Plate, analyte=Antigen, expected_conc= Dilution, fi= MFI,
               well_role="Standard", sample_id=Sample, replicate= Replicate,
               Log_Dilution= log(Dilution), Log_MFI_val= log(MFI))|>
        dplyr::select(assay_id, analyte, expected_conc, Dilution, fi, well_role, sample_id, 
                      replicate,Log_Dilution,Log_MFI_val) 
      
      
      
      
      joint_ncal_fit<- nCal::ncal(log(fi)~ expected_conc, antigen_curve_complete, return.fits=T,
                                  bcrm.fit = T, plot=F,error.model = "gh_t4")
      
      
      
      joint_ncal_params <- attr(joint_ncal_fit, "fits")
      joint_ncal_coefs<- joint_ncal_params$coefficients
      
      
      sample_norm_df<- complete_plates|>
        filter(Sample_Type%in%c("Ctrl", "TestSample"), Antigen==antigen_list[i])|>
        mutate(Log_Conc_Joint= getConc(joint_ncal_params, log(MFI))[,1])
      
      sample_log_conc_df<- rbind(sample_log_conc_df, sample_norm_df)
    }
    
    if(input=="bgMFI"){
      
      antigen_curve_complete<- std_curve_df|>
        filter(Antigen==antigen_list[i])|>
        mutate(assay_id= Plate, analyte=Antigen, expected_conc= Dilution, fi= MFI_BG,
               well_role="Standard", sample_id=Sample, replicate= Replicate,
               Log_Dilution= log(Dilution), Log_MFI_val_bg= log(MFI_BG))|>
        dplyr::select(assay_id, analyte, expected_conc, Dilution, fi, well_role, sample_id, 
                      replicate,Log_Dilution,Log_MFI_val_bg) 
      
      
      
      
      joint_ncal_fit<- nCal::ncal(log(fi)~ expected_conc, antigen_curve_complete, return.fits=T,
                                  bcrm.fit = T, plot=F,error.model = "gh_t4")
      
      
      
      joint_ncal_params <- attr(joint_ncal_fit, "fits")
      joint_ncal_coefs<- joint_ncal_params$coefficients
      
      
      sample_norm_df<- complete_plates|>
        filter(Sample_Type%in%c("Ctrl", "TestSample"), Antigen==antigen_list[i])|>
        mutate(Log_Conc_Joint_bg= getConc(joint_ncal_params, log(MFI_BG))[,1])
      
      sample_log_conc_df<- rbind(sample_log_conc_df, sample_norm_df)
    }
  }
  
  return(sample_log_conc_df)
  
}


#Normalization


