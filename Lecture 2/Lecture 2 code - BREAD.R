
require(flexfit)
# Lab 1 code --------------------------------------------------------------

#if these are in the source file, why are they present again in the lab? 
#functions should be in source file, unless we want to specifically teach how some parts work. 

  #orig function that requires dplyr and tidr
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
  
  #adjusted function in base r
  read_and_tidy_base <- function(file_name, plate_number, num_wells, antigen_names, control_samples, background_samples, standard_curve_values, bead_threshold) {
    
    data <- read.csv(file_name, stringsAsFactors = FALSE)
    
    MFI_row <- which(data == "Median", arr.ind = TRUE)[1]
    CountVal_matrix <- which(data == "Count", arr.ind = TRUE)
    BeadCount_row <- CountVal_matrix[which(CountVal_matrix[,1] > MFI_row), 1]
    
    MFI_matrix <- data[(MFI_row + 2):(MFI_row + 1 + num_wells), 1:(length(antigen_names) + 3)]
    colnames(MFI_matrix) <- as.character(data[MFI_row + 1, 1:(length(antigen_names) + 3)])
    
    # Remove 'Total Events' column if present
    MFI_matrix <- MFI_matrix[, !colnames(MFI_matrix) %in% "Total Events"]
    
    # Convert to long format using base R
    MFI_long <- reshape(
      MFI_matrix,
      varying = antigen_names,
      v.names = "MFI",
      timevar = "Antigen",
      times = antigen_names,
      idvar = "Sample",
      direction = "long"
    )
    
    BeadCount_matrix <- data[(BeadCount_row + 2):(BeadCount_row + 1 + num_wells), 1:(length(antigen_names) + 3)]
    colnames(BeadCount_matrix) <- as.character(data[BeadCount_row + 1, 1:(length(antigen_names) + 3)])
    
    BeadCount_matrix <- BeadCount_matrix[, !colnames(BeadCount_matrix) %in% "Total Events"]
    
    BeadCount_long <- reshape(
      BeadCount_matrix,
      varying = antigen_names,
      v.names = "BeadCount",
      timevar = "Antigen",
      times = antigen_names,
      idvar = "Sample",
      direction = "long"
    )
    
    # Merge MFI and BeadCount long tables
    output_df <- merge(MFI_long, BeadCount_long, by = c("Sample", "Antigen"), all = TRUE)
    
    # Add metadata and type classification
    output_df$Plate <- plate_number
    output_df$Sample_Type <- ifelse(
      output_df$Sample %in% background_samples, "BG",
      ifelse(output_df$Sample %in% control_samples, "Ctrl",
             ifelse(output_df$Sample %in% standard_curve_values$Sample, "StdCurve", "TestSample")
      )
    )
    
    # Ensure numeric values
    output_df$MFI <- as.numeric(output_df$MFI)
    output_df$BeadCount <- as.numeric(output_df$BeadCount)
    
    # Flag low bead counts
    output_df$Low_Beads <- ifelse(output_df$BeadCount < bead_threshold, 1, 0)
    
    return(output_df)
  }

  #we don't need to show the participants the test of the read and tidy for base and tidyR
  
  #this function should be done one plate at a time, and each plate should be named as a separate object here 
  
  #testing whether 2 tidy functions are equivalent
#make sure these quantities (antigen names, control samples etc. are defined within this document not just in lab 1.)
  setwd("/Users/sarahlapidus/Documents/Pipeline docs/Data for R/magpix_output/peterson_magpix")
  
  plate_1_orig_tidy<- read_and_tidy("20230307_uvira_P1_20230304_093310.csv", plate_number = "1", num_wells=96, antigen_names= antigen_names, control_samples = control_samples, background_samples = background_samples,standard_curve_values = standard_curve_df, bead_threshold = 30) 
  head(plate_1_orig_tidy)
  class(plate_1_orig_tidy)
  
  plate_1_base_tidy <- read_and_tidy("20230307_uvira_P1_20230304_093310.csv", plate_number = "1", num_wells=96, antigen_names= antigen_names, control_samples = control_samples, background_samples = background_samples,standard_curve_values = standard_curve_df, bead_threshold = 30) 
  head(plate_1_base_tidy)
  class(plate_1_base_tidy)
  
  #testing if dfs are equivalent 
  colnames(plate_1_orig_tidy) == colnames(plate_1_base_tidy)
  all.equal(plate_1_orig_tidy, plate_1_base_tidy, check.attributes = TRUE)
  

  
# Lecture 2 code ----------------------------------------------------------



# Remove background -------------------------------------------------------
#fine to leave all the options and then say that we discuss comparisons in the lecture but we don't have them do it in lab
  
  #orig function 
  rm_background_orig<- function(clean_df, method){
    
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

  rm_background_base <- function(clean_df, method) {
    
    # Ensure Antigen is character, not factor
    clean_df$Antigen <- as.character(clean_df$Antigen)
    
    # Subset background samples
    bg_df <- clean_df[clean_df$Sample_Type == "BG", ]
    
    # Compute median background MFI per Antigen (as a named vector)
    median_values <- tapply(bg_df$MFI, bg_df$Antigen, median, na.rm = TRUE)
    
    # Match median values back to full dataset using Antigen
    Median_BG <- median_values[clean_df$Antigen]
    
    # Add Median_BG column to clean_df
    clean_df$Median_BG <- as.numeric(Median_BG)
    
    # Apply background correction
    clean_df$MFI_BG <- switch(method,
                              "division"    = clean_df$MFI / clean_df$Median_BG,
                              "subtraction" = clean_df$MFI - clean_df$Median_BG,
                              stop("Invalid method: use 'division' or 'subtraction'")
    )
    
    # Ensure MFI_BG is clean numeric
    clean_df$MFI_BG <- as.numeric(clean_df$MFI_BG)
    
    return(clean_df)
  }
  
#lets show the subtraction option again make sure to do this for each plate, and still keep each plate as its own object   
plate_1_orig_bg<- rm_background_orig(plate_1, method="division")
plate_1_orig_sub<- rm_background_orig(plate_1, method="subtraction")

plate_1_base_bg<- rm_background_base(plate_1, method="division")
plate_1_base_sub<- rm_background_base(plate_1, method="subtraction")
class(plate_1_base_sub)

#this can go. 
# 2 functions are equivalent
all.equal(plate_1_orig_bg, plate_1_base_bg, check.attributes = TRUE)
all.equal(plate_1_orig_sub, plate_1_base_sub, check.attributes = TRUE)

all.equal(plate_1_orig_sub, plate_1_orig_bg, check.attributes = TRUE)
colnames(plate_1_orig_bg)
colnames(plate_1_base_bg)

#have them explore a little bit post background correction (e.g. what are the new columns etc.) if there is subtraction have them identify if anything is <= 0 (e.g. identify the rows <=0)


# Filter low bead count ---------------------------------------------------

  #orig function
  filter_low_beads_orig<- function(bg_df){
    
    filtered_bg_df<- bg_df|>
      filter((Sample_Type!="StdCurve"&Low_Beads==0)|Sample_Type=="StdCurve")
    
    return(filtered_bg_df)
  }

  #base r function
  filter_low_beads_base <- function(bg_df) {
    
    # Logical condition: (Sample_Type != "StdCurve" & Low_Beads == 0) OR Sample_Type == "StdCurve"
    keep_rows <- (bg_df$Sample_Type != "StdCurve" & bg_df$Low_Beads == 0) |
      (bg_df$Sample_Type == "StdCurve")
    
    # Subset the rows
    filtered_bg_df <- bg_df[keep_rows, ]
    
    return(filtered_bg_df)
  }

#explore how this changes the dataframe (no. of obs, per antigen etc.) what is the definition of low beads. 
  
  filter_orig_df <- filter_low_beads_orig(plate_1_orig_sub)
  filter_base_df <- filter_low_beads_base(plate_1_base_sub)

  all.equal(filter_orig_df, filter_base_df, check.attributes = TRUE)
  

##SOPHIE:: to add the BREAD function for standard curve visualization (with comparison to MFI values of samples on plates)
 
  
  # Standardization with FlexFit-----------------------------------------------------
 
  #this is pulled from FlexFit sourcecode so it can go in the source file  
  FUNinv <- function(y, par) {
    -par["Scale"]*log(((par["Aup"] - par["Alow"])/
                         (y - par["Alow"]))^(1/par["a"]) - 1) + par["Xmid"]
  }
#this is also done on a plate by plate basis so each plate needs to be handled separately the 'logging" is done within this function. 
  #since this function has different options depending on the input (raw, bg subtracted) we can have them run it once per plate on the background subtracted vs raw values for all agxs 
  #and then we can have them compare how many estimates become "NA" under each circumstance (i.e. is the curve unable to be fit because the shape is not right or is the estiamte beyond the upper/lower asymptote)
  get_concentration<- function(plate_df_norm, std_curve_values, input){
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
      
      
      
      if(input=="divMFI"){
        
        antigen_curve_complete<- std_curve_df|>
          filter(Antigen==antigen_list[i])|>
          mutate(assay_id= Plate, analyte=Antigen, expected_conc= Dilution, fi= MFI_bg,
                 well_role="Standard", sample_id=Sample, replicate= Replicate,
                 Log_Dilution= log(Dilution), Log_MFI_bg_val= log(MFI_bg))|>
          dplyr::select(assay_id, analyte, expected_conc, Dilution, fi, well_role, sample_id, 
                        replicate,Log_Dilution,Log_MFI_bg_val) 
        
        
        
        std_curveFit_antigen<- flexfit::fitStd( std= antigen_curve_complete, 
                                                xvar="Log_Dilution", yvar="Log_MFI_bg_val", interactive = F)
        
        
        
        
        if(is.null(std_curveFit_antigen$par)==F){
          sample_norm_df<- plate_df_norm|>
            filter(Sample_Type%in%c("Ctrl", "TestSample"), Antigen==antigen_list[i])|>
            mutate(Log_MFI_bg_val= log(MFI_bg),
                   Log_Conc_bg=FUNinv(log(MFI_bg),std_curveFit_antigen$par))
          
          sample_log_conc_df<- rbind(sample_log_conc_df, sample_norm_df)
          
         
        }
        else{
          sample_norm_df<- plate_df_norm|>
            filter(Sample_Type%in%c("Ctrl", "TestSample"), Antigen==antigen_list[i])|>
            mutate(Log_MFI_bg_val= log(MFI_bg), 
                  Log_Conc_bg=NA)
          
          sample_log_conc_df<- rbind(sample_log_conc_df, sample_norm_df)
          
         
        }
        
      }
      
      
      if(input=="subMFI"){
        
        antigen_curve_complete<- std_curve_df|>
          filter(Antigen==antigen_list[i])|>
          mutate(assay_id= Plate, analyte=Antigen, expected_conc= Dilution, fi= MFI_bg,
                 well_role="Standard", sample_id=Sample, replicate= Replicate,
                 Log_Dilution= log(Dilution), Log_MFI_bg_val= log(MFI_bg), )|>
          dplyr::select(assay_id, analyte, expected_conc, Dilution, fi, well_role, sample_id, 
                        replicate,Log_Dilution,Log_MFI_bg_val) 
        
        
        
        std_curveFit_antigen<- flexfit::fitStd( std= antigen_curve_complete, 
                                                xvar="Log_Dilution", yvar="Log_MFI_bg_val", interactive = F)
        
       
        
        
        if(is.null(std_curveFit_antigen$par)==F){
          sample_norm_df<- plate_df_norm|>
            filter(Sample_Type%in%c("Ctrl", "TestSample"), Antigen==antigen_list[i])|>
            mutate(Log_MFI_bg_val= log(MFI_bg), 
                   
                   Log_Conc_bg=FUNinv(log(MFI_bg),std_curveFit_antigen$par))
          
          sample_log_conc_df<- rbind(sample_log_conc_df, sample_norm_df)
         
        }
        else{
          sample_norm_df<- plate_df_norm|>
            filter(Sample_Type%in%c("Ctrl", "TestSample"), Antigen==antigen_list[i])|>
            mutate(Log_MFI_bg_val= log(MFI_bg), 
                   
                   Log_Conc_bg=NA)
          
          sample_log_conc_df<- rbind(sample_log_conc_df, sample_norm_df)
          
         
        }
        
      }

    }
    
    
    
    
    
    return(sample_log_conc_df)
    
  }

#this standardization also needs to be run on a plate by plate basis. 
  

  # Normalization with Linear Model-----------------------------------------------------
  

#this requires you to put all the plates together into a long data frame (once you've worked through BREAD up to this point let me know and i will show an example of how to do it)
#similar to above, we can have them explore with raw input, or subtracted input 
  

  get_norm_df<- function(long_df, method){
    
   
    
    
  
    if(method=="MFI"){
      
      long_df<- long_df|>
        mutate(Log_MFI= log(MFI))
      
      norm_df_full=data.frame()
      norm_methods= c("Log_MFI", "Log_Conc_bg")
      for(i in 1:length(norm_methods)){
        std_curve_fit= norm_methods[i]
        wide_df<- long_df|>
          dplyr::select( Sample, Antigen, std_curve_fit, Plate, Sample_Type)|>
          pivot_wider(id_cols=c('Plate', 'Sample', 'Sample_Type'), names_from= Antigen, values_from = std_curve_fit)|>
          drop_na(any_of(antigen_names))|> 
          filter_all(all_vars(!is.infinite(.)))
        
        
        ctrl_df<- long_df|>
          dplyr::select(Antigen, Plate, Sample, Sample_Type, std_curve_fit)|>
          filter(Sample_Type%in%c("StdCurve", "Ctrl"))|>
          drop_na(any_of(antigen_names))|>
          filter_all(all_vars(!is.infinite(.)))
        
        n_obs_per_plate<- long_df|>
          group_by(Plate)|>
          summarize(n_obs=n())
        #simple LM with ctrls
        
        lm_ctrls<- lm(ctrl_df[,std_curve_fit]~Antigen+Plate+Sample, data= ctrl_df)
        plate_coefs<- lm_ctrls$coefficients[startsWith(names(lm_ctrls$coefficients),"Plate")]
        
        long_df<- long_df|>
          mutate(plate_effect= rep(c(0, plate_coefs), c(n_obs_per_plate$n_obs)))|>
          mutate(LM_norm_MFI= long_df[,std_curve_fit]-plate_effect)
        
        norm_df<-long_df|>
          #left_join(ma_weightedReg_long, by=c('Plate', 'Sample', 'Location', 'Sample_Type', 'Antigen'))|>
          #left_join(ma_Reg_long,by=c('Plate', 'Sample', 'Sample_Type', 'Antigen'))|>
          #left_join(trim_df, by=c('Plate', 'Sample', 'Location', 'Sample_Type', 'Antigen') )|>
          #left_join(long_df, by= c('Plate', 'Sample', 'Sample_Type', 'Antigen'))|>
          dplyr::select(Plate, Sample,  Sample_Type, Antigen, 
                        std_curve_fit, LM_norm_MFI)|>
          pivot_longer(cols= !c('Plate', 'Sample', 'Sample_Type', 'Antigen'), names_to="Norm_Method", values_to="Norm_MFI")|>
          mutate(Input_Value = std_curve_fit)
        
        norm_df_full<- rbind(norm_df_full, norm_df)
      }
      
    }
    
    if(method=="bgMFI"){
      
      long_df<- long_df|>
        mutate(Log_MFI_BG= log(MFI_bg))
      
      norm_df_full=data.frame()
      norm_methods= c("Log_MFI_BG", "Log_Conc_bg")
      for(i in 1:length(norm_methods)){
        std_curve_fit= norm_methods[i]
        wide_df<- long_df|>
          dplyr::select(Sample, Antigen, std_curve_fit, Plate, Sample_Type)|>
          pivot_wider(id_cols=c('Plate', 'Sample',  'Sample_Type'), names_from= Antigen, values_from = std_curve_fit)|>
          drop_na(any_of(antigen_names))|>
          filter_all(all_vars(!is.infinite(.)))
        
        #MA loess
       
        ctrl_df<- long_df|>
          dplyr::select(Antigen, Plate, Sample, Sample_Type, std_curve_fit)|>
          filter(Sample_Type%in%c("StdCurve", "Ctrl"))|>
          drop_na(any_of(antigen_names))|>
          filter_all(all_vars(!is.infinite(.)))
        
        n_obs_per_plate<- long_df|>
          group_by(Plate)|>
          summarize(n_obs=n())
        #simple LM with ctrls
        
        lm_ctrls<- lm(ctrl_df[,std_curve_fit]~Antigen+Plate+Sample, data= ctrl_df)
        plate_coefs<- lm_ctrls$coefficients[startsWith(names(lm_ctrls$coefficients),"Plate")]
        
        long_df<- long_df|>
          mutate(plate_effect= rep(c(0, plate_coefs), c(n_obs_per_plate$n_obs)))|>
          mutate(LM_norm_MFI= long_df[,std_curve_fit]-plate_effect)
        
        norm_df<-long_df|>
          #left_join(ma_weightedReg_long, by=c('Plate', 'Sample', 'Location', 'Sample_Type', 'Antigen'))|>
          #full_join(ma_Reg_long,by=c('Plate', 'Sample',  'Sample_Type', 'Antigen'))|>
          #full_join(long_df, by= c('Plate', 'Sample',  'Sample_Type', 'Antigen'))|>
          dplyr::select(Plate, Sample,  Sample_Type, Antigen, std_curve_fit, LM_norm_MFI)|>
          pivot_longer(cols= !c('Plate', 'Sample',  'Sample_Type', 'Antigen'), names_to="Norm_Method", values_to="Norm_MFI")|>
          mutate(Input_Value = std_curve_fit)
        norm_df_full<- rbind(norm_df_full, norm_df)
      }
      
    }
    
    
    if(method=="subMFI"){
      
      
      long_df<- long_df|>
        mutate(Log_MFI_BG_corrected= log(MFI_bg))
      
      norm_df_full=data.frame()
      norm_methods= c("Log_MFI_BG_corrected")
      for(i in 1:length(norm_methods)){
        std_curve_fit= norm_methods[i]
        wide_df<- long_df|>
          dplyr::select( Sample, Antigen, std_curve_fit, Plate, Sample_Type)|>
          pivot_wider(id_cols=c('Plate', 'Sample', 'Sample_Type'), names_from= Antigen, values_from = std_curve_fit)|>
          drop_na(any_of(antigen_names))|>
          filter_all(all_vars(!is.infinite(.)))
        
       
        ctrl_df<- long_df|>
          dplyr::select(Antigen, Plate, Sample, Sample_Type, std_curve_fit)|>
          filter(Sample_Type%in%c("StdCurve", "Ctrl"))|>
          drop_na(any_of(antigen_names))|>
          filter_all(all_vars(!is.infinite(.)))
        
        n_obs_per_plate<- long_df|>
          group_by(Plate)|>
          summarize(n_obs=n())
        #simple LM with ctrls
        
        lm_ctrls<- lm(ctrl_df[,std_curve_fit]~Antigen+Plate+Sample, data= ctrl_df)
        plate_coefs<- lm_ctrls$coefficients[startsWith(names(lm_ctrls$coefficients),"Plate")]
        
        long_df<- long_df|>
          mutate(plate_effect= rep(c(0, plate_coefs), c(n_obs_per_plate$n_obs)))|>
          mutate(LM_norm_MFI= long_df[,std_curve_fit]-plate_effect)
        
        norm_df<-long_df|>
          #left_join(ma_weightedReg_long, by=c('Plate', 'Sample', 'Location', 'Sample_Type', 'Antigen'))|>
          #full_join(ma_Reg_long,by=c('Plate', 'Sample',  'Sample_Type', 'Antigen'))|>
          #left_join(trim_df, by=c('Plate', 'Sample', 'Location', 'Sample_Type', 'Antigen') )|>
          #full_join(long_df, by= c('Plate', 'Sample',  'Sample_Type', 'Antigen'))|>
          dplyr::select(Plate, Sample,  Sample_Type, Antigen, std_curve_fit, LM_norm_MFI)|>
          pivot_longer(cols= !c('Plate', 'Sample',  'Sample_Type', 'Antigen'), names_to="Norm_Method", values_to="Norm_MFI")|>
          mutate(Input_Value = std_curve_fit)
        norm_df_full<- rbind(norm_df_full, norm_df)
      }
      
    }
    return(norm_df_full)
  }
  
  #have them explore different distributions, size of data, number of NAs etc. 
  