

# Lab 1 code --------------------------------------------------------------

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

  #testing whether 2 tidy functions are equivalent
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
  
  
plate_1_orig_bg<- rm_background_orig(plate_1, method="division")
plate_1_orig_sub<- rm_background_orig(plate_1, method="subtraction")

plate_1_base_bg<- rm_background_base(plate_1, method="division")
plate_1_base_sub<- rm_background_base(plate_1, method="subtraction")
class(plate_1_base_sub)

# 2 functions are equivalent
all.equal(plate_1_orig_bg, plate_1_base_bg, check.attributes = TRUE)
all.equal(plate_1_orig_sub, plate_1_base_sub, check.attributes = TRUE)

all.equal(plate_1_orig_sub, plate_1_orig_bg, check.attributes = TRUE)
colnames(plate_1_orig_bg)
colnames(plate_1_base_bg)



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


  filter_orig_df <- filter_low_beads_orig(plate_1_orig_sub)
  filter_base_df <- filter_low_beads_base(plate_1_base_sub)

  all.equal(filter_orig_df, filter_base_df, check.attributes = TRUE)
  

# MFI transformations -----------------------------------------------------

  #orig function
  get_mfi_transforms_orig<- function(filt_df){
    
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
  
  #base r function
  get_mfi_transforms_base <- function(filt_df, antigen_names) {
    
    transformed_df <- data.frame()
    
    for (i in seq_along(antigen_names)) {
      
      ag_name <- antigen_names[i]
      
      # Filter for this antigen and valid rows
      ag_df <- filt_df[filt_df$Antigen == ag_name &
                         filt_df$Sample_Type == "TestSample" &
                         filt_df$MFI_BG >= 0, ]
      
      # Box-Cox lambda estimation
      boxcox_fit <- boxcox(MFI_BG ~ 1, plotit = FALSE, data = ag_df)
      lambda <- boxcox_fit$x[which.max(boxcox_fit$y)]
      
      # Apply transformations
      ag_df$Log_BG     <- log(ag_df$MFI_BG)
      ag_df$BoxCox_BG  <- (ag_df$MFI_BG^lambda - 1) / lambda
      ag_df$Lambda     <- lambda
      
      # Bind to result
      transformed_df <- rbind(transformed_df, ag_df)
    }
    
    # Reshape long format
    transformed_df_long <- reshape(transformed_df,
                                   varying = list(c("MFI_BG", "Log_BG", "BoxCox_BG")),
                                   v.names = "BG_Trans_Value",
                                   timevar = "Transform_Method",
                                   times = c("MFI_BG", "Log_BG", "BoxCox_BG"),
                                   direction = "long")
    
    # Add cleaned method names
    transformed_df_long$Transform_Method_Clean <- ifelse(transformed_df_long$Transform_Method == "MFI_BG", "None",
                                                         ifelse(transformed_df_long$Transform_Method == "Log_BG", "Natural Log", "Box-Cox"))
    
    # Annotate text
    annotate_df <- aggregate(Lambda ~ Antigen, data = transformed_df_long, FUN = function(x) paste0("Lambda=", round(x[1], 3)))
    annotate_df$Transform_Method_Clean <- "None"  # Just needed for merge aesthetics in ggplot
    
    # Plot
    dist_plot <- ggplot(transformed_df_long, aes(x = BG_Trans_Value, color = Transform_Method_Clean)) +
      geom_density() +
      theme_bw() +
      facet_wrap(~ Antigen, scales = "free") +
      xlab("Transformed MFI") +
      ylab("Density") +
      ggtitle(paste0("Plate ", unique(transformed_df_long$Plate))) +
      scale_color_brewer(palette = "Dark2") +
      geom_text(data = annotate_df,
                mapping = aes(x = -Inf, y = -Inf, label = Lambda),
                hjust = -0.1, vjust = -1, color = "black")
    
    return(dist_plot)
  }
  
  get_mfi_transforms_orig(filter_orig_df)
  get_mfi_transforms_base(filter_base_df)
  
  
  
  filter_orig_df <- filter_low_beads_orig(plate_1_orig_sub)
  filter_base_df <- filter_low_beads_base(plate_1_base_sub)
  
  all.equal(filter_orig_df, filter_base_df, check.attributes = TRUE)
  
  
# Lab 2 code --------------------------------------------------------------


