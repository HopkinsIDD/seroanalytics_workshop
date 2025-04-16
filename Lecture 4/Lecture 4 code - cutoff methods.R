
# Lecture 4 ---------------------------------------------------------------



# Mean + SD code ----------------------------------------------------------

#go through this with 
neg_controls <- dat$AntibodyLevel[dat$KnownStatus == "Negative"]
cutoff_2sd <- mean(neg_controls, na.rm = TRUE) + 2 * sd(neg_controls, na.rm = TRUE)
cutoff_3sd <- mean(neg_controls, na.rm = TRUE) + 3 * sd(neg_controls, na.rm = TRUE)
print(cutoff_2sd)
print(cutoff_3sd)


# mean_sd_function -------------------------------------------------------------------------

mean_sd_function <- function(antigen_mfi, pos_neg_status) {
  
  # Check for negative controls
  if (all(!grepl("^(0|neg|negative)$", pos_neg_status, ignore.case = TRUE))) {
    # If no negatives, return a formatted data frame
    return(data.frame(
      Method = c("Mean+2SD", "Mean+3SD"),
      Cutoff = c("Not calculated - no negative controls", "Not calculated - no negative controls"),
      Seropositivity = c("Not calculated - no negative controls", "Not calculated - no negative controls"),
      stringsAsFactors = FALSE
    ))
  }
  
  # Extract negatives
  neg_controls <- antigen_mfi[which(pos_neg_status == "Neg")]
  
  # Ensure no NA values in negatives
  neg_controls <- na.omit(neg_controls)
  
  # Calculate cutoff values
  cut_mean_plus2sd <- round(mean(neg_controls) + 2 * sd(neg_controls),2)
  cut_mean_plus3sd <- round(mean(neg_controls) + 3 * sd(neg_controls),2)
  
  # Calculate seropositivity percentages
  seropos_mean_plus2sd <- mean(antigen_mfi > cut_mean_plus2sd, na.rm = TRUE) * 100
  seropos_mean_plus3sd <- mean(antigen_mfi > cut_mean_plus3sd, na.rm = TRUE) * 100
  
  # Create the output data frame
  results <- data.frame(
    Method = c("Mean+2SD", "Mean+3SD"),
    Cutoff = c(cut_mean_plus2sd, cut_mean_plus3sd),
    Seropositivity = round(c(seropos_mean_plus2sd, seropos_mean_plus3sd), 1),
    stringsAsFactors = FALSE
  )
  
  return(results)
}


mean_sd_all_antigens <- function(control_dat) {
  
  antigens <- unique(control_dat$antigen)
  num_antigens <- length(antigens)
  
  # Initialize an empty data frame to store results
  all_results <- data.frame()
  
  # Initialize a data frame for plotting
  plot_data <- data.frame()
  
  # Loop through each antigen
  for (i in seq_len(num_antigens)) {
    
    antigen <- unique(control_dat$antigen[which(control_dat$antigen == antigens[i])])
    control_antigen_df <- control_dat[which(control_dat$antigen == antigens[i]),]
    
    # Extract antigen_mfi and pos_neg_status for the current antigen
    antigen_mfi <- control_antigen_df$MFI
    pos_neg_status <- control_antigen_df$pos_neg
    
    # Apply the cutoff functions
    combined_df <- mean_sd_function(antigen_mfi, pos_neg_status)
    
    # Add the antigen name to the data frame
    
    combined_df_antigen_name <- data.frame(
      antigen = antigen,  # Add the antigen name as the first column
      combined_df  # Append the rest of the data frame
    )
    # Combine with all_results
    all_results <- rbind(all_results, combined_df_antigen_name)
    
    # Add data for plotting
    plot_data <- rbind(plot_data, data.frame(
      antigen = antigen,
      antigen_mfi = antigen_mfi
    ))
  }
  
  # Generate the facetted plot
  facetted_plot <- ggplot(plot_data, aes(x = antigen_mfi)) +
    geom_histogram(bins = 50, fill = "blue", alpha = 0.6, color = "black") +
    geom_vline(
      data = all_results,
      aes(xintercept = as.numeric(Cutoff), linetype = Method),
      color = "red",
      linewidth = 0.5
    ) +
    scale_linetype_manual(values = c("dashed", "dashed")) +  # Customize line types
    facet_wrap(~ antigen, scales = "free", ncol = 4) +
    labs(
      title = "Distribution of Antigen MFI with Cutoffs",
      x = "Antigen MFI",
      y = "Frequency",
      linetype = "Cutoff Method"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Return the final combined data frame and the facetted plot
  return(list(
    Results = all_results,
    Distribution_Plot = facetted_plot
  ))
}

# Example usage
mean_plus_sd_results <- mean_sd_all_antigens(control_dat)


# Calculating sensitivity and specificity to show in class ---------------------------------


# Define cutoff value
cutoff_value <- 500

# Create Seropos variable based on MFI cutoff
dat$Seropos <- ifelse(dat$MFI > cutoff_value, 1, 0)

# Calculate confusion matrix values
TP <- sum(dat$Seropos == 1 & dat$Pos_control == 1)  # True Positives
TN <- sum(dat$Seropos == 0 & dat$Pos_control == 0)  # True Negatives
FP <- sum(dat$Seropos == 1 & dat$Pos_control == 0)  # False Positives
FN <- sum(dat$Seropos == 0 & dat$Pos_control == 1)  # False Negatives

# Calculate Sensitivity and Specificity
Sensitivity <- TP / (TP + FN)
Specificity <- TN / (TN + FP)

# Print results
cat("True Positives:", TP, "\n")
cat("True Negatives:", TN, "\n")
cat("False Positives:", FP, "\n")
cat("False Negatives:", FN, "\n")
cat("Sensitivity:", round(Sensitivity, 3), "\n")
cat("Specificity:", round(Specificity, 3), "\n")

# ROC curve to show in class ---------------------------------------------------------------

#code using pROC package
# library(pROC) 
# roc_obj <- roc(response = dat$KnownStatus, predictor = dat$AntibodyLevel)
# plot(roc_obj, print.auc = TRUE)  # Plot ROC curve with AUC
#coords(roc_obj, "best", ret = "threshold")  # Get optimal cutoff <- redo in base R

#code redone in base R
# Sample data (for demonstration)
set.seed(123)
dat <- data.frame(
  KnownStatus = c(rep(1, 50), rep(0, 50)),
  AntibodyLevel = c(rnorm(50, mean = 1.5, sd = 0.5), rnorm(50, mean = 0.5, sd = 0.5))
)

# Step 1: Create thresholds
thresholds <- sort(unique(dat$AntibodyLevel))

# Step 2: Compute TPR (sensitivity) and FPR (1 - specificity) for each threshold
roc_df <- data.frame(threshold = thresholds, TPR = NA, FPR = NA)

for (i in seq_along(thresholds)) {
  t <- thresholds[i]
  predicted <- ifelse(dat$AntibodyLevel >= t, 1, 0)
  
  TP <- sum(predicted == 1 & dat$KnownStatus == 1)
  FP <- sum(predicted == 1 & dat$KnownStatus == 0)
  FN <- sum(predicted == 0 & dat$KnownStatus == 1)
  TN <- sum(predicted == 0 & dat$KnownStatus == 0)
  
  TPR <- TP / (TP + FN)
  FPR <- FP / (FP + TN)
  
  roc_df$TPR[i] <- TPR
  roc_df$FPR[i] <- FPR
}

# Step 3: Calculate AUC using the trapezoidal rule
roc_df <- roc_df[order(roc_df$FPR), ]
auc <- sum(diff(roc_df$FPR) * (head(roc_df$TPR, -1) + tail(roc_df$TPR, -1)) / 2)

# Step 4: Plot ROC curve using ggplot2
library(ggplot2)

ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(color = "blue", size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(title = paste0("ROC Curve (AUC = ", round(auc, 3), ")"),
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme_minimal()

# Calculate Youden's Index
roc_df$Youden <- roc_df$TPR - roc_df$FPR

# Get the threshold with the highest Youden index
best_row <- roc_df[which.max(roc_df$Youden), ]
best_threshold <- best_row$threshold

# Print result
cat("Best threshold (max Youden Index):", best_threshold, "\n")




# ROC functions for students to use ---------------------------------------


#function to check for pos and neg controls
checking_pos_neg_controls <- function(pos_neg_status) {
  
  # Check if both positive and negative controls are present
  has_positives <- any(grepl("^(1|pos|positive)$", pos_neg_status, ignore.case = TRUE))
  has_negatives <- any(grepl("^(0|neg|negative)$", pos_neg_status, ignore.case = TRUE))
  
  if (!has_positives && !has_negatives) {
    # If neither positives nor negatives are present
    return(data.frame(
      Method = c("ROC_YoudenJ", "ROC_CDC_Likelihood"),
      Cutoff = c("Unable to calculate (No positive or negative controls)", 
                 "Unable to calculate (No positive or negative controls)"),
      Seropositivity = c("Unable to calculate (No positive or negative controls)", 
                         "Unable to calculate (No positive or negative controls)"),
      Sensitivity = c("Unable to calculate (No positive or negative controls)", 
                      "Unable to calculate (No positive or negative controls)"),
      Specificity = c("Unable to calculate (No positive or negative controls)", 
                      "Unable to calculate (No positive or negative controls)"),
      stringsAsFactors = FALSE
    ))
  } else if (!has_positives) {
    # If no positives are present
    return(data.frame(
      Method = c("ROC_YoudenJ", "ROC_CDC_Likelihood"),
      Cutoff = c("Unable to calculate (No positive controls)", 
                 "Unable to calculate (No positive controls)"),
      Seropositivity = c("Unable to calculate (No positive controls)", 
                         "Unable to calculate (No positive controls)"),
      Sensitivity = c("Unable to calculate (No positive controls)", 
                      "Unable to calculate (No positive controls)"),
      Specificity = c("Unable to calculate (No positive controls)", 
                      "Unable to calculate (No positive controls)"),
      stringsAsFactors = FALSE
    ))
  } else if (!has_negatives) {
    # If no negatives are present
    return(data.frame(
      Method = c("ROC_YoudenJ", "ROC_CDC_Likelihood"),
      Cutoff = c("Unable to calculate (No negative controls)", 
                 "Unable to calculate (No negative controls)"),
      Seropositivity = c("Unable to calculate (No negative controls)", 
                         "Unable to calculate (No negative controls)"),
      Sensitivity = c("Unable to calculate (No negative controls)", 
                      "Unable to calculate (No negative controls)"),
      Specificity = c("Unable to calculate (No negative controls)", 
                      "Unable to calculate (No negative controls)"),
      stringsAsFactors = FALSE
    ))
  }
  # If both positives and negatives are present, proceed with function logic
}

roc_function <- function(antigen_mfi, pos_neg_status, min_sensitivity = NULL) {
  # Ensure the pROC package is installed
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("The 'pROC' package is required but is not installed. Please install it using install.packages('pROC').")
  }
  
  # Load necessary libraries
  library(pROC)
  library(ggplot2)
  
  # Check for positive and negative controls
  check_result <- checking_pos_neg_controls(pos_neg_status)
  
  # If controls are not adequate, return the check_result directly
  if (!is.null(check_result)) {
    return(list(Results = check_result, ROC_Plot = NULL))
  }
  
  # If min_sensitivity is provided, validate and scale it
  if (!is.null(min_sensitivity)) {
    if (!is.numeric(min_sensitivity) || min_sensitivity < 0 || min_sensitivity > 100) {
      stop("min_sensitivity must be a numeric value between 0 and 100.")
    }
    if (min_sensitivity > 1) {
      min_sensitivity <- min_sensitivity / 100
    }
  } else {
    min_sensitivity <- 0.70 # Default to 70% if not specified
  }
  
  # Generate the ROC curve
  roc_curve <- roc(pos_neg_status, antigen_mfi)
  
  # Youden's J Index: Find the cutoff
  youden_index <- which.max(roc_curve$sensitivities + roc_curve$specificities - 1)
  youden_cutoff <- roc_curve$thresholds[youden_index]
  youden_sensitivity <- roc_curve$sensitivities[youden_index]
  youden_specificity <- roc_curve$specificities[youden_index]
  youden_seropositivity <- sum(antigen_mfi > youden_cutoff, na.rm = TRUE) / length(na.omit(antigen_mfi)) * 100
  
  # Maximize specificity while keeping sensitivity above the user-defined threshold
  max_specificity_cutoff <- NA
  max_specificity_sensitivity <- NA
  max_specificity_specificity <- NA
  max_specificity_value <- -Inf
  
  for (i in 1:length(roc_curve$sensitivities)) {
    if (roc_curve$sensitivities[i] >= min_sensitivity) {
      if (roc_curve$specificities[i] > max_specificity_value) {
        max_specificity_value <- roc_curve$specificities[i]
        max_specificity_cutoff <- roc_curve$thresholds[i]
        max_specificity_sensitivity <- roc_curve$sensitivities[i]
        max_specificity_specificity <- roc_curve$specificities[i]
      }
    }
  }
  
  # Calculate seropositivity for max specificity cutoff
  if (!is.na(max_specificity_cutoff)) {
    max_specificity_seropositivity <- sum(antigen_mfi > max_specificity_cutoff, na.rm = TRUE) / length(na.omit(antigen_mfi)) * 100
  } else {
    max_specificity_seropositivity <- NA
  }
  
  # Find the lowest threshold where specificity is exactly 100%
  specificity_100_cutoff <- NA
  specificity_100_sensitivity <- NA
  
  # Find all indices where specificity is 100%
  specificity_100_indices <- which(roc_curve$specificities == 1)
  
  if (length(specificity_100_indices) > 0) {
    # Take the MINIMUM threshold that gives 100% specificity
    specificity_100_cutoff <- min(roc_curve$thresholds[specificity_100_indices], na.rm = TRUE)
    specificity_100_sensitivity <- roc_curve$sensitivities[specificity_100_indices][which.min(roc_curve$thresholds[specificity_100_indices])]
  } 
  
  # Calculate seropositivity for the 100% specificity cutoff
  if (!is.na(specificity_100_cutoff)) {
    specificity_100_seropositivity <- sum(antigen_mfi > specificity_100_cutoff, na.rm = TRUE) / length(na.omit(antigen_mfi)) * 100
  } else {
    specificity_100_seropositivity <- NA
  }
  
  # # Find the highest sensitivity cutoff where specificity is exactly 100%
  # specificity_100_cutoff <- NA
  # specificity_100_sensitivity <- NA
  # 
  # for (i in 1:length(roc_curve$specificities)) {
  #   if (roc_curve$specificities[i] == 1) {  # Specificity = 100%
  #     specificity_100_cutoff <- roc_curve$thresholds[i]
  #     specificity_100_sensitivity <- roc_curve$sensitivities[i]
  #   }
  # }
  # 
  # # Calculate seropositivity for the 100% specificity cutoff
  # if (!is.na(specificity_100_cutoff)) {
  #   specificity_100_seropositivity <- sum(antigen_mfi > specificity_100_cutoff, na.rm = TRUE) / length(na.omit(antigen_mfi)) * 100
  # } else {
  #   specificity_100_seropositivity <- NA
  # }
  
  # Create a results data frame
  results <- data.frame(
    Method = c(
      "ROC_YoudenJ",
      paste0("Max_Specificity_with_Sensitivity_", round(min_sensitivity * 100), "%"),
      "100% Specificity"
    ),
    Cutoff = c(youden_cutoff, max_specificity_cutoff, specificity_100_cutoff),
    Seropositivity = c(round(youden_seropositivity, 1), round(max_specificity_seropositivity, 1), round(specificity_100_seropositivity, 1)),
    Sensitivity = c(youden_sensitivity, max_specificity_sensitivity, specificity_100_sensitivity),
    Specificity = c(youden_specificity, max_specificity_specificity, 1),
    stringsAsFactors = FALSE
  )
  
  # Return the results and the ROC plot
  return(list(
    Results = results,
    ROC_curve = roc_curve
  ))
}

# Function to create a ggplot for an ROC curve
plot_roc_ggplot <- function(roc_object, title = "ROC Curve") {
  # Extract necessary data from the roc object
  roc_data <- data.frame(
    FPR = 1 - roc_object$specificities,  # False Positive Rate
    TPR = roc_object$sensitivities       # True Positive Rate
  )
}

#function to run roc_function for all antigens
roc_all_antigens <- function(control_dat, min_sensitivity=NULL) {
  # Ensure required packages are installed
  if (!requireNamespace("pROC", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) {
    stop("The 'pROC' and 'patchwork' packages are required but are not installed. Please install them using install.packages('pROC') and install.packages('patchwork').")
  }
  
  # Load necessary libraries
  library(pROC)
  library(ggplot2)
  library(patchwork)
  
  antigens <- unique(control_dat$antigen)
  num_antigens <- length(antigens)
  
  
  # Define the expected column names for all_results
  expected_columns <- c("antigen_name", "Method", "Cutoff", "Seropositivity", "Sensitivity", "Specificity")
  
  # Initialize an empty data frame for results with the expected structure
  all_results <- data.frame(
    antigen_name = character(),
    Method = character(),
    Cutoff = numeric(),
    Seropositivity = numeric(),
    Sensitivity = numeric(),
    Specificity = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Initialize a list for ROC plots
  plot_list <- list()
  
  # Loop through each antigen
  for (i in seq_along(antigens)) {
    # Extract antigen_mfi and pos_neg_status for the current antigen
    
    antigen <- unique(control_dat$antigen[which(control_dat$antigen == antigens[i])])
    control_antigen_df <- control_dat[which(control_dat$antigen == antigens[i]),]
    
    # Extract antigen_mfi and pos_neg_status for the current antigen
    antigen_mfi <- control_antigen_df$MFI
    pos_neg_status <- control_antigen_df$pos_neg
    
    # Check if pos_neg_status has fewer than two levels
    if (length(unique(pos_neg_status)) < 2) {
      # Add a standardized error row
      error_row <- setNames(
        as.list(rep("Not calculated - missing pos or neg controls", length(expected_columns))),
        expected_columns
      )
      error_row$antigen_name <- antigen
      error_row$Method <- "ROC"
      all_results <- rbind(all_results, as.data.frame(error_row, stringsAsFactors = FALSE))
      
      # Add an error plot to the plot list
      error_plot <- ggplot() +
        geom_text(aes(x = 0.5, y = 0.5, label = "No ROC Curve\nInsufficient Controls"), size = 3, color = "red") +
        theme_void() +
        ggtitle(antigen)
      plot_list[[i]] <- error_plot
      next
    } 
    
    # Ensure there are no NA values in antigen_mfi
    antigen_mfi <- na.omit(antigen_mfi)
    
    # Skip if no valid data
    if (length(antigen_mfi) == 0) {
      
      #warning(paste("Skipping", antigen_name, "- no valid antigen_mfi values"))
      # Add a standardized error row
      error_row <- setNames(
        as.list(rep("Not calculated - no valid MFI values", length(expected_columns))),
        expected_columns
      )
      error_row$antigen_name <- antigen_name
      error_row$Method <- "ROC"
      all_results <- rbind(all_results, as.data.frame(error_row, stringsAsFactors = FALSE))
      
      # Add an error plot to the plot list
      error_plot <- ggplot() +
        geom_text(aes(x = 0.5, y = 0.5, label = "No ROC Curve\nInvalid MFI Values"), size = 3, color = "red") +
        theme_void() +
        ggtitle(antigen_name)
      plot_list[[i]] <- error_plot
      next
    }
    
    # Apply the roc_function to each antigen
    antigen_results <- tryCatch(
      roc_function(antigen_mfi, pos_neg_status, min_sensitivity),
      error = function(e) {
        list(
          Results = data.frame(
            #antigen_name = antigen_name,
            Method = "Error",
            Cutoff = NA,
            Seropositivity = NA,
            Sensitivity = NA,
            Specificity = NA,
            stringsAsFactors = FALSE
          ),
          ROC_curve = ggplot() + 
            geom_text(aes(x = 0.5, y = 0.5, label = "Error in ROC Calculation"), size = 3, color = "red") +
            theme_void() +
            ggtitle(antigen)
        )
      }
    )
    
    # Extract results and ensure consistency with expected columns
    if ("Results" %in% names(antigen_results)) {
      roc_results <- antigen_results$Results
      
      # Add antigen name
      roc_results$antigen_name <- antigen
      
      # Ensure all expected columns are present
      for (col in expected_columns) {
        if (!col %in% names(roc_results)) {
          roc_results[[col]] <- NA
        }
      }
      
      # Subset to expected columns and bind to all_results
      roc_results <- roc_results[expected_columns]
      all_results <- rbind(all_results, roc_results)
    }
    
    # Extract the ROC curve plot
    if ("ROC_curve" %in% names(antigen_results)) {
      roc_object <- antigen_results$ROC_curve
      roc_plot_dat <- plot_roc_ggplot(roc_object, title = paste("ROC Curve for", antigen))
      roc_plot <- ggplot(roc_plot_dat, aes(x = FPR, y = TPR)) +
        geom_line(color = "blue", linewidth = 1) +  # Add the ROC curve
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +  # Add the diagonal reference line
        labs(
          title = paste("ROC Curve for", antigen),
          x = "False Positive Rate (1 - Specificity)",
          y = "True Positive Rate (Sensitivity)"
        ) +
        theme_minimal()  # Apply a clean theme
      plot_list[[i]] <- roc_plot
    }
  }
  
  # Determine the grid layout dimensions
  num_antigens <- length(plot_list)
  max_columns <- 4  # Maximum number of columns
  grid_cols <- min(max_columns, ceiling(sqrt(num_antigens)))
  grid_rows <- ceiling(num_antigens / grid_cols)
  
  # Fill empty slots with blank plots
  total_slots <- grid_rows * grid_cols
  while (length(plot_list) < total_slots) {
    blank_plot <- ggplot() + theme_void() + ggtitle("")
    plot_list <- c(plot_list, list(blank_plot))
  }
  
  # Combine plots into a grid
  combined_plot <- wrap_plots(plot_list, ncol = grid_cols, nrow = grid_rows)
  
  # Return results as a list
  return(list(
    Results = all_results,
    Combined_ROC_Plot = combined_plot
  ))
}

roc_results <- roc_all_antigens(control_dat, 60)




# final histograms of data with histograms plotted ------------------------

#input should be long dataframe each antigen and cutoff methods in its own row

#get log cutoff and eliminate non-numeric values
final_df_clean <- final_df
final_df_clean$Cutoff <- as.numeric(final_df_clean$Cutoff)
final_df_clean <- final_df_clean[!is.na(final_df_clean$Cutoff), ]
final_df_clean$Log_cutoff <- log(final_df_clean$Cutoff)
final_df_clean <- final_df_clean[!is.na(final_df_clean$Log_cutoff), ]

# Create a faceted plot using all antigen data
final_plot <- ggplot() +
  geom_histogram(
    data = dat_long,
    aes(x = log(Response), y = ..density..),
    bins = 30, alpha = 0.4, fill = "gray", color = "gray",
    inherit.aes = FALSE,
    na.rm = TRUE
  ) +
  # Vertical cutoff line (if defined)
  geom_vline(
    data = final_df_clean,
    aes(xintercept = Log_cutoff, color = Method),
    linetype = "longdash",
    na.rm = TRUE
  ) +
  scale_color_brewer(palette = "Set1", name = "Cutoff Method") + # Define color scale
  facet_wrap(~ Antigen, scales = "free", ncol = 3) +
  labs(
    title = "Comparison of cutoffs",
    x = "Log Antigen MFI",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),  # Smaller font for facet labels
    axis.text = element_text(size = 6),  # Smaller font for axis text
    plot.title = element_text(size = 10)  # Smaller font for the plot title
  )        



