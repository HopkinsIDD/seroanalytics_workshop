

# Lecture 3 ---------------------------------------------------------------


# Histograms / Exploring data ---------------------------------------------

#bins <- can modify 

#check distributions of data
faceted_natural_scale <- ggplot(dat, aes(x = MFI)) +
  geom_histogram(bins = 30, color = "black", fill = "blue", alpha = 0.7) +
  facet_wrap(~ antigen, scales = "free") +
  labs(
    title = "Faceted Distribution (Natural Scale)",
    x = "Antigen MFI",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),  # Smaller font for facet labels
    axis.text = element_text(size = 10),  # Smaller font for axis text
    plot.title = element_text(size = 10)  # Smaller font for the plot title
  )



# Using standard curve to convert UI to MFI ----------------------------------


#need to convert this to base r
library(drc)  # Load package for dose-response curves

# Example dataset: standard MFI values and known IU/mL concentrations
standard_data <- data.frame(
  MFI = c(500, 1500, 5000, 15000, 45000),  # Example MFI values
  IU_mL = c(0.1, 1, 10, 100, 1000)  # Known concentrations
)

# Fit a 4PL logistic regression model
standard_curve <- drm(IU_mL ~ MFI, data = standard_data, fct = LL.4())

# Visualize the standard curve
plot(standard_curve, log = "x", xlab = "MFI", ylab = "IU/mL")



# Applying a cutoff to calculate seropositivity ---------------------------

mfi <- c(1, 5, 10)
cutoff <- 6
Seropositivity_vector <- ifelse(mfi > cutoff, 1, 0 )
#Explain how this function treats mfi variable as a vector to apply statement to each item in vector

table(Seropositivity,useNA="always")
round(prop.table(table(Seropositivity_vector,useNA="always")),3)*100

#plot cutoff on histogram

