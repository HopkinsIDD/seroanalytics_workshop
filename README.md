**Guide to the Seroanalytic Workshop files**

The following material is a training resource for learning seroanalytical methods, from understanding and pre-processing serological data to inferring transmission dynamics. The material consists of both lectures and labs, and can be followed at your own pace. To execute the code in the labs, you will need to have R, and preferably RStudio, installed on your computer before starting.

We recommend going through the files listed in the order below. We also recommend downloading the entire repository to your computer and creating a directory with the same name *seroanalytics_workshop*, as some of the R scripts will call data and code files assuming they are saved in similarly named folders. Please find a description of the course material files below.

**Part 1: Introduction to the Seroanalytics Workshop**
- Lecture 1 (pptx and pdf versions)
  - A course overview, use cases of serology, study designs in seroepidemiology, and laboratory assays, including multiplexed methods

**Part 2: Introduction to serological data analyses in R**
- Lecture 2 (pptx and pdf versions)
  - A description of R vs. Rstudio, R script vs. R markdown, ‘wide’ vs. ‘long’ dataframes, and different types of data frames used in serology
- Lab 2a_Intro to R (pdf and Rmd versions)
  - R code to practice useful functions that will be used in subsequent labs
- Lab 2b_Reading in Serological Data (pdf and Rmd versions)
  - R code to understand the format of raw multiplexed serological data

**Part 3: Pre-processing serological data**
- Lecture 3 (pptx and pdf versions)
  - A description of what is pre-processing, how to identify whether pre-processing has been effective, and a pipeline for pre-processing data
- Lab 3 (pdf and Rmd versions)
  - R code to implement a pre-processing pipeline

**Part 4: Visualizing and standardizing serological data**
- Lecture 4 (pptx and pdf versions)
  - An overview of how to visualize & standardize serological data
- Lab 4  (pdf and Rmd versions)
  - R code to visualize and standardize data

**Part 5: Determining serostatus and estimating seroprevalence**
- Lecture 5 (pptx and pdf versions)
  - Considerations for binarizing serological data, determining serostatus, and estimating seroprevalence
- Lab 5 (pdf and Rmd versions)
  - R code to binarize data to determine serostatus and estimate seroprevalence   

**Part 6: Inferring transmission dynamics from seroprevalence data**
- Lecture 6 (pptx and pdf versions)
  - A look into how to infer past transmission from age-specific seroprevalence patterns, timing and magnitude of past outbreaks in non-endemic settings, and force of infection in endemic transmission settings
- Lab 6 (pdf and Rmd versions)
  - R code to implement the serocatalytic model, infer timing and magnitude of past outbreaks, and estimate force of infection

**Source: General R functions and code used throughout this course**
- my_script.R
- utils.R

**Extra Lectures: Serodynamics and longitudinal analyses**
- Serodynamics review (pptx and pdf versions)
  - Review of methods using serological data to understand transmission dynamics and to determine disease exposure, including serocatalytic and more complex models
- Longitudinal serological analyses (pptx and pdf versions)
  - Review on using longitudinal serological data to evaluate antibody kinetics and the value of analyzing quantitative titers

**Projects: Guided worksheets to aid research question and analytic plan development** 

These worksheets can be applied to your own serosurveillance data and can be completed concurrently with the above lectures and labs.
