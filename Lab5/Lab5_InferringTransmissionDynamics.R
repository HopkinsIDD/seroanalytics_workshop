###----- Lab 5: Inferring transmission dynamics -----###
#  Megan O'Driscoll 
# email: megan.odriscoll@hug.ch 


#--- 1. Set up our environment ---#
rm(list=ls()) # clear everything in the environment

# load packages
library(ggplot2)
library(epitools)


#--- 2. Read in the data ---#
setwd("C:/Users/megan/Documents/GitHub/seroanalytics_workshop/Lab5") # need to change working dir
df <- read.csv("simulated_data.csv")

# lets take a look at our dataset
# in this lab we will focus on the dengue (DENV) and chikungunya (CHIKV) data
head(df)


#--- 3. Calculate age-specific seroprevalence ---#

# we want to calculate age-specific seroprevalence
# and we have already determined individual-level serostatuses
# let's start by assigning everyone to an age group

hist(df$age) # check the age distribution of the study population
age_groups <- c("0-4","5-9","10-14","15-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")
df$age_group <- cut(df$age, c(-Inf,4,9,14,19,29,39,49,59,69,79,Inf),
                    labels=age_groups)


# now that we have age groups we can calculate seroprevalence by age
# here we will create a new dataframe for seroprevalence

seroprev <- list()
for(p in c("CHIKV","DENV")){
  
  # create a dataframe and store it in the list
  seroprev[[p]] <- data.frame(age_group=factor(age_groups, levels=age_groups), n=NA, npos=NA)
  
  for(a in 1:length(age_groups)){ # loop through each age group
    seroprev[[p]]$n[a] <- nrow(df[df$age_group==age_groups[a], ]) # how many people in this age group
    seroprev[[p]]$npos[a] <- nrow(df[df$age_group==age_groups[a] & df[,p]==1, ]) # how many positive in this age group
  }
  
  # we can now calculate seroprevalence (npos/n)
  seroprev[[p]]$prev <- seroprev[[p]]$npos/seroprev[[p]]$n
  
  # and binomial confidence intervals
  seroprev[[p]][,c("ciL","ciU")] <- binom.exact(seroprev[[p]]$npos, seroprev[[p]]$n)[,c("lower","upper")]
  
}


#--- 4. Investigate past transmission dynamics ---#

# lets look at age-specific CHIKV trends

ggplot(seroprev$CHIKV, aes(age_group, prev))+ geom_point()+ ggtitle("Chikungunya")+ 
  geom_linerange(aes(ymin=ciL, ymax=ciU))+ ylab("Seroprevalence")+ xlab("Age group")


#*** Question: What do we think has happened with CHIKV transmission here? 
 

# lets now check the trends for DENV seroprevalence

ggplot(seroprev$DENV, aes(age_group, prev))+ geom_point()+ ggtitle("Dengue")+ ylim(0,1)+
  geom_linerange(aes(ymin=ciL, ymax=ciU))+ ylab("Seroprevalence")+ xlab("Age group")


#*** Question: what kind of transmission pattern do we see for DENV?


#--- 5. Estimating FOI (force of infection) ---#

# we will fit a simple complementary log-log (cloglog) model 
# to do this we will use the mid-points of each age group 

age_mid <- c(2, 7, 12, 17, 24.5, 34.5, 44.5, 54.5, 64.5, 74.5, 90) # age group mid-points
seroprev$DENV$age_mid <- age_mid

mod <- glm(cbind(npos, n-npos) ~ 1, offset=log(age_mid), data=seroprev$DENV, family = binomial(link = "cloglog"))

summary(mod)
log_FOI <- mod$coefficients[1]
FOI <- exp(log_FOI)

# our model estimates a 6.6% force of infection
# lets calculate the confidence intervals around this estimate

se_log_FOI <- summary(mod)$coefficients[1,2]
ci_log_FOI <- log_FOI + c(-1.96, 1.96) * se_log_FOI
ci_FOI <- exp(ci_log_FOI)


# lets now see how well this estimate fits our data

# calculate predicted seroprevalence by age group
pred_seroprev <- data.frame(age_group=seroprev$DENV$age_group,
                            pred=1-exp(-FOI*age_mid), 
                            ciL=1-exp(-ci_FOI[1]*age_mid),
                            ciU=1-exp(-ci_FOI[2]*age_mid))


# plot observed vs model estimated seroprevalence
ggplot(seroprev$DENV, aes(age_group, prev))+ geom_point()+ ggtitle("Dengue")+ ylim(0,1)+
  geom_linerange(aes(ymin=ciL, ymax=ciU))+ ylab("Seroprevalence")+ xlab("Age group")+
  geom_line(data=pred_seroprev, aes(age_group, pred, group=1), col="blue")+
  geom_ribbon(data=pred_seroprev, aes(x=age_group, y=pred, ymin=ciL, ymax=ciU, group=1), fill="blue", alpha=0.2)


# we can see that the model estimate of FOI matches our observed data pretty well!


#--- 6. Further exploration of FOI (extra exercises if time allows)

# lets see what other values of FOI can look like (feel free to add your own values here)
explore_FOIs <- c(0.01, 0.05, 0.1, 0.2, 0.3)

# we will loop through each of these values and simulate expected age-specific seroprevalence values
explore_data <- list() # create a list to store the simulated data
for(i in 1:length(explore_FOIs)){
  
  edf <- data.frame(age=seq(0,90), 
                    FOI=explore_FOIs[i],
                    prev=1-exp(-explore_FOIs[i]*seq(0,90)),
                    susc=exp(-explore_FOIs[i]*seq(0,90)))
  explore_data[[i]] <- edf
  
}

# here we unlist and bind the data together for plotting
explore_data <- do.call("rbind", explore_data)

# plotting expected age-specific prevalence for different FOIs
ggplot(explore_data, aes(age, prev, col=factor(FOI)))+ geom_line()+
  ylab("Seroprevalence")+ xlab("Age (years)")+ labs(col="FOI")

#*** Question: What proportion of 10 year olds would we expect to have been infected in each scenario?

ggplot(explore_data, aes(age, prev, col=factor(FOI)))+ geom_line()+
  ylab("Seroprevalence")+ xlab("Age (years)")+ labs(col="FOI")+
  geom_vline(aes(xintercept=10), linetype="dashed")

# now lets also see what the proportion susceptible would look like in each scenario
ggplot(explore_data, aes(age, susc, col=factor(FOI)))+ geom_line()+
  ylab("Proportion Susceptible")+ xlab("Age (years)")+ labs(col="FOI")
