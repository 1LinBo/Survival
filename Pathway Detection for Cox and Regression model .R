library(survival)  # Load the survival package for survival analysis

# Read the data file "LUNG.txt", which is tab-separated, with the first row as the header
final = read.table("LUNG.txt", sep = "\t", header = TRUE)

# Cox Model 1
# Fit a Cox proportional hazards model using overall survival time (OS.time) and event status (OS)
# The model includes 'P53' as a predictor
tes = coxph(Surv(OS.time, OS) ~ P53, data = final)

# Print the summary of the Cox model, focusing on the coefficients
summary(tes)$coefficients

# Cox Model 2
# Fit a Cox proportional hazards model including 'P53' and the logarithm of progression-free interval (PFI) time
tes = coxph(Surv(OS.time, OS) ~ P53 + log(final$PFI.time), data = final)

# Print the summary of the Cox model and extract the first coefficient (P53)
summary(tes)$coefficients[1, ]

# Cox Time-Dependent Model
# Create a data frame for the Cox model with the survival time, event status, 'P53', and logarithm of PFI time
data1 <- data.frame(
  OS.time = final$OS.time,  # Overall survival time
  OS = final$OS,  # Overall survival event status
  X = final$P53,  # 'P53' predictor variable
  log_PFI_time = log(final$PFI.time)  # Logarithm of PFI time
)

# Fit a Cox proportional hazards model using the data frame
fit = coxph(Surv(OS.time, OS) ~ ., data = data1)

# Test the proportional hazards assumption for the Cox model
zp = cox.zph(fit)
zp  # Print the test results

# Plot the Schoenfeld residuals for the log_PFI_time variable (2nd covariate)
plot(zp[2])

# Split the data into time intervals for time-dependent analysis
# The `cut` parameter specifies the time points at which the data is split
data1_group = survSplit(Surv(OS.time, OS) ~ ., data = data1, episode = "tgroup", cut = c(580, 1200))

# Fit a Cox model with time-dependent coefficients, using the split data
fit2 = coxph(Surv(tstart, OS.time, OS) ~ X + log_PFI_time:strata(tgroup), data = data1_group)
fit2  # Print the model summary

# Test the proportional hazards assumption for the time-dependent Cox model
cox.zph(fit2)

# Print the summary of the time-dependent Cox model and extract the first coefficient
summary(fit2)$coefficients[1, ]

# Weibull Regression Model 1
# Fit a Weibull regression model using 'P53' as a predictor
tes = survreg(Surv(OS.time, OS) ~ P53, dist = "weibull", data = final)

# Print the summary of the Weibull model, focusing on the coefficient of 'P53'
summary(tes)$table[2, ]

# Weibull Regression Model 2
# Fit a Weibull regression model including 'P53' and the logarithm of PFI time
tes = survreg(Surv(OS.time, OS) ~ P53 + log(final$PFI.time), dist = "weibull", data = final)

# Print the summary of the Weibull model, focusing on the coefficient of 'P53'
summary(tes)$table[2, ]
