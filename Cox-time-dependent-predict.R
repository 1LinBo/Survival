library(survival)  # Load the survival package for survival analysis
library(pROC)  # Load the pROC package for ROC curve analysis

set.seed(111)  # Set the random seed for reproducibility

# Read the data file "LUNG.txt", which is tab-separated, with the first row as the header
df = read.table("LUNG.txt", sep = "\t", header = TRUE)

# Extract the variable names from the second to the fifteenth columns
predict_variable = colnames(df)[2:15]
predict_variable  # Print the variable names

# Get the number of rows in the data frame
n <- nrow(df)

# Randomly select 65% of the data as the training set
train_indices <- sample(1:n, size = 0.65 * n)

# Split the data into training and testing sets based on the indices
train_data <- df[train_indices, ]  # Training data
test_data <- df[-train_indices, ]  # Testing data

# Initialize the matrix X for training data
X <- NULL

# Combine the predictor variables from the training data into matrix X
for (i in seq_along(predict_variable)) {	
	X <- cbind(X, train_data[[predict_variable[i]]])
}
# Convert to matrix form
X = as.matrix(X)

## time-dependent coefficients analysis

# Create a data frame for the Cox model with the survival time, event status, predictor variables,
# and the logarithm of PFI time
train <- data.frame(
  OS.time = train_data$OS.time,  # Overall survival time
  OS = train_data$OS,  # Overall survival event status
  X,  # Predictor variables
  log_PFI_time = log(train_data$PFI.time)  # Logarithm of PFI time
)

# Fit a Cox proportional hazards model using the training data
fit = coxph(Surv(OS.time, OS) ~ ., data = train)

# Test the proportional hazards assumption for the Cox model
zp = cox.zph(fit)
zp  # Print the test results

# Plot the Schoenfeld residuals for the log_PFI_time variable (15th covariate)
plot(zp[15])

# Split the training data into time intervals for time-dependent analysis
# The `cut` parameter specifies the time points at which the data is split
train_group = survSplit(Surv(OS.time, OS) ~ ., data = train, episode = "tgroup", cut = c(590, 1200))

# Fit a Cox model with time-dependent coefficients, using the split training data
fit2 = coxph(Surv(tstart, OS.time, OS) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + 
             log_PFI_time:strata(tgroup), data = train_group)
fit2  # Print the model summary


# Test the proportional hazards assumption for the time-dependent Cox model
cox.zph(fit2)

# Initialize the matrix X for testing data
X <- NULL
for (i in seq_along(predict_variable)) {	
	X <- cbind(X, test_data[[predict_variable[i]]])
}
# Convert to matrix form
X = as.matrix(X)

# Define the time interval for ROC analysis
time_start <- 0  # Start of the interval
time_end <- 365  # End of the interval (1 year)
# Define the event status for OS at the end of the time interval (1 year)
event_OS_365 <- with(test_data, ifelse(OS.time > time_start & OS.time <= time_end & OS == 1, 1, 0))

# Create a data frame for the test data with the survival time, event status,
# predictor variables, logarithm of PFI time, and the event status for 1 year OS
test <- data.frame(
  OS.time = test_data$OS.time,  # Overall survival time
  OS = test_data$OS,  # Overall survival event status
  X,  # Predictor variables
  log_PFI_time = log(test_data$PFI.time),  # Logarithm of PFI time
  event_OS_365 = event_OS_365  # Event status for 1 year OS
)

# Split the test data into time intervals for time-dependent analysis
test_group = survSplit(Surv(OS.time, OS) ~ ., data = test, episode = "tgroup", cut = c(590, 1200))

# Predict the risk scores for the test data using the time-dependent Cox model
test_predict <- data.frame(
  test_group,
  risk_predict = predict(fit2, newdata = test_group, type = "risk")
)

# Create a final test data frame with risk scores for the first time interval
test_final <- data.frame(
  test,
  risk_predict = test_predict[test_predict$tgroup == 1, ]$risk_predict
)

# Compute the ROC curve for the predicted risk scores against the 1-year OS event status
roc_curve <- roc(test_final$event_OS_365, test_final$risk_predict, quiet = TRUE)

# Plot the ROC curve
plot(roc_curve, col = "blue", main = paste("ROC Curve at Time =", 365, "days"))
abline(a = 0, b = 1, col = "red", lty = 2)  # Add a diagonal line for reference
auc_value <- auc(roc_curve)  # Calculate the area under the ROC curve (AUC)
legend("bottomleft", legend = paste("AUC =", round(auc_value, 3)), bty = "n", col = "blue")


## time-dependent covariates analysis

train2 <- data.frame(
  OS.time = train_data$OS.time,  # Overall survival time
  OS = train_data$OS,  # Overall survival event status
  X,  # Predictor variables
  PFI.time = train_data$PFI.time  # PFI time
)

library(dplyr)

# Add the ID column using the mutate function
train2 <- train2 %>%
  mutate(id = row_number())

transformed_data_list <- list()

for (i in 1:nrow(train2)) {
  patient <- train2[i, ]
  if (!is.na(patient$PFI.time)) {
    # For patients who progress, two time intervals are created
    interval1 <- patient
    interval1$tstart <- 0
    interval1$tstop <- patient$PFI.time
    interval1$PFI <- 0
    interval1$status <- 0
    interval2 <- patient
    interval2$tstart <- patient$PFI.time
    interval2$tstop <- patient$OS.time
    interval2$PFI <- 1
    interval2$status <- patient$OS
    transformed_data_list <- append(transformed_data_list, list(interval1, interval2))
  } 
}

transformed_data <- do.call(rbind, transformed_data_list)

fit3=coxph(Surv(tstart, tstop, status) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + PFI,
      data =transformed_data, cluster = id)

# Initialize the matrix X for testing data
X <- NULL
for (i in seq_along(predict_variable)) {	
  X <- cbind(X, test_data[[predict_variable[i]]])
}
# Convert to matrix form
X = as.matrix(X)

# Define the time interval for ROC analysis
time_start <- 0  # Start of the interval
time_end <- 365  # End of the interval (1 year)
# Define the event status for OS at the end of the time interval (1 year)
event_OS_365 <- with(test_data, ifelse(OS.time > time_start & OS.time <= time_end & OS == 1, 1, 0))

# Create a data frame for the test data with the survival time, event status,
# predictor variables, logarithm of PFI time, and the event status for 1 year OS
test2 <- data.frame(
  OS.time = test_data$OS.time,  # Overall survival time
  OS = test_data$OS,  # Overall survival event status
  X,  # Predictor variables
  PFI.time = test_data$PFI.time,  # PFI time
  event_OS_365 = event_OS_365  # Event status for 1 year OS
)

# Add the ID column using the mutate function
test2 <- test2 %>%
  mutate(id = row_number())

transformed_data_list <- list()

for (i in 1:nrow(test2)) {
  patient <- test2[i, ]
  if (!is.na(patient$PFI.time)) {
    # For patients who progress, two time intervals are created
    interval1 <- patient
    interval1$tstart <- 0
    interval1$tstop <- patient$PFI.time
    interval1$PFI <- 0
    interval1$status <- 0
    interval2 <- patient
    interval2$tstart <- patient$PFI.time
    interval2$tstop <- patient$OS.time
    interval2$PFI <- 1
    interval2$status <- patient$OS
    transformed_data_list <- append(transformed_data_list, list(interval1, interval2))
  } 
}

transformed_data <- do.call(rbind, transformed_data_list)

# Predict the risk scores for the test data using the time-dependent Cox model
test_predict <- data.frame(
  transformed_data,
  risk_predict = predict(fit3, newdata = transformed_data, type = "risk")
)

# Create a final test data frame with risk scores for the first time interval
test_final <- data.frame(
  test2,
  risk_predict = test_predict[(test_predict$tstop>=365 & test_predict$PFI == 0) | (test_predict$tstart<365 & test_predict$PFI == 1), ]$risk_predict
)

# Compute the ROC curve for the predicted risk scores against the 1-year OS event status
roc_curve <- roc(test_final$event_OS_365, test_final$risk_predict, quiet = TRUE)

# Plot the ROC curve
plot(roc_curve, col = "blue", main = paste("ROC Curve at Time =", 365, "days"))
abline(a = 0, b = 1, col = "red", lty = 2)  # Add a diagonal line for reference
auc_value <- auc(roc_curve)  # Calculate the area under the ROC curve (AUC)
legend("bottomleft", legend = paste("AUC =", round(auc_value, 3)), bty = "n", col = "blue")





