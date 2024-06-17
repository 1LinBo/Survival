library(survival)  # Load the survival package for survival analysis
library(timeROC)  # Load the timeROC package for time-dependent ROC curve analysis

set.seed(123)  # Set the random seed for reproducibility

# Read the data file "GBMLGG.txt", which is tab-separated, with the first row as the header
df = read.table("GBMLGG.txt", sep = "\t", header = TRUE)

# Extract the variable names from the second to the fifteenth columns
predict_variable = colnames(df)[2:15]
predict_variable  # Print the variable names

# Get the number of rows in the data frame
n <- nrow(df)

# Randomly select 60% of the data as the training set
train_indices <- sample(1:n, size = 0.6 * n)

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

# Fit a Weibull regression model using overall survival time (OS.time) and event status (OS)
# The `survreg` function is used with the Weibull distribution specified by `dist = "weibull"`
fit = survreg(Surv(train_data$OS.time, train_data$OS) ~ X, dist = "weibull")

# Extract the coefficients from the Weibull regression model, excluding the intercept
coefficients <- coef(fit)[-1]

# Initialize the matrix X for testing data
X <- NULL
for (i in seq_along(predict_variable)) {	
	X <- cbind(X, test_data[[predict_variable[i]]])
}
# Convert to matrix form
X = as.matrix(X)

# Calculate the risk score for the test data using the Weibull regression model coefficients and store the result in column 20 of `test_data`
test_data[, 20] = -X %*% coefficients

# Perform time-dependent ROC analysis on the test data
ROC.regression <- timeROC(T = test_data$OS.time,
                          delta = test_data$OS,
                          marker = test_data$V20,
                          cause = 1,
                          weighting = "marginal",
                          times = quantile(test_data$OS.time, probs = seq(0, 1, 0.01)),
                          iid = TRUE)

# Plot the AUC curve with confidence intervals
plotAUCcurve(ROC.regression, conf.int = TRUE, col = "aquamarine")
