library(timeROC)  # Load the timeROC package for time-dependent ROC curve analysis
library(survival)  # Load the survival package for survival analysis

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
# Get the number of columns in matrix X, which is the number of predictors
p = ncol(X)

# Extract survival time and event status for progression-free interval (PFI) and overall survival (OS)
T1 = train_data$PFI.time  # PFI time
event1 = train_data$PFI  # PFI event status
T2 = train_data$OS.time  # OS time
event2 = train_data$OS  # OS event status

# Set the number of random initializations for optimization(in order to quickly out of the result Set to 5 actually should be Set as larger as possible)
Randomize_num = 5

# Define a function to calculate the negative log-likelihood for the Weibull model
Weibull_predict = function(y) {
    lambda1 = exp(y[1])  # Weibull scale parameter for PFI
    lambda2 = exp(y[2])  # Weibull scale parameter for OS
    rio1 = exp(y[3])  # Weibull shape parameter for PFI
    rio2 = exp(y[4])  # Weibull shape parameter for OS
    theta = exp(y[5])  # Shared parameter between PFI and OS
    beta1 = y[(5 + 1):(5 + p)]  # Coefficients for PFI predictors
    beta2 = y[(5 + p + 1):(5 + 2 * p)]  # Coefficients for OS predictors
    l = 0  # Initialize log-likelihood
    X1 = as.vector(X %*% beta1)  # Linear combination for PFI
    X2 = as.vector(X %*% beta2)  # Linear combination for OS
    r2 = as.vector(T2^(rio2 - 1) * lambda2 * rio2)  # Baseline hazard for OS
    R1 = as.vector(T1^rio1 * lambda1)  # Baseline cumulative hazard for PFI
    R2 = as.vector(T2^rio2 * lambda2)  # Baseline cumulative hazard for OS
    l = l + sum(event2 * (log(r2) + X2)) + sum(event2 * (theta + event1)) + 
        sum((event1 + theta) * log(theta + R1 * exp(X1))) - 
        sum((event1 + event2 + theta) * log(theta + R1 * exp(X1) + R2 * exp(X2)))
    -l  # Return the negative log-likelihood
}

# Initialize the number of optimization attempts and the maximum likelihood
p_num = 0
MPL = 0

# Loop to perform random initializations and find the best optimization result
repeat {
    p0 = runif(5 + 2 * p, -1, 1)  # Generate random initial parameters
    res_Random = try(optim(p0, Weibull_predict, method = "BFGS", hessian = TRUE), silent = TRUE)  # Perform optimization
    if (class(res_Random) != "try-error") {  # Check if optimization was successful
        if (res_Random$convergence == 0 & (min(eigen(res_Random$hessian)$values) > 0)) {
            res = res_Random  # Save the result if it converged and the Hessian is positive definite
            break
        }
        if (p_num >= Randomize_num) {
            res = res_Random  # Save the result if the number of attempts is reached
            break
        }
        MPL_Random = -res_Random$value  # Calculate the negative log-likelihood
        p_num = p_num + 1  # Increment the attempt counter
        if (MPL_Random > MPL) {
            res = res_Random  # Save the best result
            MPL = MPL_Random  # Update the maximum likelihood
        }
    } 
}

# Extract the optimized parameters from the result
lambda1 = exp(res$par[1])  # Weibull scale parameter for PFI
lambda2 = exp(res$par[2])  # Weibull scale parameter for OS
rio1 = exp(res$par[3])  # Weibull shape parameter for PFI
rio2 = exp(res$par[4])  # Weibull shape parameter for OS
theta = exp(res$par[5])  # Shared parameter
beta1 = res$par[(5 + 1):(5 + p)]  # Coefficients for PFI predictors
beta2 = res$par[(5 + p + 1):(5 + 2 * p)]  # Coefficients for OS predictors

# Initialize the matrix X for testing data
X <- NULL
for (i in seq_along(predict_variable)) {	
	X <- cbind(X, test_data[[predict_variable[i]]])
}
# Convert to matrix form
X = as.matrix(X)	
X1 = as.vector(X %*% beta1)  # Linear combination for PFI in the test data
X2 = as.vector(X %*% beta2)  # Linear combination for OS in the test data
R1 = as.vector(test_data$PFI.time^rio1 * lambda1)  # Baseline cumulative hazard for PFI in the test data

# Calculate the marker for ROC analysis and store it in column 20 of the test data
test_data[, 20] = -(theta * exp(-X2) + R1 * exp(X1 - X2))

# Perform time-dependent ROC analysis on the test data
ROC.Weibull <- timeROC(T = test_data$OS.time,
                       delta = test_data$OS,
                       marker = test_data$V20,
                       cause = 1,
                       weighting = "marginal",
                       times = quantile(test_data$OS.time, probs = seq(0, 1, 0.01)),
                       iid = TRUE) 

# Plot the AUC curve with confidence intervals
plotAUCcurve(ROC.Weibull, conf.int = TRUE, col = "aquamarine")
