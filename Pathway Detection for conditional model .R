# Read the data file "GBMLGG.txt", which is tab-separated, with the first row as the header
data=read.table("GBMLGG.txt", sep = "\t", header = TRUE)

set.seed(123)  # Set the random seed for reproducibility

# Set the number of random initializations for optimization(in order to quickly out of the result Set to 20 actually should be Set as larger as possible)
Randomize_num = 20

# Exponential model
Exponential_LRT = function (T1, event1, T2, event2, X){
	X = as.matrix(X) # Convert to matrix form
	p = ncol(X) # number of variables
	# likelihood for H1
    likelihood_H1 = function(y) {
        lambda1 = exp(y[1])
        lambda2 = exp(y[2])
        theta = exp(y[3])
        beta1 = y[(3 + 1):(3 + p)]
        beta2 = y[(3 + p + 1):(3 + 2 * p)]
        l = 0
        X1 = as.vector(X %*% beta1)
        X2 = as.vector(X %*% beta2)
        r2 = as.vector(lambda2)
        R1 = as.vector(T1 * lambda1)
        R2 = as.vector(T2 * lambda2)
        l = l + sum(event2 * (log(r2) + X2)) + sum(event2 * (theta + event1)) + 
        sum((event1 + theta) * log(theta + R1 * exp(X1))) - 
        sum((event1 + event2 + theta) * log(theta + R1 * exp(X1) + R2 * exp(X2)))
        -l
	}
	# Initialize the number of optimization attempts and the maximum likelihood
	p_num = 0
	MPL = 0
	repeat {
    		p0 = runif(3 + 2 * p, -1, 1)  # Generate random initial parameters
    		res_Random = try(optim(p0, likelihood_H1, method = "BFGS", hessian = TRUE), silent = TRUE)  # Perform optimization
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
    para=res$par 
    MPL_H1=-res$value # likelihood value for H1
    # likelihood for H0
    likelihood_H0 = function(y) {
        lambda1 = exp(y[1])
        lambda2 = exp(y[2])
        theta = exp(y[3])
        beta1 = 0
        beta2 = 0
        l = 0
        X1 = as.vector(X %*% beta1)
        X2 = as.vector(X %*% beta2)
        r2 = as.vector(lambda2)
        R1 = as.vector(T1 * lambda1)
        R2 = as.vector(T2 * lambda2)
        l = l + sum(event2 * (log(r2) + X2)) + sum(event2 * (theta + event1)) + 
        sum((event1 + theta) * log(theta + R1 * exp(X1))) - 
        sum((event1 + event2 + theta) * log(theta + R1 * exp(X1) + R2 * exp(X2)))
        -l
	}
	# Initialize the number of optimization attempts and the maximum likelihood
	p_num = 0
	MPL = 0
	repeat {
    		p0 = runif(3, -1, 1)  # Generate random initial parameters
    		res_Random = try(optim(p0, likelihood_H0, method = "BFGS", hessian = TRUE), silent = TRUE)  # Perform optimization
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
    MPL_H0=-res$value # likelihood value for H0
	LR=2*(MPL_H1-MPL_H0) # LRT value
	c(LR,para)
}

para1=NULL

tes=Exponential_LRT(T1=data$PFI.time,event1=data$PFI,T2=data$OS.time,event2=data$OS,X=data$WNT)
	P1=1-pchisq(tes[1],df=2)
	beta1=tes[5]
	beta2=tes[6]
	para1=tes[2:6]
	if(beta2>=0 & beta1<=beta2){
		mono1=1
	}
	if(beta2>=0 & beta1>beta2){
		mono1=2
	}
	if(beta2<0 & beta1>=beta2){
		mono1=3
	}
	if(beta2<0 & beta1<beta2){
		mono1=4
	}

# P-value
print(P1)

# parameter
print(para1)

# monotonicity
print(mono1)


set.seed(123)  # Set the random seed for reproducibility

# Set the number of random initializations for optimization(in order to quickly out of the result Set to 20 actually should be Set as larger as possible)
Randomize_num = 20

# Weibull model
Weibull_LRT = function (T1, event1, T2, event2, X){
	X = as.matrix(X) # Convert to matrix form
	p = ncol(X) # number of variables
	# likelihood for H1
    likelihood_H1 = function(y) {
        lambda1 = exp(y[1])
        lambda2 = exp(y[2])
        rio1 = exp(y[3])
        rio2 = exp(y[4])
        theta = exp(y[5])
        beta1 = y[(5 + 1):(5 + p)]
        beta2 = y[(5 + p + 1):(5 + 2 * p)]
        l = 0
        X1 = as.vector(X %*% beta1)
        X2 = as.vector(X %*% beta2)
        r2 = as.vector(T2^(rio2-1) * lambda2 * rio2)
        R1 = as.vector(T1^(rio1) * lambda1)
        R2 = as.vector(T2^(rio2) * lambda2)
        l = l + sum(event2 * (log(r2) + X2)) + sum(event2 * (theta + event1)) + 
        sum((event1 + theta) * log(theta + R1 * exp(X1))) - 
        sum((event1 + event2 + theta) * log(theta + R1 * exp(X1) + R2 * exp(X2)))
        -l
	}
	# Initialize the number of optimization attempts and the maximum likelihood
	p_num = 0
	MPL = 0
	repeat {
    		p0 = runif(5 + 2 * p, -1, 1)  # Generate random initial parameters
    		res_Random = try(optim(p0, likelihood_H1, method = "BFGS", hessian = TRUE), silent = TRUE)  # Perform optimization
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
    para=res$par
    MPL_H1=-res$value # likelihood value for H1
    # likelihood for H0
    likelihood_H0 = function(y) {
        lambda1 = exp(y[1])
        lambda2 = exp(y[2])
        rio1 = exp(y[3])
        rio2 = exp(y[4])
        theta = exp(y[5])
        beta1 = 0
        beta2 = 0
        l = 0
        X1 = as.vector(X %*% beta1)
        X2 = as.vector(X %*% beta2)
        r2 = as.vector(lambda2)
        R1 = as.vector(T1 * lambda1)
        R2 = as.vector(T2 * lambda2)
        l = l + sum(event2 * (log(r2) + X2)) + sum(event2 * (theta + event1)) + 
        sum((event1 + theta) * log(theta + R1 * exp(X1))) - 
        sum((event1 + event2 + theta) * log(theta + R1 * exp(X1) + R2 * exp(X2)))
        -l
	}
	# Initialize the number of optimization attempts and the maximum likelihood
	p_num = 0
	MPL = 0
	repeat {
    		p0 = runif(5, -1, 1)  # Generate random initial parameters
    		res_Random = try(optim(p0, likelihood_H0, method = "BFGS", hessian = TRUE), silent = TRUE)  # Perform optimization
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
    MPL_H0=-res$value # likelihood value for H0
	LR=2*(MPL_H1-MPL_H0) # LRT value
	c(LR,para)
}

para2=NULL
tes= Weibull_LRT(T1=data$PFI.time,event1=data$PFI,T2=data$OS.time,event2=data$OS,X=data$WNT)
	P2=1-pchisq(tes[1],df=2)
	beta1=tes[7]
	beta2=tes[8]
	para2=tes[2:8]
	if(beta2>=0 & beta1<=beta2){
		mono2=1
	}
	if(beta2>=0 & beta1>beta2){
		mono2=2
	}
	if(beta2<0 & beta1>=beta2){
		mono2=3
	}
	if(beta2<0 & beta1<beta2){
		mono2=4
	}

# P-value
print(P2)

# parameter
print(para2)

# monotonicity
print(mono2)







