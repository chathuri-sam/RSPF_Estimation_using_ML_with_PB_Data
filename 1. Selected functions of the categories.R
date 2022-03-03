##Part 1
#Observed Logistic parametric functions fitted with logistic functions
#Each category has 3 functions, linear, quadratic and cubic.


#####################################################################
#####################################################################
#Category 1
#Both RSPF and local certainty conditions satisfied.

#Linear 
logistic_3 <- function(x) plogis(0.606 - 3.64*x)
#Quadratic
Qlogistic_3 <- function(x) plogis(0.606 - 3.64*x + 1.26*x^2)
#Cubic
Clogistic_4 <- function(x) plogis(-0.284 -2.298*x- 0.232*x^2 + 7.269*x^3)


#####################################################################
#####################################################################

#Category 2
#Only local certainty condition satisfied
#Linear
semilogistic <- function(x) 8.3*plogis(-4+2*x)
#Quadratic
gaussian <- function(x) 0.99*exp(-(4*x-2)^2)
#Cubic
exponent <- function(x) exp(-0.5 + 0.5*x - 0.6*x^2 - 0.8*x^3)

#####################################################################
#####################################################################
#Category 3
#Only RSPF Condition satisfied.

#Linear
logistic_1 <- function(x) plogis(0.5855 + 1.064*x)
#Quadratic  
Qlogistic_1 <- function(x) plogis(0.5855 + 1.064*x- 0.218*x^2)
#Cubic 
Clogistic_1 <- function(x) plogis(0.5855 + 1.064*x- 0.218*x^2- 1.81*x^3)

#####################################################################
#####################################################################
#####################################################################
#####################################################################
##Part 2
#Observed cloglog parametric functions fitted with logistic functions 
#(Mis-specfied models)
#Each category has 3 functions, linear, quadratic and cubic.

#####################################################################
#####################################################################
#Category 1
#Both RSPF and local certainty conditions satisfied.

#Linear 
cloglog_1 <- function(x) binomial("cloglog")$linkinv(0.5855 + 1.064*x)
#Quadratic
Qcloglog_1 <- function(x) binomial("cloglog")$linkinv(0.5855 + 1.064*x- 0.218*x^2)
#Cubic
Ccloglog_2 <- function(x) binomial("cloglog")$linkinv(0.37 + 1.56*x - 1.5*x^2 + 3.267*x^3)


#####################################################################
#####################################################################
#Category 2
#Only RSPF conditions satisfied.

#Linear 
cloglog_5 <- function(x) binomial("cloglog")$linkinv(-0.886 - 1.16*x)
#Quadratic
Qcloglog_2 <- function(x) binomial("cloglog")$linkinv(0.37 + 1.56*x - 1.5*x^2)
#Cubic
Ccloglog <- function(x)binomial("cloglog")$linkinv(0.1 - 0.064*x - 0.85*x^2 - 0.81*x^3)

#####################################################################
#####################################################################
#Hastie and Fithian (2013) functions used
#Full logistic
hastie <- function(x) plogis(-1+x)
#Scaled logistic
hastie_half <- function(x) 0.5*plogis(-1+x)
