library('mgcv')
library('gamCopula')
library('VineCopula')

##### A gamVine in 3 dimensions
set.seed(100)
## Simulation parameters
# Sample size
n <- 1e3
# Copula families
familyset <- c(1,301,401)
# Define a 3-dimensional R-vine tree structure matrix
d <- 3
Matrix <- c(1,3,2,0,3,2,0,0,2)
Matrix <- matrix(Matrix,d,d)
nnames <- paste("X", 1:d, sep = "")
# if we want directly to specify the parameters and the copula families we proceed this way:
# family <- c( 0, 1, 4,  0, 0, 1,0, 0, 0)
# family <- matrix(family, 3, 3)
# par <- c(0, 0.7, 3.9,0, 0, 0.9,0, 0, 0)
# par <- matrix(par, 3, 3)
# RVM <- RVineMatrix(Matrix = Matrix, family = family,par = par,names = c("V1", "V2", "V3"))
# contour(RVM)

# create a calibration surface:
# it is composed of three different functions-> linear, sinusoidal and exponential.
eta0 <- 1
calib.surf <- list(
  calib.quad <- function(t, Ti = 0, Tf = 1, b = 8) {
    Tm <- (Tf - Ti)/2
    a <- -(b/3) * (Tf^2 - 3 * Tf * Tm + 3 * Tm^2)
    return(a + b * (t - Tm)^2)},
  calib.sin <- function(t, Ti = 0, Tf = 1, b = 1, f = 1) {
    a <- b * (1 - 2 * Tf * pi/(f * Tf * pi +
                                 cos(2 * f * pi * (Tf - Ti))
                               - cos(2 * f * pi * Ti)))
    return((a + b)/2 + (b - a) * sin(2 * f * pi * (t - Ti))/2)},
  calib.exp <- function(t, Ti = 0, Tf = 1, b = 2, s = Tf/8) {
    Tm <- (Tf - Ti)/2
    a <- (b * s * sqrt(2 * pi)/Tf) * (pnorm(0, Tm, s) - pnorm(Tf, Tm, s))
    return(a + b * exp(-(t - Tm)^2/(2 * s^2)))})


## Create the model
# Define gam-vine model list
count <- 1
model <- vector(mode = "list", length = d*(d-1)/2)
sel <- seq(d,d^2-d, by = d)
# First tree
set.seed(100)
for (i in 1:(d-1)) {
  # Select a copula family
  family <- sample(familyset, 1)
  model[[count]]$family <- family
  # Use the canonical link and a randomly generated parameter
  if (is.element(family,c(1,2))) {
    model[[count]]$par <- tanh(rnorm(1)/2)
    if (family == 2) {model[[count]]$par2 <- 2+exp(rnorm(1))
    }
  } else {
    if (is.element(family,c(401:404))) {
      rr <- rnorm(1)
      model[[count]]$par <- sign(rr)*(1+abs(rr))
    } else {
      model[[count]]$par <- rnorm(1)
    }
    model[[count]]$par2 <- 0
  }
  count <- count + 1
}

# A dummy dataset
set.seed(100)
data <- data.frame(u1 = runif(1e2), u2 = runif(1e2), matrix(runif(1e2*d),1e2,d))
# Trees 2 to (d-1)
for(j in 2:(d-1)){
  for(i in 1:(d-j)){
    # Select a copula family
    family <- sample(familyset, 1)
    # Select the conditiong set and create a model formula
    cond <- nnames[sort(Matrix[(d-j+2):d,i])]
    tmpform <- paste("~",paste(paste("s(", cond, ", k=10, bs='cr')",
                                     sep = ""), collapse=" + "))
    l <- length(cond)
    temp <- sample(3, l, replace = TRUE)
    # Spline approximation of the true function
    m <- 1e2
    x <- matrix(seq(0,1,length.out=m), nrow = m, ncol = 1)
    if(l != 1){
      tmp.fct <- paste("function(x){eta0+",
                       paste(sapply(1:l, function(x)
                         paste("calib.surf[[",temp[x],"]](x[",x,"])",
                               sep="")), collapse="+"),"}",sep="")
      tmp.fct <- eval(parse(text = tmp.fct))
      x <- eval(parse(text = paste0("expand.grid(",
                                    paste0(rep("x",l), collapse = ","),")",
                                    collapse = "")))
      y <- apply(x,1,tmp.fct)
    }else{
      tmp.fct <- function(x) eta0+calib.surf[[temp]](x)
      colnames(x) <- cond
      y <- tmp.fct(x)
    }
    # Estimate the gam model
    form <- as.formula(paste0("y", tmpform))
    dd <- data.frame(y, x)
    names(dd) <- c("y", cond)
    b <- gam(form, data = dd)
    plot(b) #to check which function we have selected
    #plot(x[,1],(y-fitted(b))/y)
    # Create a dummy gamBiCop object
    tmp <- gamBiCopFit(data = data, formula = form, family = 1, n.iters = 1)$res
    # Update the copula family and the model coefficients
    attr(tmp, "model")$coefficients <- coefficients(b)
    attr(tmp, "model")$smooth <- b$smooth
    attr(tmp, "family") <- family
    if (family == 2) {
      attr(tmp, "par2") <- 2+exp(rnorm(1))
    }
    model[[count]] <- tmp
    count <- count+1
  }
}

# Create the gamVineCopula object, notice that we could add covariates here.
GVC <- gamVine(Matrix=Matrix,model = model,names=nnames)
print(GVC)

# we are going to simulate data from this model, and lately we will try
# to fit it with copulas with simplied and non-simplified assumptions, and we will 
# compare its performances.
sim <- gamVineSimulate(n, GVC)

# fitting with all the information:
fit1GVC <- gamVineSeqFit(sim, GVC, verbose = TRUE)
plot(fit1GVC)
nobs(fit1GVC)
# fitting only with the matrix information:
fit2GVC <- gamVineCopSelect(sim, Matrix, verbose = TRUE)
# fitting non simplified copulas without any information:
fit3GVC <- gamVineStructureSelect(sim, lin.covs = NULL, smooth.covs = NULL,
                                  simplified = FALSE,verbose= TRUE, type = 0,
                                  familyset <- c(1:2,301:304,401:404))

# fitting simplified copulas with the matrix information:
RVM1 <- RVineCopSelect(sim,familyset = NA,Matrix,selectioncrit = "AIC",indeptest = FALSE,
  level = 0.05,trunclevel = NA,weights = NA,rotations = TRUE,se = FALSE,presel = TRUE,
  method = "mle",cores = 1)
RVM1$AIC

# fitting simplified copulas without any information:
RVM2 <- RVineStructureSelect(sim,familyset = NA,
                           type = 0, selectioncrit = "AIC",indeptest = FALSE,
                           level = 0.05,trunclevel = NA,progress = FALSE,weights = NA,
                           treecrit = "tau",rotations = TRUE,se = FALSE,presel = TRUE,
                           method = "mle",cores = 1)
RVM2$AIC



# Now we will study the logLikelihood our non-simplified model evaluated at our sim. data
# for the 1st copula: Clayton we have:
#term1<-log(prod(BiCopPDF(X1, X2, family=3, par=2, par2=0, obj=NULL, check.pars=TRUE)))
# for the 2nd copula: Gumbell we have:
#term2<-log(prod(BiCopPDF(X2, X3, family=4, par=1.5, par2=0, obj=NULL, check.pars=TRUE)))
# we have to compute the conditional values:
#X1_2<-BiCopHfunc(X1, X2, family=3, 2, par2 = 0, obj = NULL, check.pars = TRUE)$hfunc2
#X3_2<-BiCopHfunc(X3, X2, family=3, 2, par2 = 0, obj = NULL, check.pars = TRUE)$hfunc2
# compute the rho parameter of the gaussian copula for each x2:
#rho<-sin(pi/2*(exp(1-pi+pi*sin(X2*pi))-1)/(exp(1-pi+pi*sin(X2*pi))+1))
# for the 3rd copula: Gaussian we have:
#term3<-prod(BiCopPDF(X1_2 , X3_2, family=1, par=rho, par2=0, obj=NULL, check.pars=TRUE))