library('mgcv')
library('gamCopula')
library('VineCopula')
library('ggplot2')
library('latex2exp')
library('skimr')
library('kdecopula')
library('pacotest')
library('plot3D')
library('plotly')

########################################
# dimension of the copula that we want to construct.
d <- 3
# matrix to specify the vine copula model.
Matrix <- c(1,3,2,0,3,2,0,0,2)
Matrix <- matrix(Matrix,d,d)
nnames <- paste("X", 1:d, sep = "")

# Possible functions that we will use to define the relation of the gam:
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

# we will proceed first with the sinusoidal function.
model <- vector(mode = "list", length = d*(d-1)/2)
sel <- seq(d,d^2-d, by = d)
# specifying the first copula -> Clayton.
model[[1]]$family <- 301
model[[1]]$par <- 2
model[[1]]$par2 <- 0
# specifying the second copula -> Gumbell.
model[[2]]$family <- 401
model[[2]]$par <- 1.5
model[[2]]$par2 <- 0
# specifying the third copula -> Gaussian.
# a dummy data set that we will use to fit the gam model.
set.seed(50) # reproducibility
data <- data.frame(u1 = runif(1e2), u2 = runif(1e2), matrix(runif(1e2*d),1e2,d))

cond <- nnames[2]
tmpform <- paste("~",paste(paste("s(", cond, ", k=10, bs='cr')",sep = ""), collapse=" + "))
m <- 1e2
x <- matrix(seq(0,1,length.out=m), nrow = m, ncol = 1)
tmp.fct <- function(x) eta0+calib.surf[[2]](x) # we first proceed with the sinusoidal.
colnames(x) <- cond
y <- tmp.fct(x)
form <- as.formula(paste0("y", tmpform))
dd <- data.frame(y, x)
names(dd) <- c("y", cond)
b1 <- gam(form, data = dd)
plot(b1,main='fitted Sinusoidal function') 
tmp <- gamBiCopFit(data = data, formula = form, family = 1, n.iters = 1)$res
attr(tmp, "model")$coefficients <- coefficients(b1)
attr(tmp, "model")$smooth <- b1$smooth
attr(tmp, "family") <- 1
model[[3]] <- tmp

# we define our Gam Vine Copula:
GVC <- gamVine(Matrix=Matrix,model = model,names=nnames)
summary(GVC)
##########################################
# once the model is defined, we are going to simulate some data from it:
n=1000
set.seed(50) #reproducibility.
sim <- gamVineSimulate(n, GVC)
cdata = as.copuladata(sim)
pairs(cdata)
X1<-sim[,1]; X2<-sim[,2]; X3<-sim[,3];
sim <- data.frame(sim)
# first we plot the data
par(mfrow=c(1,1))
# plotting Clayton Copula
ggplot(sim, aes(x = sim[,1], y = sim[,2]))+geom_point()+
  labs(title = TeX('Clayton Copula with $\\theta$=2'), x = "X1",y = "X2")
# plotting Gumbell Copula
ggplot(sim, aes(x = sim[,2], y = sim[,3]))+geom_point()+
  labs(title = TeX('Gumbell Copula with $\\alpha$=1.5'), x = "X2",y = "X3")
# plotting the Margins x1,x3
ggplot(sim, aes(x = sim[,1], y = sim[,3]))+geom_point()+
  labs(title = TeX('Margins X1 vs X3 sinusoidal dep. on X2'), x = "X1",y = "X3")
# 3d plot
plot_ly(x=sim[,1], y=sim[,2], z=sim[,3], type="scatter3d", mode="markers",size=0.1)


logLik<-sum(log(gamVinePDF(GVC, sim)))
print(logLik)
# But we are also interested in visualising the copula of the second tree
# In order to do that, we have to compute psuedoobservations:
X1cond2<- BiCopHfunc1(X2,X1,3,2,0,obj=NULL,check.pars=TRUE)
X3cond2<- BiCopHfunc1(X2,X3,4,1.5,0,obj=NULL,check.pars=TRUE)
pseudoobs=data.frame(X1cond2,X3cond2)
#Plotting the psuedo observations:
ggplot(pseudoobs, aes(x=pseudoobs[,1], y=pseudoobs[,2]))+geom_point()+
labs(title = TeX('Pseudoobservations, Gaussian Copula'), 
     x = TeX('X_{1|2}'),y = TeX('X_{3|2}'))
cdata2= as.copuladata(pseudoobs)
pairs(cdata2)

#############################
# Then with this data we are going to fit several models and see their performances.
# from best to worst:
# knowing the whole structure:
system.time(fit1GVC <- gamVineSeqFit(sim, GVC, verbose = TRUE))
summary(fit1GVC)
logLik1<-sum(log(gamVinePDF(fit1GVC, sim)))
print(logLik1) # we obtain better loglik, but it is normal if we think about it

# knowing only the Matrix that defines the vine copula:
system.time(fit2GVC <- gamVineCopSelect(sim, Matrix, verbose = TRUE))
summary(fit2GVC)
logLik2<-sum(log(gamVinePDF(fit2GVC, sim)))
print(logLik2)

# knowing nothing
system.time(fit4GVC <- gamVineStructureSelect(sim, lin.covs = NULL, smooth.covs = NULL,
            simplified = FALSE,verbose= TRUE, type = 0,familyset <- c(1:2,301:304,401:404)))
summary(fit4GVC)
logLik4<-sum(log(gamVinePDF(fit4GVC, sim)))
print(logLik4)

# Now we use the simplified assumption to fit this data:
# If we use the copula families involved in the previous model, we get:
system.time(RVM1<-RVineStructureSelect(sim, familyset = c(1:4),
                type = 0, selectioncrit = "logLik",indeptest = FALSE,level = 0.05,
                trunclevel = NA,progress = FALSE,weights = NA,treecrit = "tau",
                rotations = TRUE,se = FALSE,presel = TRUE,method = "mle",cores = 1))
summary(RVM1)
contour(RVM1)
logLik5<-RVM1$logLik
print(logLik5)

# If we do not make any assumptions about the families that are presented:
system.time(RVM2<-RVineStructureSelect(sim,familyset = NA,type = 0,
                  selectioncrit = "logLik",indeptest = FALSE,level = 0.05,
                  trunclevel = NA,progress = FALSE,weights = NA, treecrit = "tau",
                  rotations = TRUE,se = FALSE,presel = TRUE,method = "mle",cores = 1))
summary(RVM2)
logLik6<-RVM2$logLik
print(logLik6)

# assuming that we know the structure, the results are far better:
system.time(RVM3<-RVineCopSelect(sim, familyset = c(1:4), Matrix, selectioncrit = "logLik",
               indeptest = FALSE,level = 0.05,trunclevel = NA,weights = NA,
               rotations = TRUE,se = FALSE,presel = TRUE,method = "mle",
               cores = 1))
summary(RVM3)
print(RVM3$logLik)

#########################################
# now we will proceed using some goodness of fit test.
# In particular we are going to use the 
# Probability Integral Transform that we studied in class.
data.test1<-RVinePIT(sim,RVM1)
m<-qnorm(data.test1[,1],0,1)^2+qnorm(data.test1[,2],0,1)^2+
  qnorm(data.test1[,3],0,1)^2
ks.test(m,'pchisq',3,alternative="two.sided")
# but we will doa robust analysis of it
# NOT RUN, NOT GOOD RESULTS.
p=c();
n=1000
for (i in 1:100){
set.seed(i) #reproducibility.
sim <- gamVineSimulate(n, GVC)
RVM1<-RVineStructureSelect(sim, familyset = c(1:4),
      type = 0, selectioncrit = "logLik",indeptest = FALSE,level = 0.05,
      trunclevel = NA,progress = FALSE,weights = NA,treecrit = "tau",
      rotations = TRUE,se = FALSE,presel = TRUE,method = "mle",cores = 1)
data.test1<-RVinePIT(sim,RVM1)
m<-qnorm(data.test1[,1],0,1)^2+qnorm(data.test1[,2],0,1)^2+
  qnorm(data.test1[,3],0,1)^2
mm<-ks.test(m,'pchisq',3,alternative="two.sided")
p[i]=mm$p.value
}
hist(p,breaks = 15)
length(which(p>0.05))/length(p)


# Moreover we will aply the pacotest to test simplified assumptions:
pacotestOptions=pacotestset(testType='CCC',  grouping ='TreeCCC', withEstUncert = FALSE,
              groupedScatterplots=TRUE, estUncertWithRanks = FALSE, decisionTreePlot=TRUE)
test=pacotestRvineSeq(sim, RVM3, pacotestOptions 
                 ,level = 0.05, illustration = 2, stopIfRejected = TRUE)
labels=c()
for(i in 1:1000){
  if(X2[i]<=0.0935){labels[i]='a'}
  if(X2[i]>0.0935 && X2[i]<=0.4841){labels[i]='b'}
  if(X2[i]>0.4841 && X2[i]<=0.621){labels[i]='c'}
  if(X2[i]>0.621){labels[i]='d'}
}
pseudoobs=data.frame(X1cond2,X3cond2,labels)

# all the pseudo observations
ggplot(pseudoobs, aes(x=pseudoobs[,1], y=pseudoobs[,2],color = labels)) +
  geom_point(show.legend = FALSE)+ scale_color_manual(values=c('#000000','#000000', '#000000','#000000'))
+labs(title = TeX('Pseudoobservations, Gaussian Copula'))

# 2nd tree:
ggplot(pseudoobs, aes(x=pseudoobs[,1], y=pseudoobs[,2],color = labels)) +
  geom_point(show.legend = FALSE)+ scale_color_manual(values=c('#000000','#000000', '#D3D3D3','#D3D3D3'))
+labs(title = TeX('Pseudoobservations, Gaussian Copula'))

ggplot(pseudoobs, aes(x=pseudoobs[,1], y=pseudoobs[,2],color = labels)) +
  geom_point(show.legend = FALSE)+ scale_color_manual(values=c('#D3D3D3','#D3D3D3', '#000000','#000000'))
+labs(title = TeX('Pseudoobservations, Gaussian Copula'))

#3rd tree:
ggplot(pseudoobs, aes(x=pseudoobs[,1], y=pseudoobs[,2],color = labels)) +
  geom_point(show.legend = FALSE)+ scale_color_manual(values=c('#000000','#D3D3D3', '#D3D3D3','#D3D3D3'))
+labs(title = TeX('Pseudoobservations, Gaussian Copula'))

ggplot(pseudoobs, aes(x=pseudoobs[,1], y=pseudoobs[,2],color = labels)) +
  geom_point(show.legend = FALSE)+ scale_color_manual(values=c('#D3D3D3','#000000', '#D3D3D3','#D3D3D3'))
+labs(title = TeX('Pseudoobservations, Gaussian Copula'))

ggplot(pseudoobs, aes(x=pseudoobs[,1], y=pseudoobs[,2],color = labels)) +
  geom_point(show.legend = FALSE)+ scale_color_manual(values=c('#D3D3D3','#D3D3D3','#000000','#D3D3D3'))
+labs(title = TeX('Pseudoobservations, Gaussian Copula'))

ggplot(pseudoobs, aes(x=pseudoobs[,1], y=pseudoobs[,2],color = labels)) +
  geom_point(show.legend = FALSE)+ scale_color_manual(values=c('#D3D3D3','#D3D3D3','#D3D3D3','#000000'))
+labs(title = TeX('Pseudoobservations, Gaussian Copula'))


# We are going to try to find dependence between kendalls and the value of X2
M=seq(0,1,length.out=11)
correl=c()
for(i in 1:length(M)-1){
  correl[i]=cor(pseudoobs[which(X2>M[i] & X2<M[i+1]),], method="kendall")[1,2]
}
plot(M[1:length(M)-1],correl)
correlations=data.frame(M[1:length(M)-1],correl)
ggplot(correlations, aes(x=M[1:length(M)-1], y=correl))+
  geom_point()+labs(title = TeX('Value of Kendalls tau depending on $X_{2}$'), x = "X2",y = "tau")


# We lastly make some visual test for the fitted models obtained previously:
n=1000
set.seed(50) #reproducibility.
simnonsimply<-gamVineSimulate(n, fit2GVC)
X1<-simnonsimply[,1]
X2<-simnonsimply[,2]
X3<-simnonsimply[,3]
X1cond2nonsimply<- BiCopHfunc1(X2,X1,3,2.16,0,obj=NULL,check.pars=TRUE)
X3cond2nonsimply<- BiCopHfunc1(X2,X3,4,1.51,0,obj=NULL,check.pars=TRUE)
pseudoobsnonsimply=data.frame(X1cond2nonsimply,X3cond2nonsimply)
ggplot(pseudoobsnonsimply, aes(x=pseudoobsnonsimply[,1], y=pseudoobsnonsimply[,2]))+
  geom_point()+labs(title = TeX('Pseudoobservations, Gaussian Copula'), 
       x = TeX('X_{1|2}'),y = TeX('X_{3|2}'))
cdata2nonsimply= as.copuladata(pseudoobsnonsimply)
pairs(cdata2nonsimply)


n=1000
set.seed(50) #reproducibility.
simsimply<-RVineSim(n, RVM3)
X1<-simsimply[,1]; X2<-simsimply[,2]; X3<-simsimply[,3]
X1cond2simply<- BiCopHfunc1(X2,X1,3,2.16,0,obj=NULL,check.pars=TRUE)
X3cond2simply<- BiCopHfunc1(X2,X3,4,1.51,0,obj=NULL,check.pars=TRUE)
pseudoobssimply=data.frame(X1cond2simply,X3cond2simply)
ggplot(pseudoobssimply, aes(x=pseudoobssimply[,1], y=pseudoobssimply[,2]))+geom_point()+
  labs(title = TeX('Pseudoobservations, t copula'), 
       x = TeX('X_{1|2}'),y = TeX('X_{3|2}'))
cdata2simply= as.copuladata(pseudoobssimply)
pairs(cdata2simply)


#####################################
# now we are going to check the effect of the correlationsbetween
# the copulas on the first tree:
# low correlated case:
d <- 3
# matrix to specify the vine copula model.
Matrix <- c(1,3,2,0,3,2,0,0,2)
Matrix <- matrix(Matrix,d,d)
nnames <- paste("X", 1:d, sep = "")

model <- vector(mode = "list", length = d*(d-1)/2)
sel <- seq(d,d^2-d, by = d)
# specifying the first copula -> Clayton.
model[[1]]$family <- 301
model[[1]]$par <- 0.25
model[[1]]$par2 <- 0
# specifying the second copula -> Gumbell.
model[[2]]$family <- 401
model[[2]]$par <- 1.1
model[[2]]$par2 <- 0

set.seed(50) # reproducibility
data <- data.frame(u1 = runif(1e2), u2 = runif(1e2), matrix(runif(1e2*d),1e2,d))

cond <- nnames[2]
tmpform <- paste("~",paste(paste("s(", cond, ", k=10, bs='cr')",sep = ""), collapse=" + "))
m <- 1e2
x <- matrix(seq(0,1,length.out=m), nrow = m, ncol = 1)
tmp.fct <- function(x) eta0+calib.surf[[2]](x) # we first proceed with the sinusoidal.
colnames(x) <- cond
y <- tmp.fct(x)
form <- as.formula(paste0("y", tmpform))
dd <- data.frame(y, x)
names(dd) <- c("y", cond)
b1 <- gam(form, data = dd)
plot(b1,main='fitted Sinusoidal function') 
tmp <- gamBiCopFit(data = data, formula = form, family = 1, n.iters = 1)$res
attr(tmp, "model")$coefficients <- coefficients(b1)
attr(tmp, "model")$smooth <- b1$smooth
attr(tmp, "family") <- 1
model[[3]] <- tmp

# we define our Gam Vine Copula:
GVClowcor <- gamVine(Matrix=Matrix,model = model,names=nnames)
summary(GVClowcor)

n=1000
set.seed(10) #reproducibility.
simlowcor <- gamVineSimulate(n, GVClowcor)
cdatalowcor = as.copuladata(simlowcor)
pairs(cdatalowcor)
X1<-simlowcor[,1]; X2<-simlowcor[,2]; X3<-simlowcor[,3];
simlowcor <- data.frame(simlowcor)
# plotting Clayton Copula
ggplot(simlowcor, aes(x = simlowcor[,1], y = simlowcor[,2]))+geom_point()+
  labs(title = TeX('Clayton Copula with $\\theta$=0.25'), x = "X1",y = "X2")
# plotting Gumbell Copula
ggplot(simlowcor, aes(x = simlowcor[,2], y = simlowcor[,3]))+geom_point()+
  labs(title = TeX('Gumbell Copula with $\\theta$=1.1'), x = "X2",y = "X3")
# plotting the Margins x1,x3
ggplot(simlowcor, aes(x = simlowcor[,1], y = simlowcor[,3]))+geom_point()+
  labs(title = TeX('Margins X1 vs X3 sinusoidal dep. on X2'), x = "X1",y = "X3")
logLik<-sum(log(gamVinePDF(GVClowcor, simlowcor)))
print(logLik)


X1cond2<- BiCopHfunc1(X2,X1,3,0.25,0,obj=NULL,check.pars=TRUE)
X3cond2<- BiCopHfunc1(X2,X3,4,1.1,0,obj=NULL,check.pars=TRUE)
pseudoobslowcor=data.frame(X1cond2,X3cond2)
#Plotting the psuedo observations:
ggplot(pseudoobslowcor, aes(x=pseudoobslowcor[,1], y=pseudoobslowcor[,2]))+geom_point()+
  labs(title = TeX('Pseudoobservations, Gaussian Copula'), 
       x = TeX('X_{1|2}'),y = TeX('X_{3|2}'))
cdatalowcor= as.copuladata(pseudoobslowcor)
pairs(cdatalowcor)

# fitting models

system.time(fitGVClowcor <- gamVineCopSelect(simlowcor, Matrix, verbose = TRUE))
summary(fitGVClowcor)
logLik2<-sum(log(gamVinePDF(fitGVClowcor, simlowcor)))
print(logLik2)


system.time(RVMlowcor<-RVineCopSelect(simlowcor, familyset = c(1:4), Matrix, selectioncrit = "logLik",
                                 indeptest = FALSE,level = 0.05,trunclevel = NA,weights = NA,
                                 rotations = TRUE,se = FALSE,presel = TRUE,method = "mle",
                                 cores = 1))
summary(RVMlowcor)
print(RVMlowcor$logLik)



#######################
#highly correlated case


d <- 3
# matrix to specify the vine copula model.
Matrix <- c(1,3,2,0,3,2,0,0,2)
Matrix <- matrix(Matrix,d,d)
nnames <- paste("X", 1:d, sep = "")

model <- vector(mode = "list", length = d*(d-1)/2)
sel <- seq(d,d^2-d, by = d)
# specifying the first copula -> Clayton.
model[[1]]$family <- 301
model[[1]]$par <- 8
model[[1]]$par2 <- 0
# specifying the second copula -> Gumbell.
model[[2]]$family <- 401
model[[2]]$par <- 10
model[[2]]$par2 <- 0

set.seed(50) # reproducibility
data <- data.frame(u1 = runif(1e2), u2 = runif(1e2), matrix(runif(1e2*d),1e2,d))

cond <- nnames[2]
tmpform <- paste("~",paste(paste("s(", cond, ", k=10, bs='cr')",sep = ""), collapse=" + "))
m <- 1e2
x <- matrix(seq(0,1,length.out=m), nrow = m, ncol = 1)
tmp.fct <- function(x) eta0+calib.surf[[2]](x) # we first proceed with the sinusoidal.
colnames(x) <- cond
y <- tmp.fct(x)
form <- as.formula(paste0("y", tmpform))
dd <- data.frame(y, x)
names(dd) <- c("y", cond)
b1 <- gam(form, data = dd)
plot(b1,main='fitted Sinusoidal function') 
tmp <- gamBiCopFit(data = data, formula = form, family = 1, n.iters = 1)$res
attr(tmp, "model")$coefficients <- coefficients(b1)
attr(tmp, "model")$smooth <- b1$smooth
attr(tmp, "family") <- 1
model[[3]] <- tmp

# we define our Gam Vine Copula:
GVChighcor <- gamVine(Matrix=Matrix,model = model,names=nnames)
summary(GVChighcor)


n=1000
set.seed(100) #reproducibility.
simhighcor <- gamVineSimulate(n, GVChighcor)
cdatahighcor = as.copuladata(simhighcor)
pairs(cdatahighcor)
X1<-simhighcor[,1]; X2<-simhighcor[,2]; X3<-simhighcor[,3];
simhighcor <- data.frame(simhighcor)
# plotting Clayton Copula
ggplot(simhighcor, aes(x = simhighcor[,1], y = simhighcor[,2]))+geom_point()+
  labs(title = TeX('Clayton Copula with $\\theta$=8'), x = "X1",y = "X2")
# plotting Gumbell Copula
ggplot(simhighcor, aes(x = simhighcor[,2], y = simhighcor[,3]))+geom_point()+
  labs(title = TeX('Gumbell Copula with $\\theta$=10'), x = "X2",y = "X3")
# plotting the Margins x1,x3
ggplot(simhighcor, aes(x = simhighcor[,1], y = simhighcor[,3]))+geom_point()+
  labs(title = TeX('Margins X1 vs X3 sinusoidal dep. on X2'), x = "X1",y = "X3")
logLik<-sum(log(gamVinePDF(GVChighcor, simhighcor)))
print(logLik)


X1cond2<- BiCopHfunc1(X2,X1,3,8,0,obj=NULL,check.pars=TRUE)
X3cond2<- BiCopHfunc1(X2,X3,4,10,0,obj=NULL,check.pars=TRUE)
pseudoobshighcor=data.frame(X1cond2,X3cond2)
#Plotting the psuedo observations:
ggplot(pseudoobshighcor, aes(x=pseudoobshighcor[,1], y=pseudoobshighcor[,2]))+geom_point()+
  labs(title = TeX('Pseudoobservations, Gaussian Copula'), 
       x = TeX('X_{1|2}'),y = TeX('X_{3|2}'))
cdatahighcor= as.copuladata(pseudoobshighcor)
pairs(cdatahighcor)


system.time(fitGVChighcor <- gamVineCopSelect(simhighcor, Matrix, verbose = TRUE))
summary(fitGVChighcor )
logLik2<-sum(log(gamVinePDF(fitGVChighcor, simhighcor)))
print(logLik2)


system.time(RVMhighcor<-RVineCopSelect(simhighcor, familyset = c(1:4), Matrix, selectioncrit = "logLik",
                                      indeptest = FALSE,level = 0.05,trunclevel = NA,weights = NA,
                                      rotations = TRUE,se = FALSE,presel = TRUE,method = "mle",
                                      cores = 1))
summary(RVMhighcor)
print(RVMhighcor$logLik)
