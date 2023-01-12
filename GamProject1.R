library('copula')
library('VineCopula')
library('mgcv')


# first we have to define the formula we are working with:
#smooth terms within gam model formulae:
s(..., k=-1,fx=FALSE,bs="tp",m=NA,by=NA,xt=NULL,id=NULL,sp=NULL,pc=NULL)
# Define tensor product smooths  in GAM formulae
te(..., k=NA,bs="cr",m=NA,d=NA,by=NA,fx=FALSE,
   np=TRUE,xt=NULL,id=NULL,sp=NULL,pc=NULL)
# choose.k Choosing the basis dimension,
# and checking the choice, when using penalized regression smoothers
# large enough to represent the underlying truth, small enough for comp. efficiency.
#tensor product interactions
ti(..., k=NA,bs="cr",m=NA,d=NA,by=NA,fx=FALSE,
   np=TRUE,xt=NULL,id=NULL,sp=NULL,mc=NULL,pc=NULL)
example: y ~ ti(x) + ti(z) + ti(v) + ti(x,z) + ti(x,v) + ti(z,v)+ti(x,z,v).
# byvariables
y ~ s(x,by=z)
# Linear functional terms
s(X,Z,by=L)
#random efffects:
s(...,bs="re")
# penalising the model
gam(y ~ X -1,paraPen=list(X=list(S1,S2)))



gam(formula,family=gaussian(),data=list(),weights=NULL,subset=NULL,
    na.action,offset=NULL,method="GCV.Cp",
    optimizer=c("outer","newton"),control=list(),scale=0,
    select=FALSE,knots=NULL,sp=NULL,min.sp=NULL,H=NULL,gamma=1,
    fit=TRUE,paraPen=NULL,G=NULL,in.out,drop.unused.levels=TRUE,
    drop.intercept=NULL,discrete=FALSE,...)

gamBiCop(family, model, par2 = 0, tau = TRUE)
# where model is an gamObject as return by the gam function from the mgcv package.




