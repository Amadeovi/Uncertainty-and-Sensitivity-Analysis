n=1000
set.seed(7) #reproducibility.
sim <- gamVineSimulate(n, GVC)
X1<-sim[,1]; X2<-sim[,2]; X3<-sim[,3];
sim <- data.frame(sim)
X1cond2<- BiCopHfunc1(X2,X1,3,2,0,obj=NULL,check.pars=TRUE)
X3cond2<- BiCopHfunc1(X2,X3,4,1.5,0,obj=NULL,check.pars=TRUE)
pseudoobs=data.frame(X1cond2,X3cond2)
#Plotting the psuedo observations:
ggplot(pseudoobs, aes(x=pseudoobs[,1], y=pseudoobs[,2]))+geom_point()+
  labs(title = TeX('Pseudoobservations, Gaussian Copula'), 
       x = TeX('X_{1|2}'),y = TeX('X_{3|2}'))
cdata2= as.copuladata(pseudoobs)
pairs(cdata2)


