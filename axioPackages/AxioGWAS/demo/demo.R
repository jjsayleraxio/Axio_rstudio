#######################
#This is the power function for continuous outcome to detect interaction term b3
#Y= b0+b1Z+b2G+b3G*Z
#two arm: intervention and placebo
#N is total sample size
#r is proportion in intervention
#ey0, ey1 is the observed outcome mean in the control arm and in the intervention arm
#sy0, sy1 is the observed outcome SD in the control arm and in the intervention arm
#p is the minor allele frequency
###################
################

ey1<-1.09
ey0<-0.84
sy1<-1.29
sy0<-0.78
p<-0.3

N<-372

getPower(sy1,sy0, p, b2=0,b3=0.3, alpha=0.05, N=1000, r=0.5)
getSample(sy1,sy0, p, b2=0,b3=0.3, alpha=0.05, power=0.829, r=0.5)

###############################

N00<-400
N01<-100
N10<-300
N11<-200

###Power curve for MAF
p.mat<-matrix(seq(0.05,0.5,by=0.005),ncol=1)
power.pmat<-apply(p.mat,1,getPowerCC,ORg=1,ORgxt=1.5,N00=N00,N01=N01,N10=N10,N11=N11,alpha=0.01)

###Power curve for Interaction OR
OR.mat<-matrix(seq(1,2.5,by=0.05),ncol=1)
power.omat<-apply(OR.mat,1,getPowerCC,ORg=1,N00=N00,N01=N01,N10=N10,N11=N11,p=0.3,alpha=0.01)

###sensitivity test of power accross a range of genetic effect in control
sen.mat<-matrix(seq(0.5,2.5,by=0.05),ncol=1)
power.senmat<-apply(sen.mat,1,getPowerCC,ORgxt=1.5,N00=N00,N01=N01,N10=N10,N11=N11,p=0.3,alpha=0.00)


plot(OR.mat[,1],power.omat, xlab="Per allele OR in trt arm/Per allele OR in control arm", ylab="power")
plot(p.mat[,1],power.pmat, xlab="MAF", ylab="power")
plot(sen.mat[,1],power.senmat, xlab="Per allele OR in control arm", ylab="power")

################################

Pow.or <- computePowerCohort(1-(1-0.2)^2, 0.1, seq(from=1.1, to=1.8, by=0.025), c(1000, 2000, 5000), 10^(-6))
plot(Pow.or$OR, Pow.or[,2], type="l", ylim=c(0,1), xlab="Odds Ratio", ylab="Power")
for(i in 3:(dim(Pow.or)[2])){
	lines(Pow.or$OR, Pow.or[,i], lty=i)
}

N<-cbind(c(500, 1000, 1000), c(1000, 1000, 2000))
Pow.or<-computePowerCC(1-(1-0.2)^2, 0.1, seq(from=1.1, to=1.8, by=0.025), N, 10^(-6))
plot(Pow.or$OR, Pow.or[,2], type="l", ylim=c(0,1), xlab="Odds Ratio", ylab="Power")
for(i in 3:(dim(Pow.or)[2])){
	lines(Pow.or$OR, Pow.or[,i], lty=i)
}
