########################################LOAD DATA######################################
library(mvtnorm)
source("./load.data.R")

gen<-read.csv(file="./HuyangMap-Genotype.csv", header=T);
phe1<-read.csv("./OCK_length.csv");
phe2<-read.csv("./Salt_length.csv");
phe3<-read.csv("./OCK_number.csv");
phe4<-read.csv("./Salt_number.csv");
times<-seq(13,78,by=5)


dat1<-load_data(datt=gen,pheno=phe1[seq(1,dim(phe1)[1]-2,by=3),],nstart=1,nlen=8305,times)
dat2<-load_data(datt=gen,pheno=phe1[seq(2,dim(phe1)[1]-1,by=3),],nstart=1,nlen=8305,times)

sdat1<-load_data(datt=gen,pheno=phe2[seq(1,dim(phe2)[1]-2,by=3),],nstart=1,nlen=8305,times)
sdat2<-load_data(datt=gen,pheno=phe2[seq(2,dim(phe2)[1]-1,by=3),],nstart=1,nlen=8305,times)



##########################################FITTING##########################################
library(mvtnorm)
library(deSolve)
source("./IE.ode.R")



#s.par<-c(31.36955876 ,0.08151995,-0.27062018,0.02151682,0.05686985,
#         19.01929746,0.07426594,0.11060744,2.38008408,-0.16858089)#0.8597233
s.par1<-c(31.36955876 ,0.08151995,-0.0007,0.02151682,0.05686985,
         19.01929746,0.07426594,0.00043,2.38008408,-0.16858089)

H0 <- get_est_param(dat1,dat2,s.par,
                    x1=colMeans(dat1$pheno)[1],x2=colMeans(dat2$pheno)[1])
RSS <- H0[1]
predict<-ind.get_mu(par=H0[-1], times=dat1$sample_times,x1=colMeans(dat1$pheno)[1],x2=colMeans(dat2$pheno)[1])
TSS<-sum((predict-mean(predict))^2)
R_squard<-1-RSS/TSS
RSS
R_squard

#s.par<-c(33.49258284,0.09225487,-1.01534770,0.02208013, 0.03332821,
#         22.49663161,0.03682584,-0.41867501,0.89444496,-0.07420519)#0.55308552
s.par2<-c(33.49258284,0.09225487,-0.002796,0.02208013, 0.03332821,
         22.49663161,0.03682584,-0.000685,0.89444496,-0.07420519)
H0 <- get_est_param(sdat1,sdat2,s.par,
                    x1=colMeans(sdat1$pheno)[1],x2=colMeans(sdat2$pheno)[1])
RSS <- H0[1]
predict<-ind.get_mu(par=H0[-1], times=sdat1$sample_times,x1=colMeans(sdat1$pheno)[1],x2=colMeans(sdat2$pheno)[1])
TSS<-sum((predict-mean(predict))^2)
R_squard<-1-RSS/TSS
RSS
R_squard


##################################################################################
library(deSolve)
library(mvtnorm)

source("IE.ode.R")
source("IE.covar.R")
source("IE.curve.R")
source("IE.optim.R")


init.par<-IE.H0(dat1,dat2,parin2=c(0.5,13.1850,0.5,6.6132,0.3,s.par))
ret1 <- IE.est1(dat1,dat2,
                parin2=c(0.5,13.1850,0.5,6.6132,0.3,s.par),interval = c(1,20))

init.par<-IE.H0(sdat1,sdat2,parin2=c(0.5,7.2421,0.5,4.0028,0.3,s.par))
ret1 <- IE.est1(sdat1,sdat2,
                   parin2=c(0.5,7.2421,0.5,4.0028,0.3,s.par),interval = c(1,20))




