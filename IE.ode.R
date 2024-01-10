get_est_param<-function(dat1,dat2,s.par,x1,x2)
{
  allpheno <- cbind(dat1$pheno,dat2$pheno)
  allpheno[which(is.na(allpheno)),1]<-0
  times <- dat1$sample_times
  
  mpheno <- as.numeric(colMeans(allpheno, na.rm=TRUE))
  
  res <- optim(s.par,ind.mle,s.y=mpheno,s.t=times,x1=x1,x2=x2,
               method="BFGS",control=list(trace=T,maxit=1000))
  value <- res$value
  allpar <- res$par
  
  return(c(value,allpar));
}

ind.mle <- function(s.par,s.y,s.t,x1,x2){
  A <- sum((s.y - ind.get_mu(s.par,s.t,x1,x2))^2 )
  A
}


ind.get_mu <- function(par, times,x1,x2)
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      K1 = par[1],
      r1 = par[2],
      a  = par[3],
      c1 = par[4],
      d1 = par[5],
      K2 = par[6],
      r2 = par[7],
      b =  par[8],
      c2 = par[9],
      d2 = par[10]);
  }
  
  state0 <- c(X=x1, Y=x2);
  y <- COMP.f( par0, state0, times );
  
  
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:3] ) );
}


ind.get_mu1 <- function(par, times,x1,x2)
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      K1 = par[1],
      r1 = par[2],
      K2 = par[6],
      r2 = par[7]);
  }
  
  state0 <- c(X=x1, Y=x2);
  y <- COMP.f1( par0, state0, times );
  
  
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:3] ) );
}


ind.get_mu2 <- function(par, times,x1,x2)
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      K1 = par[1],
      r1 = par[2],
      a  = par[3],
      K2 = par[6],
      r2 = par[7],
      b =  par[8]);
  }
  
  state0 <- c(X=x1, Y=x2);
  y <- COMP.f2( par0, state0, times );
  
  
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:3] ) );
}


COMP.f <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            
            dX <- r1*X*(1-(X/K1))+a*X*Y+c1*exp(d1*t)
            dY <- r2*Y*(1-(Y/K2))+b*X*Y+c2*exp(d2*t)
            #dX <- r1*X*(1-(X/K1))+r1*a*X*Y/K1+c1*exp(d1*t)
            #dY <- r2*Y*(1-(Y/K2))+r2*b*X*Y/K2+c2*exp(d2*t)
            
            list(c(dX, dY))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}

COMP.f1 <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX <- r1*X*(1-(X/K1))
            dY <- r2*Y*(1-(Y/K2))
            
            list(c(dX, dY))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}



COMP.f2 <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            
            dX <- r1*X*(1-(X/K1))+a*X*Y
            dY <- r2*Y*(1-(Y/K2))+b*X*Y
            
            list(c(dX, dY))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}
