#' @import MASS
#' @import Matrix
#' @import matrixcalc
#' @import corrplot
#' @import pdSpecEst
#' @import depthTools
NULL
#> NULL


#' @description
#' evol function generates modulus for AR2 processes
#'
#' @param ini_mod hello
#' @param incre hello
#' @param trials hello
#' @param n_indpt hello
#' @export
evol <- function(ini_mod, incre, trials, n_indpt){
  module <- matrix(0, n_indpt, trials)
  module[,1] <- ini_mod
  for (tri in 2:trials){
    module[,tri] <- module[,tri-1] + incre
  }
  return(module)
}


#' @description
#' tran_AR2 function transfers module to AR2 parameters
#'
#' @param mod hello
#' @param phase hello
tran_AR2 <- function(mod, phase){
  phi1 <- (1/mod)*cos(phase)*2
  phi2 <- -1/(mod^2)
  return(c(phi1, phi2))
}

#'
#' @description
#' Time series generation
#'
#' @param module hello
#' @param sigma hello
#' @param Sigma hello
#' @param M      hello
#' @param n_channel  hello
#' @param n_indpt    hello
#' @param trials    hello
#' @param n_time   hello
#' @export
generate_time <- function(module, sigma, Sigma, M, n_channel, n_indpt, trials, n_time){
  con_Y <- array(0, dim = c(n_channel, n_time, trials))
  con_State <- array(0,  dim = c(n_indpt, n_time, trials))
  for (tri in 1:trials){
    phi1 = tran_AR2(module[1, tri], (2/n_time)*2*pi) #delta
    phi2 = tran_AR2(module[2, tri], (6/n_time)*2*pi) #theta
    phi3 = tran_AR2(module[3, tri], (9/n_time)*2*pi) #alpha
    phi4 = tran_AR2(module[4, tri], (15/n_time)*2*pi) #beta
    phi5 = tran_AR2(module[5, tri], (32/n_time)*2*pi) #gamma

    s1 <- arima.sim(model = list(ar = phi1), n = n_time,
                    innov = rnorm(n_time) * sigma[tri])
    s2 <- arima.sim(model = list(ar = phi2), n = n_time,
                    innov = rnorm(n_time) * sigma[tri])
    s3 <- arima.sim(model = list(ar = phi3), n = n_time,
                    innov = rnorm(n_time) * sigma[tri])
    s4 <- arima.sim(model = list(ar = phi4), n = n_time,
                    innov = rnorm(n_time) * sigma[tri])
    s5 <- arima.sim(model = list(ar = phi5), n = n_time,
                    innov = rnorm(n_time) * sigma[tri])



    S <- rbind(s1/sd(s1), s2/sd(s2), s3/sd(s3), s4/sd(s4), s5/sd(s5))
    con_State[,,tri] <- S
    con_Y[,,tri] <- M %*% S + matrix(mvrnorm(n_time, rep(0,nrow(M)), Sigma), nrow(M), n_time, byrow = TRUE)
  }
  return(list("realization" = con_Y, "latent" = con_State))
}


#' @description
#' Theoretical spectrum
#'
#' @param phi  hello
#' @param sigma  hello
spec <- function(phi, sigma){
  temp = 2*phi[1]*(phi[2] - 1)*cos(2*pi*seq(-500:500)/1000) - 2*phi[2]*cos(4*pi*seq(-500:500)/1000)+1+sum(phi^2)
  return(sigma^2/temp)
}




#' @description
#' Calculate Sum on the diagonal
#'
#' @param L hello
#' @param n hello
CalculateSumOnDiag = function(L, n){
  sum=0
  for(k in 1:(n-1)){
    sum = sum + L[n,k]*Conj(L[n,k])
  }
  return(sum)
}


#' @description
#' Calculate sum off the diagonal
#'
#' @param L hello
#' @param i hello
#' @param j hello
CalculateSumOffDiag = function(L, i, j){
  sum=0
  if(j>1){
    for(k in 1:(j-1)){
      sum = sum + L[i,k]*Conj(L[j,k])
    }
  }
  return(sum)
}



#' @description
#' Compute Complex Cholesky
#'
#' @param A hello
ComplexCholesky = function(A){
  n = dim(A)[1]
  L=matrix(0,n,n)
  L[1,1] = sqrt(A[1,1])
  for(i in 2:n){
    for(j in 1:(i-1)){
      L[i,j] = (A[i,j] - CalculateSumOffDiag(L,i,j))/L[j,j]
    }
    L[i,i] = sqrt(A[i,i] - CalculateSumOnDiag(L,i))
  }
  return(L)
}


#' @description
#' Smoothing window
#'
#' @param rawperiod hello
#' @param L hello
SmoothingWindow=function(rawperiod, L){

  ## Input
  ## rawperiod	=	raw periodograms
  ## 			=     (1/T)(abs(fft(timeseries)))^2
  ## L  		= any positive integer; span is 2*L + 1
  ## Output
  ## smoothedcrossper 		= smoothed cross periodogram defined from (0, 2*pi]

  T = length(rawperiod);
  temp = rawperiod;
  tempre = Re(temp);
  tempim = Im(temp);

  realtemp = c(rev(tempre[2:(L+1)]), tempre[1:T], rev(tempre[(T-1):(T-L)]));
  imagtemp = c(-rev(tempim[2:(L+1)]), tempim[1:T], -rev(tempim[(T-1):(T-L)]));

  realsmoothedper = c(1:T);
  imagsmoothedper = c(1:T);

  for (k in c(1:T)){

    startindex = k;
    endindex =  startindex + 2*L;
    realsmoothedper[k] = mean(realtemp[startindex:endindex]);
    imagsmoothedper[k] = mean(imagtemp[startindex:endindex]);
  }

  smoothedcrossper = complex(real=realsmoothedper, imag=imagsmoothedper);
  return(smoothedcrossper)
}


#' @description
#' Estimate spectrum density at particular frequency band
#'
#' @param rawPdg hello
#' @param span hello
#' @param p hello
#' @param T hello
Estimatef.omega = function(rawPdg,span,p,T){
  smoothPdg =  array(0,dim=c(T,p,p))

  for(i in 1:p){
    for(j in 1:i){
      smoothPdg[,i,j] = SmoothingWindow(rawPdg[,i,j],span)
    }
  }

  return(smoothPdg)
}



#' @description
#' Compute Spectral estimate
#'
#' @param X.sim time series
#' @param Time time point
#' @param span spanning window size
#' @param n.taper Number of taper
#' @export
SpectrumEstimate = function(X.sim, Time, span, n.taper){
  p=dim(X.sim)[1];
  #J = floor((Time-K+1)/2)-1;
  J = floor((Time+1)/2)-1;

  H.taper = matrix(0,n.taper,Time)
  Time.seq=seq(-(Time/2 - 1),Time/2,1)
  for(i in 1:n.taper){
    for(j in 1:Time){
      #for(j in Time.seq){
      H.taper[i,j] = sqrt(2/(Time+1))*sin((pi*j*(i-1+1))/(Time+1))
    }
  }

  #### With tapper
  D = array(0,dim=c(n.taper,dim(X.sim)[1],Time))
  for(i in 1:n.taper){
    X.sim.taper = X.sim%*%diag(H.taper[i,])
    for(j in 1:p){
      D[i,j,] = fft(X.sim.taper[j,])
    }
  }


  I = array(0,dim = c(Time,dim(X.sim)[1],dim(X.sim)[1]))
  for(omega in 1:Time){
    for(i in 1:n.taper){
      temp = ((D[i,,omega])%*%Conj(t(D[i,,omega])))
      I[omega,,] = temp + I[omega,,]
    }
    I[omega,,] = I[omega,,]/n.taper
  }


  R = array(0,dim = c(Time,dim(X.sim)[1],dim(X.sim)[1]))
  for(omega in 1:Time){
    R[omega,,] = ComplexCholesky(I[omega,,])
  }

  ###Correct R
  temp = 1:p
  delta0 = sqrt(2)*gamma(n.taper/2 - temp/2 + 1)/(sqrt(n.taper)*gamma(n.taper/2-temp/2+1/2))
  Delta0 = diag(delta0)

  #delta = gamma(n.taper - temp + 3/2)/(sqrt(n.taper)*gamma(n.taper-temp+1))
  delta = gamma(n.taper - temp + 3/2)/(sqrt(n.taper)*gamma(n.taper-temp+1))
  Delta = diag(delta)
  for(omega in 1:Time){
    if((omega==1)|(omega==Time/2)|(omega==Time)){
      R[omega,,] = R[omega,,]%*%solve(Delta0)
    }else{
      R[omega,,] = R[omega,,]%*%solve(Delta)
    }
  }

  ##Smothing Cholesky Components
  R.smooth = Estimatef.omega(R,span,p,Time)
  f.smooth = array(0,dim=c(Time,p,p))
  for(i in 1:Time){
    f.smooth[i,,] = R.smooth[i,,]%*%Conj(t(R.smooth[i,,]))
  }

  return(f.smooth)
}


#' @description
#' Calculate the spectrum power at particular frequency band
#'
#' @param X spectrum power
#' @param freqband frequency band
#' @export
CoherenceEstimate = function(X, freqband){
  #dim(X) = Time, p, p
  p = dim(X)[2]
  Time = dim(X)[1]
  result= matrix(1,p,p)

  temp=0
  for(j in freqband){
    temp = temp + X[j,,]
  }

  result = temp/(length(freqband))

  coh = matrix(1,p,p)
  for(i in 2:p){
    for(j in 1:(i-1)){
      coh[i,j] = (abs(result[i,j]))^2/(Re(result[i,i])*Re(result[j,j]))
      ##Fisher transformation does not give positive definite.
      ##coh[i,j,] = log((1 + coh[i,j,])/(1 - coh[i,j,]))/2
      coh[j,i] = coh[i,j]
    }
  }
  ##Power at particular frequency band
  return(coh)
}



#' @description
#' Combinat function
#'
#' @param n hello
#' @param p hello
combinat=function(n,p){
  if (n<p){combinat=0}
  else {combinat=exp(lfactorial(n)-(lfactorial(p)+lfactorial(n-p)))}
}


#' @description
#' Compute GMVD using log Euclidean metric
#'
#' @param data hello
fGSB=function(data){
  pdDepth(X=data, metric = 'logEuclidean')
}



#' @description
#' Compute MVD using Euclidean metric
#'
#' @param data hello
fMBD=function(data){
  p=dim(data)[1]
  n=dim(data)[2]
  rmat=apply(data,1,rank)
  down=rmat-1
  up=n-rmat
  (rowSums(up*down)/p+n-1)/combinat(n,2)
}

#surface boxplot
#fit: p1 by p2 by n data array, n is the number of surfaces
#method: BD2, MBD
#' @description
#' Geometric surface boxplot
#'
#' @param fit p1 by p2 by n data array, n is the number of surfaces
#' @param method metric used for constructing surface boxplot.
#' @export
GSBplot=function(fit,grid=NULL,method='GSB',factor=1.5,...){
  size=dim(fit)
  p1=size[1]
  p2=size[2]
  n=size[3]
  tp=p1*p2

  dmat=matrix(as.vector(fit),tp,n)

  if (method=='GSB') {depth=fGSB(fit)}
  else if (method=='SB') {depth=fMBD(dmat)}

  dp_s=sort(depth,decreasing=T)
  index=order(depth,decreasing=T)
  med=index[1]
  prob=0.5
  m=ceiling(n*prob)#at least 50%
  center=dmat[,index[1:m]]
  out=dmat[,index[(m+1):n]]
  infenv=apply(center,1,order)[1,]
  supenv=apply(center,1,order)[m,]

  inf=apply(center,1,min)
  sup=apply(center,1,max)

  dist=factor*(sup-inf)
  upper=sup+dist
  lower=inf-dist
  outly=(dmat<lower)+(dmat>upper)
  outcol=colSums(outly)
  remove=(outcol>0)
  ##outlier column
  colum=1:n
  outpoint=colum[remove==1]
  woout=dmat
  good=dmat[,(remove==0),drop=FALSE]
  mm=dim(good)[2]
  maxenv=apply(good,1,order)[mm,]
  minenv=apply(good,1,order)[1,]

  return(list(depth=depth,outpoint=outpoint,medindex=med,indexOrder=index,
              centerupper=matrix(supenv,p1,p2),centerlower=matrix(infenv,p1,p2),
              maxenv=matrix(maxenv,p1,p2),minenv=matrix(minenv,p1,p2)))
}

#' @description
#' Compute Modified volume depth
#'
#' @param x hello
#' @param xreff hello
#' @param method hello
ModifiedVolumeDepth = function(x,xreff,method='MBD',factor=1.5){
  n1 = dim(x)[3]
  n2 = dim(xreff)[3]
  p = dim(x)[1]
  N=n1+n2
  x.merge = array(0, dim = c(p,p,N))
  x.merge[,,1:n1] = x[,,1:n1]
  x.merge[,,(n1+1):N] = xreff[,,1:n2]
  result = GSBplot(x.merge, method='SB')
  location = match(1:n1,result$indexOrder)
  x.depth = result$depth[location]
  return(x.depth)
}

#' @description
#' Pointwise rank sum test using MVD metric
#'
#' @param x the first sample of matrices
#' @param y the second sample of matrices
#' @param n sample size of the first sample
#' @param m sample size of the second sample
PointwiseRankSumTest = function (x, y, n, m){
  set.seed(333)
  p <- dim(x)[1]
  N <- dim(x)[3]
  M <- dim(y)[3]
  u <- matrix(runif(N), 1, N)
  w <- matrix(runif(M), 1, M)
  I <- order(u)
  J <- order(w)
  n0 <- max(N - n, M - m)
  if (n0 <= max(n, m)) {
    stop("Incorrect sample sizes")
  }
  g1 <- x[,,I[1:n]]
  g2 <- y[,,J[1:m]]
  if ((N - n) >= (M - m)) {
    gref <- x[,,I[(n + 1):N]]
  }else{
    gref <- y[,,J[(m + 1):M]]
  }
  r <- ModifiedVolumeDepth(g1, gref)
  s <- ModifiedVolumeDepth(g2, gref)
  z <- GSBplot(gref, method = 'SB')$depth
  n1 <- length(r)
  n2 <- length(s)
  g <- cbind(r, s)
  I <- order(g)
  u <- which(I > n1)
  w <- which(I <= n1)
  p.value <- wilcox.test(u, w)$p.value
  Rx <- c()
  Ry <- c()
  for (i in 1:n) {
    Rx[i] <- sum(z <= r[i])/length(z)
  }
  for (i in 1:m) {
    Ry[i] <- sum(z <= s[i])/length(z)
  }
  R <- order(c(Rx, Ry))
  W <- sum(R[(n + 1):(n + m)])
  return(list(p.value = p.value, statistic = W))
}

#' @description
#' Plot the coherence matrix
#'
#' @param cohX the coherence matrix
#' @param col1 color map
#' @export
plotCohereneMatrix = function(cohX, col1){
  corrplot::corrplot(cohX, method="shade",col=col1(100),
           col.lim=c(0,1),title = "",mar=c(0,0,1,0),
           tl.cex = 2, cl.cex = 2)
}

#' @description
#' Rank sum test for two independent samples of SPD matrices
#'
#' @param cohX the first sample
#' @param cohY the second sample
#' @param nchannel a number of channel
#' @param trials a number of trials/epoch
#' @export
rankSumTest = function(cohX, cohY, n_channel, trials){
  result = pdSpecEst::pdRankTests(array(c(cohX, cohY), dim=c(n_channel,n_channel,2*trials)), sample_sizes = c(trials,trials), test = 'rank.sum', metric = 'logEuclidean')
  return(list(p.value = result$p.value,
              test.stat = result$statistic,
              depth = result$depth.values))
}


#' @description
#' Rank signed test for two dependent samples of SPD matrices
#'
#' @param cohX the first sample
#' @param cohY the second sample
#' @param nchannel a number of channel
#' @param trials a number of trials/epoch
#' @export
rankSignedTest = function(cohX, cohY, n_channel, trials){
  result = pdSpecEst::pdRankTests(array(c(cohX, cohY), dim=c(n_channel,n_channel,2*trials)), test = 'signed.rank', metric = 'logEuclidean')
  return(list(p.value = result$p.value,
              test.stat = result$statistic,
              depth = result$depth.values))
}


#' @description
#' Pointwise rank sum test using MVD metric
#' @param cohX the first sample
#' @param cohY the second sample
#' @param n1 initial the first sample size
#' @param n2 initial the second sample size
#' @export
rankSumTestMVD = function(cohX, cohY, n1 = 40, n2 = 40){
  result = PointwiseRankSumTest(cohX, cohY, n1, n2)
  return(result)
}



