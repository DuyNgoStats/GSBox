devtools::install_local("GSBox_0.1.0.tar.gz")

library(GSBox)
library(LSTS)
##color scale for coherence
col1 <- colorRampPalette(c("white","blue","blue","yellow", "red"))

################################################################
################################################################
## Time point
n_time = 200
## Number of trials  
trials = 100
## Numer of EEG channels
n_channel = 15
## Numer of frequency band.
n_indpt = 5   
## Number of replicated data
n.sim = 100
## Compute the AR2 coefficient 
module = evol(rep(1.001, 5), 0.00001, trials, n_indpt)

## A multivariate time series (MTS) X(t) = [X_1(t), X_2(t), X_3(t)]
## Each cluster X_i(t) have nM elements.
nM = 5
## vector of zero
m0 = rep(0, nM)
## Mixing matrix for the MTV X(t)
M = rbind(cbind(rep(1, nM),  rep(0.1, nM),   m0,   m0,   m0),
          cbind(rep(1, nM),  rep(0.1, nM),   m0,     m0,    m0),
          cbind(m0, rep(0.1, nM),  rep(0.5, nM),   m0,  m0))

##### Mixing matrix for outliers
M1 = rbind(cbind(rep(0.3, nM),  rep(0.1, nM),   m0,   m0,   m0),
           cbind(rep(1, nM),  rep(0.5, nM),   m0,     m0,    m0),
           cbind(rep(1, nM), m0,  rep(0.5, nM),   m0,  m0))

## Design covariance matrix
rho12=0; rho13=0; rho23=0; r1= 0.01; r2= 0.01; r3=0.01
R1 = matrix(r1, 5, 5)
R2 = matrix(r2, 5, 5)
R3 = matrix(r3, 5, 5)
R12 = matrix(rho12, 5, 5)
R13 = matrix(rho13, 5, 5)
R23 = matrix(rho23, 5, 5)
Sigma = rbind(cbind(R1, R12, R13),
              cbind(R12, R2, R23),
              cbind(R13, R23, R3))
diag(Sigma) = 1

## Generate MTV
X = generate_time(module, sigma=rep(.05, trials), Sigma, 
                  M, n_channel, n_indpt, trials, n_time)
## Generate outlying at the last 5 trials 
X.out = generate_time(module, sigma=rep(.01, 5), Sigma, 
                      M1, n_channel, n_indpt, 5, n_time)$realization
#################################################################
### Plot the true TS and powers of latent sources and 
### corresponding MTS X(t)
#################################################################
## Plot true latent sources
X.latent = X$latent[c(1,2,3),,5]

plot(X.latent[1,1:100], type="l", xlab = "Time", ylab="",
     cex.lab=1.6, cex.axis=1.7, cex.main=1.7, cex.sub=1.7,
     lwd=2, ylim = c(min(X.latent)-0.1, max(X.latent)+0.1))
lines(X.latent[2,1:100], lty=2, lwd=2)
lines(X.latent[3,1:100], lty=3, lwd=2)
legend('top', legend=c(expression(paste(Z[1],'(t', ',', 'r)')),
                       expression(paste(Z[2],'(t', ',', 'r)')),
                       expression(paste(Z[3],'(t', ',', 'r)'))),
       horiz=TRUE, inset=c(0, -0.175), xpd=TRUE,
       lty=1:3, lwd = 4, cex=1.7,title = '', text.font=4, box.lty=0)


## Plot the spectrum powers of latent sources
peri_true = matrix(0, 3, 200)
for(i in 1:3){
  peri_true[i,]= scale(abs(fft(X.latent[i,])/sqrt(200))^2)
}
plot(peri_true[1,1:30], type="l", xlab = "Frequency in Hz", ylab="Power",
     cex.lab=1.6, cex.axis=1.7, cex.main=1.7, cex.sub=1.7,
     lwd=2, ylim = c(0, max(peri_true)+0.5))
lines(peri_true[2,1:30], lty=2, lwd=2)
lines(peri_true[3,1:30], lty=3, lwd=2)
legend(20, 10, legend=c(expression(paste(Z[1],'(t', ',', 'r)')),
                        expression(paste(Z[2],'(t', ',', 'r)')),
                        expression(paste(Z[3],'(t', ',', 'r)'))),
       lty=1:3, lwd=4, cex=1.9,title = 'Latents', text.font=4, box.lty=0)

## Smooth periodogram for MTS X(t) at trial 1
X.ts = X$realization[,,1]
f.smooth = matrix(0, dim(X.ts)[1], dim(X.ts)[2])
for(i in 1:dim(X.ts)[1]){
  f.smooth[i, ]= smooth.periodogram(X.ts[i,], plot = FALSE, spar = 0.5)$smooth.periodogram
}

plot(f.smooth[1,1:30], type="l", xlab = "Frequency in Hz", ylab="Power",
     cex.lab=1.6, cex.axis=1.7, cex.main=1.7, cex.sub=1.7,
     lwd=2, ylim = c(0, max(f.smooth)+0.5))
for(i in 2:15){
  if(i <=5 )  lines(f.smooth[i,1:30], lty=1, lwd=2)
  if(i>5 & i <=10)  lines(f.smooth[i,1:50], lty=2, lwd=2)
  if(i > 10 )  lines(f.smooth[i,1:30], lty=3, lwd=2)
}
legend(20, max(f.smooth), legend=c(expression(paste(X[1],'(t', ',', 'r)')),
                       expression(paste(X[2],'(t', ',', 'r)')),
                       expression(paste(X[3],'(t', ',', 'r)'))),
       lty=1:3, lwd=4, cex=1.9,title = 'Clusters', text.font=4, box.lty=0)


#################################################################
### Plot coherence matrix
#################################################################
## Compute the coherence matrix at delta band 
cohX = array(0,dim = c(n_channel, n_channel, trials+5))
for(j in 1:trials){
  X.sim = X$realization[,,j]
  f.smooth = SpectrumEstimate(X.sim, n_time, span=30, n_channel)
  cohX[,,j] = CoherenceEstimate(f.smooth,c(1:3))
}
## Plot the coherence matrix at the 1st trial
plotCohereneMatrix(cohX[,,1], col1)

## Adding outlying coherence matrix at the last 5 trials
for(j in (trials + 1): (trials + 5)){
  X.sim = X.out[,,j-100]
  f.smooth = SpectrumEstimate(X.sim, n_time, span=30, n_channel)
  cohX[,,j] = CoherenceEstimate(f.smooth,c(1:3))
}
## Plot the outlying coherence matrix 
plotCohereneMatrix(cohX[,,101], col1)

## Perform GSBox
out.s = GSBplot(cohX, method = 'GSB')
## The median coherence matrix
out.s$medindex
plotCohereneMatrix(cohX[,,out.s$medindex], col1)
## The outlying coherence matrices
out.s$outpoint
plotCohereneMatrix(cohX[,,out.s$outpoint[1]], col1)



