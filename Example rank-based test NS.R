devtools::install_local("GSBox_0.1.0.tar.gz")

library(GSBox)
################################################################
########### Rank-based tests for high dimension multivariate time
########### series.
################################################################
## Time point
n_time = 200
## Number of trials
trials = 100
## Numer of EEG channels
n_channel = 15
## Numer of frequency band.
## delta 1-3Hz, theta 4-8Hz, alpha 9-12Hz, beta 13-30Hz, gamma 30-50Hz.
n_indpt = 5   
## Number of replicated data
n.sim = 100
## Compute the AR2 coefficient 
module = evol(rep(1.001, 5), 0.00001, trials, n_indpt)

## A multivariate time series (MTS) X(t) = [X_1(t), X_2(t), X_3(t)]
## Each cluster X_i(t) and Y_i(t) have nM elements.
nM = 5
m0 = rep(0, nM) 
nM2 = 5
m02 = rep(0, nM2) 
##### Mixing matrix M1 (Equal case)
M11 = rbind(cbind(rep(1, nM),  rep(0.2, nM),   m0,   m0,   m0),
            cbind(rep(1, nM),  rep(0.5, nM),   m0,     m0,    m0),
            cbind(rep(0.1, nM2), rep(0.2, nM2),  rep(0.5, nM2),   m02,  m02))
M12 = rbind(cbind(rep(1, nM),  rep(0.2, nM),   m0,   m0,   m0),
            cbind(rep(1, nM),  rep(0.5, nM),   m0,     m0,    m0),
            cbind(rep(0.1, nM2), rep(0.2, nM2),  rep(0.5, nM2),   m02,  m02))
M1 = rbind(M11, M12)


##### Mixing matrix M2 (Non-equal case)
M21 = rbind(cbind(rep(1, nM),  rep(0.2, nM),   m0,   m0,   m0),
            cbind(rep(1, nM),  rep(0.5, nM),   m0,     m0,    m0),
            cbind(rep(0.1, nM2), rep(0.2, nM2),  rep(0.5, nM2),   m02,  m02))
M22 = rbind(cbind(rep(1, nM),  rep(0.2, nM),   m0,   m0,   m0),
            cbind(rep(0.2, nM),  m0,   m0,     rep(0.2, nM),    m0),
            cbind(m02, m02,  rep(0.5, nM2),   rep(0.2, nM2),  m02))
M2 = rbind(M21, M22)


####### Abrupt change
M31 = rbind(cbind(rep(0.5, nM),  rep(0.2, nM),   m0,   m0,   m0),
            cbind(rep(0.5, nM),  rep(0.5, nM),   m0,     m0,    m0),
            cbind(rep(0.2, nM2), rep(0.2, nM2),  rep(0.2, nM2),   m02,  m02))
M32 = rbind(cbind(rep(0.5, nM),  rep(0.2, nM),   m0,   m0,   m0),
            cbind(rep(0.5, nM),  rep(0.5, nM),   m0,     m0,    m0),
            cbind(rep(0.2, nM2), rep(0.2, nM2),  rep(0.2, nM2),   m02,  m02))
M3 = rbind(M31, M32)

## Design covariance matrix
Sigma1 = matrix(0.55, n_channel, n_channel)
Sigma2 = matrix(0.7, n_channel, n_channel)
## Rho can only take values [0, 0.6]
## 0: independent between X and Y.
Rho = matrix(0.4, n_channel, n_channel)
Sigma = rbind(cbind(Sigma1, Rho),
              cbind(Rho, Sigma2))
diag(Sigma) = 1

################################################################
################ 2 MTV X and Y have similar latent sources
################################################################
X.ts1 = generate_time(module, sigma=rep(.05, trials), Sigma, M1, 2*n_channel, 
                       n_indpt, trials, n_time)$realization

X.ts2 = generate_time(module, sigma=rep(.05, trials), Sigma, M3, 2*n_channel, 
                       n_indpt, trials, n_time)$realization

X.ts =  abind(X.ts1[,1:100,],X.ts2[,101:200,],along=2)

######## Plot non-stationary X(t)
X.sim = X.ts[c(1, 20), , 1]
plot(X.sim[1,], type="l", xlab = "Time", ylab="", 
     cex.lab=1.6, cex.axis=1.7, cex.main=1.7, cex.sub=1.7, 
     lwd=2, ylim = c(min(X.sim)-0.5, max(X.sim)+0.5))
lines(X.sim[2,], lwd=2, lty = 2)
legend('top', legend=c( expression(paste('X'^(1),'(t, 1)')), 
                        expression(paste('X'^(2),'(t, 1)'))), 
       horiz=TRUE, inset=c(0, -0.17), xpd=TRUE,
       lty=1:2, cex=1.55,title = '', text.font=2, box.lty=0)

acf(X.sim[1,],cex.lab=1.6, cex.axis=1.7, 
    main = '', 
    cex.main=1.7, cex.sub=1.7,lwd=2)

####### Perform tests on two samples of coherence matrices
cohX = array(0,dim = c(n_channel, n_channel, trials))
for(j in 1:trials){
  X.sim = X.ts[1:15,,j]
  f.smooth = SpectrumEstimate(X.sim, n_time, span=30, n_channel)
  cohX[,,j] = CoherenceEstimate(f.smooth,c(1:3))
}

cohY = array(0,dim = c(n_channel, n_channel, trials))
for(j in 1:trials){
  X.sim = X.ts[16:30,,j]
  f.smooth = SpectrumEstimate(X.sim, n_time, span=30, n_channel)
  cohY[,,j] = CoherenceEstimate(f.smooth,c(1:3))
}

result = matrix(0, 3, 2)
rownames(result) = c('Rank Signed Test', 'Rank sum Test',  'Rank Sum Test MVD')
colnames(result) = c('p-value', 'test statistics')
rsdt = rankSignedTest(cohX, cohY, n_channel, trials)
rst = rankSumTest(cohX, cohY, n_channel, trials)
rstmvd = rankSumTestMVD(cohX, cohY)
result[, 1] = c(rsdt$p.value, rst$p.value, rstmvd$p.value)
result[, 2] = c(rsdt$test.stat, rst$test.stat, rstmvd$statistic)
round(result, 2)

################################################################
################ 2 MTV X and Y have different latent sources
################################################################
X.ts1 = generate_time(module, sigma=rep(.05, trials), Sigma, M2, 2*n_channel, 
                       n_indpt, trials, n_time)$realization

X.ts2 = generate_time(module, sigma=rep(.05, trials), Sigma, M3, 2*n_channel, 
                       n_indpt, trials, n_time)$realization

X.ts =  abind(X.ts1[,1:100,],X.ts2[,101:200,],along=2)

cohX = array(0,dim = c(n_channel, n_channel, trials))
for(j in 1:trials){
  X.sim = X.ts[1:15,,j]
  f.smooth = SpectrumEstimate(X.sim, n_time, span=30, n_channel)
  cohX[,,j] = CoherenceEstimate(f.smooth,c(1:3))
}

cohY = array(0,dim = c(n_channel, n_channel, trials))
for(j in 1:trials){
  Y.sim = X.ts[16:30,,j]
  f.smooth = SpectrumEstimate(Y.sim, n_time, span=30, n_channel)
  cohY[,,j] = CoherenceEstimate(f.smooth,c(1:3))
}

result = matrix(0, 3, 2)
rownames(result) = c('Rank Signed Test', 'Rank sum Test',  'Rank Sum Test MVD')
colnames(result) = c('p-value', 'test statistics')
rsdt = rankSignedTest(cohX, cohY, n_channel, trials)
rst = rankSumTest(cohX, cohY, n_channel, trials)
rstmvd = rankSumTestMVD(cohX, cohY)
result[, 1] = c(rsdt$p.value, rst$p.value, rstmvd$p.value)
result[, 2] = c(rsdt$test.stat, rst$test.stat, rstmvd$statistic)
round(result, 2)

##### Check running time
start_time <- Sys.time()
rsdt = rankSignedTest(cohX, cohY, n_channel, trials)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
rst = rankSumTest(cohX, cohY, n_channel, trials)
end_time <- Sys.time()
end_time - start_time

