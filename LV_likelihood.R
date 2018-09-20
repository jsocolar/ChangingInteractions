library(gllvm)          # for fitting Latent Variable Models

User = "Jacob"
setwd(paste0("/Users/", User, "/Dropbox/Work/Sierra_nevada"))
load("bugs_data2.Rdata")
load("species.Rdata")

setwd(paste0("/Users/", User, "/Dropbox/Work/Sierra_nevada/Community_stability"))

nsp <- dim(bugs.data$X)[3]-bugs.data$nzeroes
nst <- dim(bugs.data$X)[1]

# Get the total number of detections of each species at each site
sum_rm <- function(x){
  return(sum(x, na.rm = T))
}
site_sums <- apply(bugs.data$X, c(1,3), sum_rm)

# Get species detected at least once in both the historical and modern eras
both_sp <- species[which(colSums(site_sums[which(bugs.data$era == 0), ] > 1) > 0 & 
                           colSums(site_sums[which(bugs.data$era == 1), ] > 1) > 0)]

# Trim data to just those species
X2 <- bugs.data$X[,,which(species %in% both_sp)]

nsp2 <- dim(X2)[3]

# Create data structure
ssm2 <- as.data.frame(matrix(data=NA, nrow=nst, ncol=nsp2))
names(ssm2) <- species[which(species %in% both_sp)]
for(i in 1:nst){
  for(j in 1:nsp2){
    ssm2[i,j] <- as.numeric(sum(X2[i,,j], na.rm=T)>0)
  }
}

ssmH <- ssm2[which(bugs.data$era==0),] # Historical data only
ssmM <- ssm2[which(bugs.data$era==1),] # Modern data only
ssmA <- rbind(ssmH, ssmM) # All data, reordered to ensure that all historical data precedes 
                          # all modern (this should already be the case, but just to be sure...)
covars <- c(rep(0,64), rep(1,71)) # A covariate vector giving the era for the rows of ssmA

# Fitting the latent variable models using variational approximation. The fits are approximate
# and the gllvm manual recommends fitting each model several times and selecting the fit with
# the highest log-likelihood.
Hmod <- gllvm(ssmH, family='binomial', num.lv=2)
Mmod <- gllvm(ssmM, family='binomial', num.lv=2)
Amod <- gllvm(as.matrix(ssmA), X=as.data.frame(covars), family='binomial', num.lv=2)
Hlik <- as.vector(Hmod$logL)
Mlik <- as.vector(Mmod$logL)
Alik <- as.vector(Amod$logL)

for(i in 2:10){
  print(i)
  Hmod1 <- gllvm(ssmH, family='binomial', num.lv=2)
  Mmod1 <- gllvm(ssmM, family='binomial', num.lv=2)
  Amod1 <- gllvm(as.matrix(ssmA), X=as.data.frame(covars), family='binomial', num.lv=2)
  Hlik[i] <- Hmod1$logL
  Mlik[i] <- Mmod1$logL
  Alik[i] <- Amod1$logL
  if(Hmod1$logL > Hmod$logL){Hmod <- Hmod1}
  if(Mmod1$logL > Mmod$logL){Mmod <- Mmod1}
  if(Amod1$logL > Amod$logL){Amod <- Amod1}
}

max(Hlik) - min(Hlik)
max(Mlik) - min(Mlik)
max(Alik) - min(Alik)

# The likelihood ratio test:
pchisq(2*(Hmod$logL + Mmod$logL - Amod$logL), 225, lower.tail = F)
# Note that the more complex model is Hmod + Mmod, so their log-likelihoods add.
# The less complex model is Amod.  The test statistic is -2(log(simpleL/complexL)),
# which is equal to 2(log(complexL/simpleL) = 2 (logL(complex) - logL(simple)).
# The difference in degrees of freedom is two times the number of species minus 1, i.e. 225.
# I confirmed this by back-calculating the number of degrees of freedom in each model
# based on the log-likelihoods and AIC scores provided by gllvm.
# The p-value against the null-hypothesis is then 1 minus the lower-tail chi-squared probability, 
# which we can also get by supplying the $lower.tail = F$ argument. A significant test statistic
# rejects the null hypothesis that the simpler model is correct, and implies that not all
# restrictions on the simpler model hold.

# Plotting the resulting ordinations
plot(Hmod$params$theta[,1], Hmod$params$theta[,2], pch="", xlab='LV1', ylab = 'LV2', main = 'Historical')
text(Hmod$params$theta[,1], Hmod$params$theta[,2], labels = both_sp, cex=.5)

plot(Mmod$params$theta[,1], Mmod$params$theta[,2], pch="", xlab='LV1', ylab = 'LV2', main = 'Modern')
text(Mmod$params$theta[,1], Mmod$params$theta[,2], labels = both_sp, cex=.5)

plot(Amod$params$theta[,1], Amod$params$theta[,2], pch="", xlab='LV1', ylab = 'LV2', main = 'Combined')
text(Amod$params$theta[,1], Amod$params$theta[,2], labels = both_sp, cex=.5)


## Same as above, but only for common species.
# Common species are defined somewhat arbitrarily, using the cutoff of Tingley et al 2009,
# as species detected at at least 9 sites in both eras. This yielded 48 species for analysis.

comm_sp <- species[which(colSums(site_sums[which(bugs.data$era == 0), ] > 1) > 8 & 
                           colSums(site_sums[which(bugs.data$era == 1), ] > 1) > 8)]

# Trim data to just those species
X3 <- bugs.data$X[,,which(species %in% comm_sp)]

nsp3 <- dim(X3)[3]

# Create data structure
ssm3 <- as.data.frame(matrix(data=NA, nrow=nst, ncol=nsp3))
names(ssm3) <- species[which(species %in% comm_sp)]
for(i in 1:nst){
  for(j in 1:nsp3){
    ssm3[i,j] <- as.numeric(sum(X3[i,,j], na.rm=T)>0)
  }
}

ssmH.c <- ssm3[which(bugs.data$era==0),] # Historical data only
ssmM.c <- ssm3[which(bugs.data$era==1),] # Modern data only
ssmA.c <- rbind(ssmH.c, ssmM.c) # All data, reordered to ensure that all historical data precedes 
# all modern (this should already be the case, but just to be sure...)
covars <- c(rep(0,64), rep(1,71)) # A covariate vector giving the era for the rows of ssmA

# Fitting the latent variable models using variational approximation. The fits are approximate
# and the gllvm manual recommends fitting each model several times and selecting the fit with
# the highest log-likelihood.
Hmod.c <- gllvm(ssmH.c, family='binomial', num.lv=2)
Mmod.c <- gllvm(ssmM.c, family='binomial', num.lv=2)
Amod.c <- gllvm(as.matrix(ssmA.c), X=as.data.frame(covars), family='binomial', num.lv=2)
Hlik.c <- as.vector(Hmod.c$logL)
Mlik.c <- as.vector(Mmod.c$logL)
Alik.c <- as.vector(Amod.c$logL)

for(i in 2:10){
  print(i)
  Hmod1 <- gllvm(ssmH.c, family='binomial', num.lv=2)
  Mmod1 <- gllvm(ssmM.c, family='binomial', num.lv=2)
  Amod1 <- gllvm(as.matrix(ssmA.c), X=as.data.frame(covars), family='binomial', num.lv=2)
  Hlik.c[i] <- Hmod1$logL
  Mlik.c[i] <- Mmod1$logL
  Alik.c[i] <- Amod1$logL
  if(Hmod1$logL > Hmod.c$logL){Hmod.c <- Hmod1}
  if(Mmod1$logL > Mmod.c$logL){Mmod.c <- Mmod1}
  if(Amod1$logL > Amod.c$logL){Amod.c <- Amod1}
}

max(Hlik) - min(Hlik)
max(Mlik) - min(Mlik)
max(Alik) - min(Alik)

# The likelihood ratio test:
pchisq(2*(Hmod.c$logL + Mmod.c$logL - Amod.c$logL), 95, lower.tail = F)

# Plotting the resulting ordinations
plot(Hmod.c$params$theta[,1], Hmod.c$params$theta[,2], pch="", xlab='LV1', ylab = 'LV2', main = 'Historical')
text(Hmod.c$params$theta[,1], Hmod.c$params$theta[,2], labels = comm_sp, cex=.5)

plot(Mmod.c$params$theta[,1], Mmod.c$params$theta[,2], pch="", xlab='LV1', ylab = 'LV2', main = 'Modern')
text(Mmod.c$params$theta[,1], Mmod.c$params$theta[,2], labels = comm_sp, cex=.5)

plot(Amod.c$params$theta[,1], Amod.c$params$theta[,2], pch="", xlab='LV1', ylab = 'LV2', main = 'Combined')
text(Amod.c$params$theta[,1], Amod.c$params$theta[,2], labels = comm_sp, cex=.5)







