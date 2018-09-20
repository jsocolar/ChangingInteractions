library(gllvm)          # for fitting Latent Variable Models

User = "Jacob"

setwd(paste0("/Users/", User, "/Dropbox/Work/Sierra_nevada"))
load("bugs_data2.Rdata")
load("species.Rdata")
setwd(paste0("/Users/", User, "/Dropbox/Work/Sierra_nevada/Community_stability"))

nsp <- dim(bugs.data$X)[3]-bugs.data$nzeroes
nst <- dim(bugs.data$X)[1]

sum_rm <- function(x){
  return(sum(x, na.rm = T))
}
site_sums <- apply(bugs.data$X, c(1,3), sum_rm)

comm_sp <- species[which(colSums(site_sums[which(bugs.data$era == 0), ] > 1) > 0 & 
                           colSums(site_sums[which(bugs.data$era == 1), ] > 1) > 0)]

comm_sp

X2 <- bugs.data$X[,,which(species %in% comm_sp)]

nsp2 <- dim(X2)[3]


ssm2 <- as.data.frame(matrix(data=NA, nrow=nst, ncol=nsp2))
names(ssm2) <- species[which(species %in% comm_sp)]

for(i in 1:nst){
  for(j in 1:nsp2){
    ssm2[i,j] <- as.numeric(sum(X2[i,,j], na.rm=T)>0)
  }
}

ssmH <- ssm2[which(bugs.data$era==0),]
ssmM <- ssm2[which(bugs.data$era==1),]
ssmA <- rbind(ssmH, ssmM)
covars <- c(rep(0,64), rep(1,71))
Hmod <- gllvm(ssmH, family='binomial', num.lv=2)
Mmod <- gllvm(ssmM, family='binomial', num.lv=2)
Amod <- gllvm(as.matrix(ssmA), X=as.data.frame(covars), family='binomial', num.lv=2)

# The likelihood ratio test:
pchisq(2*(Hmod$logL + Mmod$logL - Amod$logL), 225, lower.tail = F)
# Note that the more complex model is Hmod + Mmod, so their log-likelihoods add.
# The less complex model is Amod.  The test statistic is -2(log(simpleL/complexL)),
# which is equal to 2(log(complexL/simpleL) = 2 (logL(complex) - logL(simple)).
# The difference in degrees of freedom is two times the number of species minus 1, i.e. 225.
# I confirmed this by back-calculating the number of degrees of freedom in each model
# based on the log-likelihoods and AIC scores provided by gllvm.

dimnames(Hmod$params$theta)

plot(Hmod$params$theta[,1], Hmod$params$theta[,2], pch="", xlab='LV1', ylab = 'LV2', main = 'Historical')
text(Hmod$params$theta[,1], Hmod$params$theta[,2], labels = comm_sp, cex=.5)

plot(-Mmod$params$theta[,1], -Mmod$params$theta[,2], pch="", xlab='LV1', ylab = 'LV2', main = 'Modern')
text(-Mmod$params$theta[,1], -Mmod$params$theta[,2], labels = comm_sp, cex=.5)

plot(Amod$params$theta[,1], Amod$params$theta[,2], pch="", xlab='LV1', ylab = 'LV2', main = 'Combined')
text(Amod$params$theta[,1], Amod$params$theta[,2], labels = comm_sp, cex=.5)

