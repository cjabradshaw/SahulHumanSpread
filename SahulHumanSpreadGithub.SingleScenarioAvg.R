###################################################################
## Cellular automaton model of human spread through Sahul        ##
## SINGLE-SCENARIO AVERAGES (ITERATED)                           ##
## March 2020                                                    ##
## Corey JA Bradshaw, Flinders University                        ##
## corey.bradshaw@flinders.edu.au                                ##
## http://GlobalEcologyFlinders.com                              ##
## http://github.com/cjabradshaw                                 ##
###################################################################

## remove everything
rm(list = ls())

## libraries
library(sp)
library(rgdal)
library(raster)
#library(Scale)
library(oceanmap)
library(insol)
library(OceanView)
library(abind)
library(pracma)
library(binford)
library(rgl)
library(scatterplot3d) 
library(spatstat)
library(spatialEco)
library(SpatialPack)
library(foreach)
library(doSNOW)


## functions
# gradient function for immigration
# pI ~ a / (1 + (b * exp(-c * Rk)))
aI <- 0.95; bI <- 5000; cI <- 3
pI.func <- function(Rk) {
  pI <- aI / (1 + (bI * exp(-cI * Rk)))
  return(pI)}
pI.func(4)
xI <- seq(1,10,0.01)
yI <- pI.func(xI)
plot(xI,yI,type="l",xlab="Rk",ylab="pI")

# gradient function for emigration
aE <- 1; bE <- -3.2
pE.func <- function(Rk) {
  pE <- aE * exp(bE * Rk)
  return(pE)}
pE.func(0.5)
xE <- seq(0.01,1,0.01)
yE <- pE.func(xE)
plot(xE,yE,type="l",xlab="Rk",ylab="pE")
pE.out <- data.frame(xE,yE)
setwd("~/Documents/Papers/Palaeo/Sahul/Aus human spread/out/output data/")
write.table(pE.out,file="pE.csv",sep=",", row.names = F, col.names = T)

# stochastic beta sampler (single sample)
stoch.beta.func <- function(mu, var) {
  Sx <- rbeta(length(mu), (((1 - mu) / var - 1 / mu) * mu ^ 2), ((((1 - mu) / var - 1 / mu) * mu ^ 2)*(1 / mu - 1)))
  return(params=Sx)
}

# stochastic beta sampler (n samples)
stoch.n.beta.func <- function(n, mu, var) {
  Sx <- rbeta(n, (((1 - mu) / var - 1 / mu) * mu ^ 2), ((((1 - mu) / var - 1 / mu) * mu ^ 2)*(1 / mu - 1)))
  return(params=Sx)
}

# dynamics model function
Nproj.func <- function(Nt, rm, K) {
  Nt1 <- round(Nt * exp(rm*(1-(Nt/K))), 0)
  return(Nt1)
}

# rescale a range
rscale <- function (x, nx1, nx2, minx, maxx) {
  nx = nx1 + (nx2 - nx1) * (x - minx)/(maxx - minx)
  return(nx)
}

# matrix rotation
rot.mat <- function(x) t(apply(x, 2, rev))

# matrix poisson resampler
rpois.fun <- function(x,y,M) {
  rpois(1,M[x,y])
}
rpois.vec.fun <- Vectorize(rpois.fun,vectorize.args = c('x','y'))


####################################################
## set grids
####################################################

## relative density, carrying capacity & feedbacks
## NPP (LOVECLIM)
npp <- read.table("NppSahul(0-140ka).csv", header=T, sep=",") # 0.5 deg lat resolution
head(npp[,1:3])
not.na <- which(is.na(npp[,3:dim(npp)[2]]) == F, arr.ind=T)
upper.row <- as.numeric(not.na[1,1])
lower.row <- as.numeric(not.na[dim(not.na)[1],1])
min.lat <- max(npp[not.na[,1], 1])  
max.lat <- min(npp[not.na[,1], 1])
min.lon <- min(npp[not.na[,1], 2])
max.lon <- max(npp[not.na[,1], 2])

sahul.sub <- rep(0, dim(npp)[1])
for (n in 1:dim(npp)[1]) {
  sahul.sub[n] <- ifelse(npp[n,1] <= min.lat & npp[n,1] >= max.lat & npp[n,2] >= min.lon & npp[n,2] <= max.lon, 1, 0)
}  
sah.keep <- which(sahul.sub == 1)
npp.sah <- npp[sah.keep,]
head(npp.sah[,1:3])

## import sea level & palaeo-lakes layer (1 = land; 2 = palaeo-lake)
sll <- read.table("SeaLevelSahul(0-140ka)&LakeC.csv", header=T, sep=",") # 0.5 deg lat resolution
sll.sah <- sll[sah.keep,]
head(sll.sah)
range(sll.sah[,3:dim(sll.sah)[2]], na.rm=T)

## import ruggedness (average elevational slope) per cell
rug <- read.table("RuggednessSahul(0-140ka).csv", header=T, sep=",") # 0.5 deg lat resolution
sah.keep <- which(sahul.sub == 1)
rug.sah <- rug[sah.keep,]
head(rug.sah)
range(rug.sah[,3:dim(rug.sah)[2]], na.rm=T)

## set base human Leslie matrix from Bradshaw et al. (2014; PNAS)
source("~/Documents/Papers/Other/Global human population/ms/PNAS/R1/matrixOperators.r")

# Siler hazard h(x) (Gurven et al. 2007)
# average hunter-gatherer
a1 <- 0.422 # initial infant mortality rate (also known as αt)
b1 <- 1.131 # rate of mortality decline (also known as bt)
a2 <- 0.013 # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 1.47e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.086 # rate of mortality increase
longev <- 80
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model

l.inf <- exp(-a1/b1) # survival at infinite time
T.m <- 1/b1 # time constant at which maturity is approached
h.m <- a2 # hazard for mature animals
l.m <- exp(-a2*x) # survival
h.s <- a3*exp(b3*x) # hazard for senescence
l.s <- exp((a3/b3)*(1 - exp(b3*x))) # survival for senescence
f.x <- a3*exp(b3*x)*exp((a3/b3)/(1-exp(b3*x))) # probability density function
(log(a3) - log(a1)) / a3
T.s <- (1/b3) # modal survival time

## survival
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
Sx <- 1 - qx
sx <- lx[2:len.lx]/lx[1:(len.lx-1)]
mx <- 1 - sx
Lx <- (lx[1:(len.lx-1)] + lx[2:len.lx])/2
ex <- rev(cumsum(rev(Lx)))/lx[-len.lx]
ex.avg <- ex + x[-len.lx]

# set SD for Sx
Sx.sd <- 0.05 # can set to any value

# fertility (Walker et al. 2006)
primiparity.walker <- c(17.7,18.7,19.5,18.5,18.5,18.7,25.7,19,20.5,18.8,17.8,18.6,22.2,17,16.2,18.4)
prim.mean <- round(mean(primiparity.walker),0)
prim.lo <- round(quantile(primiparity.walker,probs=0.025),0)
prim.hi <- round(quantile(primiparity.walker,probs=0.975),0)
print(c(prim.lo, prim.mean, prim.hi))
dat.world13 <- read.table("world2013lifetable.csv", header=T, sep=",")
fert.world13 <- dat.world13$m.f
fert.trunc <- fert.world13[1:(longev+1)]
pfert.trunc <- fert.trunc/sum(fert.trunc)
fert.bentley <- 4.69/2 # Bentley 1985 for !Kung
fert.vec <- fert.bentley * pfert.trunc

## construct matrix
stages <- len.lx
popmat <- matrix(0,nrow=stages,ncol=stages)
colnames(popmat) <- x
rownames(popmat) <- x

## populate matrix
popmat[1,] <- fert.vec
diag(popmat[2:stages,]) <- Sx
popmat[stages,stages] <- 0 # Sx[stages-1]
popmat.orig <- popmat ## save original matrix

## matrix properties
max.lambda(popmat) ## 1-yr lambda
r.ann <- max.r(popmat) # rate of population change, 1-yr
r.ann
#stable.stage.dist(popmat) ## stable stage distribution
R.val(popmat,stages) # reproductive value
gen.l <- G.val(popmat,stages) # mean generation length

## for r.max (set Sx=1)
Sx.1 <- Sx
Sx.1[] <- 1
popmat.max <- popmat.orig
diag(popmat.max[2:stages,]) <- Sx.1
popmat.max[stages,stages] <- 0 # Sx[stages-1]
max.lambda(popmat.max) ## 1-yr lambda
rm.ann <- max.r(popmat.max) # rate of population change, 1-yr
rm.ann

#stable.stage.dist(popmat) ## stable stage distribution
R.val(popmat,stages) # reproductive value
gen.l <- G.val(popmat,stages) # mean generation length

### population dynamics parameters
# dynamical model
# Nt+1 = Nt * exp(rm*(1-(Nt/K))) - (E - I)
lambda.ann <- exp(r.ann) # annual lambda
r.max.NEE <- 2 * log(lambda.ann^gen.l) # human rmax at generational scale (arbitrarily double)
lambda.max.ann <- exp(rm.ann)
rm.max.NEE <- log(lambda.max.ann^gen.l) # human rmax at generational scale (from Sx=1 Leslie matrix)

# Cole's allometric calculation (high)
alpha.Ab <- 15
a.Cole <- -0.16
a.lo.Cole <- -0.41
a.up.Cole <- 0.10
a.sd.Cole <- mean(c((a.Cole - a.lo.Cole)/1.96, (a.up.Cole - a.Cole)/1.96))
b.lo.Cole <- -1.2
b.up.Cole <- -0.79
b.Cole <- -0.99
b.sd.Cole <- mean(c((b.Cole - b.lo.Cole)/1.96, (b.up.Cole - b.Cole)/1.96))
r.max.Cole <- 10^(a.Cole + b.Cole*log10(alpha.Ab)) # from Hone et al. 2003-JApplEcol
r.max.Cole
r.max.lo.Cole <- 10^(a.lo.Cole + b.lo.Cole*log10(alpha.Ab))
r.max.lo.Cole
r.max.up.Cole <- 10^(a.up.Cole + b.up.Cole*log10(alpha.Ab))
r.max.up.Cole

r.max.up.gen.Cole <- log((exp(r.max.up.Cole))^gen.l)
r.max.gen.Cole <- log((exp(r.max.Cole))^gen.l)
r.max.lo.gen.Cole <- log((exp(r.max.lo.Cole))^gen.l)
print(c(r.max.up.gen.Cole, r.max.gen.Cole, r.max.lo.gen.Cole))
r.max.gen.Cole.sd <- mean(c((r.max.gen.Cole - r.max.lo.gen.Cole)/1.96, (r.max.up.gen.Cole - r.max.gen.Cole)/1.96))

## dispersal calculations
# natal dispersal distance function (from Sutherland et al. 2000 Conserv Biol)
# median
a.mid <- 1.45
b.mid <- 0.54
a.lo <- a.mid - 1.05
a.up <- a.mid + 1.05
b.lo <- b.mid - 0.01
b.up <- b.mid + 0.01
M <- 50 # average mass of hunter-gatherer adult
D.med.lo <- a.lo*M^b.lo
D.med.up <- a.up*M^b.up
D.vec <- 1:round((10*111/2), 0)
Pr.Dmed.lo <- exp(-D.vec/D.med.lo) 
Pr.Dmed.up <- exp(-D.vec/D.med.up)

# max
A.mid <- 3.31
B.mid <- 0.65
A.lo <- A.mid - 1.17
A.up <- A.mid + 1.17
B.lo <- B.mid - 0.05
B.up <- B.mid + 0.05
D.max.lo <- A.lo*M^B.lo
D.max.up <- A.up*M^B.up
Pr.Dmx.lo <- exp(-D.vec/D.max.lo) 
Pr.Dmx.up <- exp(-D.vec/D.max.up)

par(mfrow=c(2,2))
plot(D.vec, Pr.Dmed.up, type="l", xlab="distance (km)", ylab="Pr(D)")
lines(D.vec, Pr.Dmed.lo)
lines(D.vec, Pr.Dmx.lo, col="red")
lines(D.vec, Pr.Dmx.up, col="red")

plot(log10(D.vec), Pr.Dmed.up, type="l", xlab="log10 distance (km)", ylab="Pr(D)")
lines(log10(D.vec), Pr.Dmed.lo)
lines(log10(D.vec), Pr.Dmx.lo, col="red")
lines(log10(D.vec), Pr.Dmx.up, col="red")

cellD <- D.vec/(111/2)
plot(cellD, Pr.Dmed.up, type="l", xlab="distance (0.5-deg cells)", ylab="Pr(D)")
lines(cellD, Pr.Dmed.lo)
lines(cellD, Pr.Dmx.lo, col="red")
lines(cellD, Pr.Dmx.up, col="red")

plot(log10(cellD), Pr.Dmed.up, type="l", xlab="log10 distance (0.5-deg cells)", ylab="Pr(D)")
lines(log10(cellD), Pr.Dmed.lo)
lines(log10(cellD), Pr.Dmx.lo, col="red")
lines(log10(cellD), Pr.Dmx.up, col="red")
par(mfrow=c(1,1))

disp.max.out <- data.frame(cellD,Pr.Dmx.lo,Pr.Dmx.up)

## Hiscock rainfall-territory size relationsip (2008)
terr.rain <- read.table("territory.rainfall.Hiscock.csv", header=T, sep=",")
plot(terr.rain$annrainmm, terr.rain$terrkm2, pch=19)
plot(log10(terr.rain$annrainmm), log10(terr.rain$terrkm2), pch=19)
fit.terr.rain <- lm(log10(terr.rain$terrkm2) ~ log10(terr.rain$annrainmm))
abline(fit.terr.rain, lty=2)
range(terr.rain$terrkm2)
median(terr.rain$terrkm2)

plot(log10(terr.rain$annrainmm), log10(terr.rain$terr.rkm), pch=19, ylim=c(-0.1,2), xlab="log10 annual rainfall (mm)", ylab="log10 territory radius (km)")
fit.dispkm.rain <- lm(log10(terr.rain$terr.rkm) ~ log10(terr.rain$annrainmm))
abline(fit.dispkm.rain, lty=2)
abline(h=log10(D.max.lo), lty=2, col="dark red")
abline(h=log10(D.max.up), lty=2, col="dark red")
abline(h=log10(D.med.lo), lty=3, col="red")
abline(h=log10(D.med.up), lty=3, col="red")

## make rainfall relative to minimum
hist(terr.rain$annrainmm/min(terr.rain$annrainmm), main="")
plot(log10(terr.rain$annrainmm/min(terr.rain$annrainmm)), log10(terr.rain$terr.rkm), pch=19, xlim=c(-1,1.6), ylim=c(-0.1,2), xlab="log10 relative annual rainfall", ylab="log10 territory radius (km)")
fit.dispkm.rain <- lm(log10(terr.rain$terr.rkm) ~ log10(terr.rain$annrainmm/min(terr.rain$annrainmm)))
abline(fit.dispkm.rain, lty=2)
abline(v=0, lty=2, col="black")
yDmaxup <- (log10(D.max.up) - as.numeric(coef(fit.dispkm.rain)[1]))/as.numeric(coef(fit.dispkm.rain)[2])
yDminlo <- (log10(D.med.lo) - as.numeric(coef(fit.dispkm.rain)[1]))/as.numeric(coef(fit.dispkm.rain)[2])
abline(v=yDmaxup, lty=2, col="black")
abline(h=log10(D.max.lo), lty=2, col="dark red")
abline(h=log10(D.max.up), lty=2, col="dark red")
abline(h=log10(D.med.lo), lty=3, col="red")
abline(h=log10(D.med.up), lty=3, col="red")

hiscock.out <- data.frame(terr.rain$annrainmm/min(terr.rain$annrainmm), log10(terr.rain$terr.rkm))

## Binford's environmental and hunter-gatherer frames of reference (Binford 2001)
## to estimate effect of rugosity on dispersal
binforddat = NULL
binforddat$annual_move <- LRB$dismov
binforddat$annual_rain <- LRB$bio.12
binforddat$ID <- seq(1:length(binforddat$annual_rain))
binford.dat <- as.data.frame(binforddat)
binford.dat$altitude_max <- LRB$h25
binford.dat$altitude_min <- LRB$l25
binford.dat$altitude_dif <- binford.dat$altitude_max - binford.dat$altitude_min
binford.dat <- binford.dat[complete.cases(binford.dat),]
binford.dat <- binford.dat [binford.dat$altitude_dif>0,]  # remove the one case with negative altitude difference
plot(binford.dat$altitude_dif, binford.dat$annual_move, pch=19, xlab="altitude difference", ylab="annual movement")

# cube root
binford.dat$annual_move_tr <- binford.dat$annual_move^(1/3)
binford.dat$annual_rain_tr <- binford.dat$annual_rain^(1/3)
binford.dat$altitude_dif_tr <- binford.dat$altitude_dif^(1/3)

res <- lm(annual_move_tr ~ + annual_rain + altitude_dif_tr, data = binford.dat)
summary(res)
#plot(res)

# 3D scatterplot with droplines
s3d <- scatterplot3d(binford.dat$annual_rain, binford.dat$altitude_dif_tr, binford.dat$annual_move_tr, pch=16, highlight.3d=F,
                     type="h", main="", angle=210, asp=4, grid=T, 
                     xlab="annual rainfall (mm)", ylab="cube root altitude difference (Δm)", zlab="cube root annual movement (km)")
s3d$plane3d(res, draw_polygon = T)
altdiff.vec <- seq(from=range(binford.dat$altitude_dif)[1], to=range(binford.dat$altitude_dif)[2], by=1)
altdiff.st <- altdiff.vec / max(altdiff.vec)
annmov.pred <- (coef(res)[1] + (altdiff.vec^(1/3) * coef(res)[3]^3) + (mean(binford.dat$annual_rain, na.rm=T) * coef(res)[2]))^3
annmov.st <- annmov.pred / max(annmov.pred)
plot(altdiff.vec, annmov.pred, type="l", xlab="altitude difference (m)", ylab="predicted annual movement (km)")
plot(altdiff.st, annmov.st, type="l", xlab="standardised altitude difference (m)", ylab="standardised predicted annual movement (km)")
altmov.dat <- data.frame(altdiff.st, annmov.st)

# exponential decay function fit
# y=a+b(x)^(1/3)
param.init <- c(1, -0.01)
fit.expd <- nls(annmov.st ~ a + b*(altdiff.st)^(1/3), 
                data = altmov.dat,
                algorithm = "port",
                start = c(a = param.init[1], b = param.init[2]),
                trace = TRUE,      
                nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
plot(altdiff.st, annmov.st, pch=19, cex=0.001, xlab="standardised altitude difference (m)", ylab="standardised predicted annual movement (km)", ylim=c(0.98,1))
altdiff.st.vec <- seq(0,1,0.01)
pred.annmov.st <- as.numeric(coef(fit.expd)[1]) + as.numeric(coef(fit.expd)[2]) * (altdiff.st.vec)^(1/3)
range(pred.annmov.st)
lines(altdiff.st.vec, pred.annmov.st, lty=2, lwd=3, col="red")

## distance to water (based on modern-day water distribution; Damien O'grady)
dir.tmp <- getwd()
rastlist <- list.files(path = dir.tmp, pattern='.tif$', all.files=TRUE, full.names=FALSE)
allrasters <- lapply(rastlist, raster)
allrasters
allrasters[[1]]
plot(allrasters[[1]], col=(topo.colors(50, 0.9)))

d2w <- read.table("Distance2Freshwater.csv", header=T, sep=",") # 0.5 deg lat resolution / units in dd
d2w.sah <- d2w[sah.keep,]
head(d2w.sah)
range(d2w.sah[,3:dim(d2w.sah)[2]], na.rm=T)

## data from Fagan & Holmes to estimate > mortality rate < MVP size
max.r.decl <- c(-3.24, -1.09, -1.96, -2.28, -0.69, -0.68, -1.88, -1.76, -2.05)
max.lam.decl <- exp(max.r.decl)
1-mean(max.lam.decl)



    ########################################################
    ########################################################
    ########################################################
    ## start diffusion model - iterate for average output ##
    ########################################################
    ########################################################
    ########################################################
    
    ### date of first colonisation?
    entry.date <- 50000
    
    ### how many generations to run?
    gen.run <- 300 # number of generations to run
    round(gen.run*gen.l, 0) # how many years is that?
    entry.date-round(gen.run*gen.l,0) # this means will run until year ...
    
    ### choose linear relationship between NPP and K, or 180 deg-rotated parabolic relationship
    #K.NPP <- "linear"
    #K.NPP <- "rotated parabolic"
    K.NPP <- "rotated quadratic yield density"
    
    ## K reduction scalar
    #modify.K <- "yes"
    modify.K <- "no"
    if (modify.K == "yes") {
      Kreduce <- 0.75 # if yes, by how much?
    } 
    
    ### for unknown SDs, choose % of mean (i.e., for M.cat.sd, pmov.sd, rm.max.NEE.sd)
    #SD.prop.xbar <- 0.05 # 5%
    SD.prop.xbar <- 0.10 # 10%
    
    # impose higher extinction probability below MVP size
    small.pop.ext <- "yes"
    #small.pop.ext <- "no"
    MVP.thresh <- 100 # increase mortality in cells with N < MVP.thresh
    ltMVP.red <- 0.2 # mean (beta resampled) additional mortality expected for cells meeting criteria
    
    ### choose low (2*NEE estimate) or high (generationally scaled Cole's estimate from alpha) rmax
    rmax.sel <- "low" # more defensible
    #rmax.sel <- "high" # seems unrealistically high
    
    ### if a lake is present, how much to reduce K in that grid (0 = 0 K; 1 = full K)?
    lake.red <- 0.01 # cannot be zero
    
    ### max long-distance dispersal modifier
    ldp.mult <- 1 # modify to deviate from theoretical expectations (0 - ∞)
    
    ## modifier of maximum ruggedness effect on movement
    rugmov.mod <- 1
    
    ## set resistance of landscape for distance-to-water function
    watmovresist <- 3 # from Saltré et al. 2009 Ecol Model
    
    ### apply catastrophe function (Reed et al. 2005)?
    catastrophe <- "yes"
    #catastrophe <- "no"
    pop.adjacency <- 0 # to account for 'population' area of influence, in incrementing neighbourhood (1 = immediate neighbours; 2 = 2 cells away in every direction, ...)
    cat.pr <- 0.14/((2*pop.adjacency+1)^2) # probability of catastrophe per generation (Reed et al. 2003) = 0.14, modified by adjacency from above
    M.cat <- 0.75 # mean mortality during a catastrophe event
    M.cat.sd <- SD.prop.xbar*M.cat # sd mortality during a catastrophe event
    
    ### are catastrophe's spatially clustered?
    spatial.cluster <- "yes"
    #spatial.cluster <- "no"
    
    # generate a random point pattern (Thomas cluster process)
    rpp.scale <- 0.015 # controls intensity of clustering (lower values = more clustering)
    kappa.mod <- 1 # modifies Thomas cluster process kappa upward or downward; ~ simulates changes to Pr(catastrophe)
    
    ### print each iteration's map (progression)?
    #print.map <- "yes"
    print.map <- "no"
    
    ### save each iteration's map as a jpg?
    save.map <- "no"
    #save.map <- "yes"
    
    # save output images?
    #save.outputs <- "yes"
    save.outputs <- "no"
    
    ### pick entry cell
    #start.col1 <- 45; start.row1 <- 1 # north of Bird's head, N Guinea
    #start.col1 <- 40; start.row1 <- 4 # west of Bird's head, N Guinea
    #start.col1 <- 38; start.row1 <- 21 # N Sahul shelf
    start.col1 <- 31; start.row1 <- 27 # S Sahul shelf
    #start.col1 <- 8; start.row1 <- 58 # SWA
    #start.col1 <- 87; start.row1 <- 59 # SNSW
    
    if (entry.date > 75000 & start.col1 == 31) {
      start.col1 <- 32}
    
    ### add secondary colonisation event?
    second.col <- "no"
    #second.col <- "yes"
    
    ### lag between first and second colonisation events
    lag.l <- 1000 # lag (between 1st & 2nd colonisation events) length?
    lag.n <- 2 # how many lag increments?
    start.time2 <- ifelse((lag.n * round(lag.l/gen.l)) == 0, 1, (lag.n * round(lag.l/gen.l))) # generations after first colonisation event (increments of ~ lag years)
    
    if (second.col == "yes") {
      #start.col2 <- 38; start.row2 <- 21 # N Sahul shelf
      #start.col2 <- 31; start.row2 <- 27 # S Sahul shelf
      start.col2 <- 40; start.row2 <- 4 # west of Bird's head, N Guinea
    }  
    
    if (entry.date > 75000 & start.col2 == 31) {
      start.col2 <- 32}
    
    # add this to jpg file names for scenario description
    name.add <- paste(entry.date/1000, "ka.", gen.run, "gen.", ifelse(K.NPP=="rotated quadratic yield density", "rqydK", ifelse(K.NPP=="linear", "linK", "parK")), ".rmax", ifelse(rmax.sel == "low","-lo","-hi"), ifelse(ldp.mult != 1, paste(".ldispmod", ldp.mult, sep=""), ""), ifelse(catastrophe=="yes",".CAT", ".NOCAT"), M.cat*100, ifelse(spatial.cluster=="yes",paste("cl",rpp.scale,sep=""),""), ".neigh-adj", pop.adjacency, ".", "intro", ifelse(second.col=="no",1,2), ifelse(start.row1==1, "nBH", ifelse(start.row1==4, "wBH", "SSh")), ifelse(second.col=="yes", ifelse(second.col=="yes" & start.row2==21,"-SSh", ifelse(start.row2==4, "-wBH", "-nBH")), ""), ifelse(second.col=="yes", paste("-lag", round(lag.l*lag.n/gen.l), "g", sep=""), ""), sep="") 
    
    ### minimum viable population size to consider a cell 'colonised' (for first-appearance map)
    min.pop.size <- 100
    
    ### start founding population in 1 cell
    N.found.mod <- 1
    N.found.lo <- 1300*N.found.mod; N.found.up <- 1500*N.found.mod
    
    ## estimate SD of rmax according to SD proportions et above
    rm.max.NEE.sd <- SD.prop.xbar * rm.max.NEE
    
    ### dispersal parameters
    pmov.mean <- 1/3 # mu for beta function indicating proportion of people immigrating/emigrating
    pmov.sd <- SD.prop.xbar*pmov.mean # sd for beta function
    NK.emig.min <- 0.3; NK.emig.max <- 0.7 # N/K when pmov.mean emigrates
    
    # update ruggedness movement function with rugmov.mod
    rugmovmod.func <- function(x) {
      rugmovmod <- as.numeric(coef(fit.expd)[1]) + rugmov.mod*as.numeric(coef(fit.expd)[2]) * (x)^(1/3)
      return(rugmovmod=rugmovmod)
    }
    
    # npp @ entry date ka
    sub.entry <- which(colnames(npp.sah) == paste("X",entry.date,sep=""))
    npp.sah.entry <- npp.sah[,c(1,2,sub.entry)]
    head(npp.sah.entry)
    
    # plot raster
    coordinates(npp.sah.entry) = ~ Lon + Lat
    proj4string(npp.sah.entry)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
    gridded(npp.sah.entry) = TRUE
    npp.entry = raster(npp.sah.entry)
    lim.exts <- 5
    plot(npp.entry, xlab="longitude", ylab="latitude", xlim=c(min.lon-lim.exts, max.lon+lim.exts), ylim=c(max.lat-lim.exts, min.lat+lim.exts))
    
    # transform to array
    lz <- dim(npp.sah)[2] - 2
    npp.array <- array(data=NA, dim=c(dim(raster2matrix(npp.entry)),lz))
    for (k in 3:(lz+2)) {
      npp.sah.k <- npp.sah[,c(1,2,k)] 
      coordinates(npp.sah.k) = ~ Lon + Lat
      proj4string(npp.sah.k)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
      gridded(npp.sah.k) = TRUE
      npp.k = raster(npp.sah.k)
      npp.array[,,k-2] <- raster2matrix(npp.k)
      image(npp.array[,,k-2], col=rev(grey(1:100/100)))
      print(k-2)
    }
    # example plots (entry & another)
    image(npp.array[,,sub.entry-3],col=rev(grey(1:100/100)))
    image(npp.array[,,12],col=rev(grey(1:100/100)))
    
    ## calculate all Ks as relative to current
    npp.sah.rel <- npp.array
    for (z in 1:lz) {
      npp.sah.rel[,,z] <- npp.array[,,z] / npp.array[,,3]
    }
    npp.sah.rel[,,3] <- npp.array[,,1]
    
    # npp to K
    hum.dens.med <- 6.022271e-02
    hum.dens.lq <- 3.213640e-02
    hum.dens.uq <- 1.439484e-01
    hum.dens.max <- 1.152206e+00
    hum.dens.min <- 1.751882e-02
    cell.area <- (111.12/2)*(111.12/2) # km2
    
    # create vector of K reduction scalars for projection interval
    if (modify.K == "yes") {
      Kreduce.vec <- stoch.n.beta.func(lz, Kreduce, 0.05*Kreduce)
    }
    if (modify.K == "no") {
      Kreduce <- rep(1,lz)
    } 
    
    # modify underlying K magnitude by modifying NPP across the board
    K.array <- npp.sah.rel
    for (z in 1:lz) {
      K.array[,,z] <- rscale(npp.array[,,z], round(hum.dens.min*cell.area, 0), round(hum.dens.max*cell.area, 0), min(npp.array[,,z], na.rm=T), max(npp.array[,,z], na.rm=T))
    }
    image(K.array[,,13], col=rev(grey(1:100/100)))
    
    # 180-deg rotated parabola
    # y = a(x - h)^2 + k
    # h = median NPP; k = max K; a = negative for 180 flipped
    k.Kmax <- max(K.array, na.rm=T)/2
    h.NPPmed <- mean(npp.array, na.rm=T)
    h.NPPmed <- mean(range(npp.array, na.rm=T))
    Kmin <- min(K.array, na.rm=T)
    NPP.seq <- seq(min(npp.array, na.rm=T), max(npp.array, na.rm=T), 0.01)
    K.parab.pred <- (-3 * (NPP.seq - h.NPPmed)^2) + k.Kmax
    K.parab.pred.rescale <- rscale(K.parab.pred, round(hum.dens.min*cell.area, 0), 0.5*round(hum.dens.max*cell.area, 0), min(K.parab.pred), max(K.parab.pred))
    plot(NPP.seq, K.parab.pred.rescale, type="l", xlab="NPP", ylab="K")
    
    # slow exponential increase combined with peak
    # reciprocal quadratic yield density
    # y = x/(a + b*x + c*x^2)
    # y = K, x = NPP
    a.rqyd <- 200; b.rqyd <- 0.60; c.rqyd <- 0.2
    K.rqyd.pred <- a.rqyd * exp(-(NPP.seq-b.rqyd)^2/(2*c.rqyd^2))
    K.rqyd.pred.rescale <- rscale(K.rqyd.pred, round(hum.dens.min*cell.area, 0), 0.5*round(hum.dens.max*cell.area, 0), min(K.rqyd.pred), max(K.rqyd.pred))
    plot(NPP.seq, K.rqyd.pred.rescale, type="l", xlab="NPP", ylab="K")
    
    K.lin.x <- c(min(npp.array, na.rm=T), max(npp.array, na.rm=T))
    K.lin.y <- c(min(K.array, na.rm=T), max(K.array, na.rm=T))
    plot(K.lin.x, K.lin.y, xlab="NPP", ylab="K", pch=19)
    fit.K.lin <- lm(K.lin.y ~ K.lin.x)
    abline(fit.K.lin, lty=2)
    K.lin.pred <- as.numeric(coef(fit.K.lin)[1]) + as.numeric(coef(fit.K.lin)[2])*NPP.seq
    plot(NPP.seq, K.lin.pred, type="l", lty=2, xlab="NPP", ylab="K", pch=19)
    lines(NPP.seq, K.parab.pred.rescale, lty=2)
    lines(NPP.seq, K.rqyd.pred.rescale, lty=2)
    abline(h=0.5*max(K.array, na.rm=T), lty=3)
    
    # rotated parabolic
    K.array.parab <- (-3 * (npp.array - h.NPPmed)^2) + k.Kmax
    K.array.parab.rescale <- K.array.parab
    for (z in 1:lz) {
      K.array.parab.rescale[,,z] <- rscale(K.array.parab[,,z], round(hum.dens.min*cell.area, 0), round(hum.dens.max*cell.area, 0), min(K.array.parab[,,z], na.rm=T), max(K.array.parab[,,z], na.rm=T))
    }
    image(K.array.parab.rescale[,,13], col=rev(grey(1:100/100)))
    
    # reciprocal quadratic yield density
    K.array.rqyd <- a.rqyd * exp(-(npp.array - b.rqyd)^2 / (2*c.rqyd^2))
    K.array.rqyd.rescale <- K.array.rqyd
    for (z in 1:lz) {
      K.array.rqyd.rescale[,,z] <- rscale(K.array.rqyd[,,z], round(hum.dens.min*cell.area, 0), round(hum.dens.max*cell.area, 0), min(K.array.rqyd[,,z], na.rm=T), max(K.array.rqyd[,,z], na.rm=T))
    }
    image(K.array.rqyd.rescale[,,13], col=rev(grey(1:100/100)))
    
    
    # rescale so that parabolic total K = linear total K
    hist.K.array <- hist(K.array, br=12)
    hist.K.array.dat <- data.frame(hist.K.array$mids, hist.K.array$density)
    sum(K.array, na.rm=T)
    hist(K.array.parab.rescale, br=12)
    sum(K.array.parab.rescale, na.rm=T)
    sum(K.array.parab.rescale, na.rm=T)/sum(K.array, na.rm=T)
    # rescale K.array.parab.rescale to same sum as K.array (distribution of Ks = same total)
    K.array.parab.rescale2 <- K.array.parab.rescale / (sum(K.array.parab.rescale, na.rm=T)/sum(K.array, na.rm=T))
    sum(K.array.parab.rescale2, na.rm=T)
    
    K.parab.pred.rescale2 <- K.parab.pred.rescale / (sum(K.array.parab.rescale, na.rm=T)/sum(K.array, na.rm=T))
    hist.K.parab.pred.rescale2 <- hist(K.parab.pred.rescale2,br=12)
    hist.K.parab.pred.rescale2.dat <- data.frame(hist.K.parab.pred.rescale2$mids, hist.K.parab.pred.rescale2$density)
    
    # rescale so that rotated quadratic yield density total K = linear total K
    hist(K.array.rqyd.rescale, br=12)
    sum(K.array.rqyd.rescale, na.rm=T)
    sum(K.array.rqyd.rescale, na.rm=T)/sum(K.array, na.rm=T)
    # rescale K.array.parab.rescale to same sum as K.array (distribution of Ks = same total)
    K.array.rqyd.rescale2 <- K.array.rqyd.rescale / (sum(K.array.rqyd.rescale, na.rm=T)/sum(K.array, na.rm=T))
    sum(K.array.rqyd.rescale2, na.rm=T)
    
    K.rqyd.pred.rescale2 <- K.rqyd.pred.rescale / (sum(K.array.rqyd.rescale, na.rm=T)/sum(K.array, na.rm=T))
    hist.K.rqyd.pred.rescale2 <- hist(K.rqyd.pred.rescale2,br=12)
    hist.K.rqyd.pred.rescale2.dat <- data.frame(hist.K.rqyd.pred.rescale2$mids, hist.K.rqyd.pred.rescale2$density)
    
    plot(NPP.seq, K.lin.pred, type="l", lty=2, xlab="NPP", ylab="K", pch=19)
    lines(NPP.seq, K.parab.pred.rescale2, lty=2)
    lines(NPP.seq, K.rqyd.pred.rescale2, lty=2)
    abline(h=0.5*max(K.array, na.rm=T), lty=3)
    NPP.K.out <- data.frame(NPP.seq,K.lin.pred,K.parab.pred.rescale2,K.rqyd.pred.rescale2)
    colnames(NPP.K.out) <- c("NPP","K.lin","K.para","K.rqyd")
  
    write.table(NPP.K.out,file="NPP.K.csv",sep=",", row.names = F, col.names = T)
    write.table(hist.K.rqyd.pred.rescale2.dat,file="hist.K.rqyd.pred.rescale2.csv",sep=",", row.names = F, col.names = T)
    
    # rotate matrix -90 & renumber from oldest to youngest
    if (K.NPP == "linear") {
      K.array.use <- K.array
    }
    if (K.NPP == "rotated parabolic") {
      K.array.use <- K.array.parab.rescale
    }
    if (K.NPP == "rotated quadratic yield density") {
      K.array.use <- K.array.rqyd.rescale
    }
    
    K.rot.array <- array(data=NA, c(dim(K.array.use)[2], dim(K.array.use)[1], lz))
    for (z in 1:lz) {
      K.rot.array[,,z] <- apply(t(K.array.use[,,142-z]),2,rev)
    }
    range(K.rot.array[,,1], na.rm=T)
    image(rot.mat(K.rot.array[,,81]), col=rev(grey(1:100/100)))
    
    if (modify.K == "yes") {
      for (z in 1:lz) {
        K.rot.array[,,z] <- K.rot.array[,,z] * Kreduce.vec[z]
      }
    }
    
    ## block out Indonesia & never-connected islands (make NA)
    K.rot.array[1:20, 1:39, ] <- NA
    K.rot.array[6:20, 40:44, ] <- NA
    K.rot.array[9:17, 45:46, ] <- NA
    K.rot.array[21:22, 24:33, ] <- NA
    K.rot.array[27:28, 24:26, ] <- NA
    K.rot.array[34:35, 77:83, ] <- NA
    K.rot.array[22:23, 85:87, ] <- NA
    K.rot.array[18, 84:85, ] <- NA
    K.rot.array[9:12, 77:86, ] <- NA
    K.rot.array[2:6, 70:87, ] <- NA
    K.rot.array[20, 71, ] <- NA
    K.rot.array[2, 52, ] <- NA
    
    # block passage to Tasmania (70 to 67K; 60 to 46K; 43-42K cannot cross)
    # K.rot.array[80, 68:70, c(71:74, 81:95, 98:99)] <- NA
    # K.rot.array[79, 74:77, c(71:74, 81:95, 98:99)] <- NA
    K.rot.array[80, 68:77, c(71:101)] <- NA
    K.rot.array[79, 68:77, c(71:101)] <- NA
    image(rot.mat(K.rot.array[,,96]), col=rev(grey(1:100/100)))
    
    # interpolate between 1000-year NPP intervals per human generation
    # interpolate Ks at gen.l intervals between 1000-yr slices
    mill.gen.div <- round(1000/gen.l, 0)
    subtentry <- dim(K.rot.array)[3] - (entry.date/1000)
    Kentry.array <- K.rot.array[,,subtentry:(dim(K.rot.array)[3])]
    image(rot.mat(Kentry.array[,,11]),col=rev(grey(1:100/100)))
    
    K.start <- entry.date
    K.end <- 0
    yr.proj.vec <- seq(K.start, K.end, -round((1000/mill.gen.div),0))
    Kentry.interp.array <- array(data=0, dim=c(dim(Kentry.array[,,1])[1],dim(Kentry.array[,,1])[2],length(yr.proj.vec)))
    mill.vec <- seq(K.start,K.end,-1000)
    
    for (i in 1:dim(Kentry.array)[1]) {
      for (j in 1:dim(Kentry.array)[2]) {
        if (length(which(is.na(Kentry.array[i,j,]))==T) <  (dim(Kentry.array)[3]-1))
          Kentry.interp.array[i,j,] <- approx(mill.vec, Kentry.array[i,j,], xout=yr.proj.vec, method="linear")$y
        else {
          Kentry.interp.array[i,j,] <- rep(NA, length(yr.proj.vec))
        }
      }
    }
    image(rot.mat(Kentry.interp.array[,,400]), col=rev(grey(1:100/100)))
    image(rot.mat(Kentry.interp.array[,,1]), col=rev(grey(1:100/100)))
    
    
    # transform sea level & palaeo-lakes file (sll) to an array
    lz <- dim(sll.sah)[2] - 2
    sll.array <- array(data=NA, dim=c(dim(raster2matrix(npp.entry)),lz))
    for (k in 3:(lz+2)) {
      sll.sah.k <- sll.sah[,c(1,2,k)] 
      coordinates(sll.sah.k) = ~ Lon + Lat
      proj4string(sll.sah.k)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
      gridded(sll.sah.k) = TRUE
      sll.k = raster(sll.sah.k)
      sll.array[,,k-2] <- raster2matrix(sll.k)
      image(sll.array[,,k-2], col=rev(grey(1:100/100)))
      print(k-2)
    }
    # example plots (entry & another)
    par(mfrow=c(1,2))
    image(sll.array[,,sub.entry],col=rev(grey(1:100/100)))
    image(sll.array[,,12],col=rev(grey(1:100/100)))
    par(mfrow=c(1,1))
    
    # rotate matrix -90 & renumber from oldest to youngest
    sll.rot.array <- array(data=NA, c(dim(sll.array)[2], dim(sll.array)[1], lz))
    for (z in 1:lz) {
      sll.rot.array[,,z] <- apply(t(sll.array[,,142-z]),2,rev)
    }
    range(sll.rot.array[,,1], na.rm=T)
    image(rot.mat(sll.rot.array[,,100]), col=(rainbow(2,1)))
    
    ## copy values between 1000-yr intervals
    sllentry.array <- sll.rot.array[,,subtentry:(dim(sll.rot.array)[3])]
    sllentry.copy.array <- array(data=0, dim=c(dim(sllentry.array[,,1])[1],dim(sllentry.array[,,1])[2],length(yr.proj.vec)))
    
    for (i in 1:dim(sllentry.array)[1]) {
      for (j in 1:dim(sllentry.array)[2]) {
        if (length(which(is.na(sllentry.array[i,j,]))==T) < (dim(sllentry.array)[3]-1))
          sllentry.copy.array[i,j,] <- approx(mill.vec, sllentry.array[i,j,], xout=yr.proj.vec, method="linear")$y
        else {
          sllentry.copy.array[i,j,] <- rep(NA, length(yr.proj.vec))
        }
      }
    }
    image(rot.mat(sllentry.copy.array[,,240]), col=(rainbow(2,1)))
    hist(sllentry.copy.array)
    
    # transform distance to water file (d2w) to an array
    lz <- dim(d2w.sah)[2] - 2
    d2w.array <- array(data=NA, dim=c(dim(raster2matrix(npp.entry)),lz))
    for (k in 3:(lz+2)) {
      d2w.sah.k <- d2w.sah[,c(1,2,k)] 
      coordinates(d2w.sah.k) = ~ Lon + Lat
      proj4string(d2w.sah.k)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
      gridded(d2w.sah.k) = TRUE
      d2w.k = raster(d2w.sah.k)
      d2w.array[,,k-2] <- raster2matrix(d2w.k)
      image(d2w.array[,,k-2], col=rev(grey(1:100/100)))
      print(k-2)
    }
    # just use matrix
    d2w1.mat <- (as.matrix(d2w.array[,,1]))

    # rotate matrix -90
    d2w.mat <- matrix(data=NA, nrow=dim(d2w1.mat)[2], ncol=dim(d2w1.mat)[1])
    d2w.mat <- apply(t(d2w1.mat),2,rev)
    range(d2w.mat, na.rm=T)
    image(rot.mat(d2w.mat), col=rev(grey(1:100/100)))

    # transform ruggedness file to array
    lz <- dim(rug.sah)[2] - 2
    rug.array <- array(data=NA, dim=c(dim(raster2matrix(npp.entry)),lz))
    for (k in 3:(lz+2)) {
      rug.sah.k <- rug.sah[,c(1,2,k)] 
      coordinates(rug.sah.k) = ~ Lon + Lat
      proj4string(rug.sah.k)=CRS("+proj=longlat +datum=WGS84") # set it to lat-long
      gridded(rug.sah.k) = TRUE
      rug.k = raster(rug.sah.k)
      rug.array[,,k-2] <- raster2matrix(rug.k)
      image(rug.array[,,k-2], col=rev(grey(1:100/100)))
      print(k-2)
    }
    # example plots (entry & another)
    par(mfrow=c(1,2))
    image(rug.array[,,sub.entry],col=rev(grey(1:100/100)))
    image(npp.array[,,sub.entry],col=rev(grey(1:100/100)))
    par(mfrow=c(1,1))
    par(mfrow=c(1,2))
    image(rug.array[,,12],col=rev(grey(1:100/100)))
    image(npp.array[,,12],col=rev(grey(1:100/100)))
    par(mfrow=c(1,1))
    
    # rescale rugosity from 0-1
    rug.array.rescale <- rug.array
    for (z in 1:lz) {
      rug.array.rescale[,,z] <- rscale(rug.array[,,z], 0, 1, min(rug.array[,,z], na.rm=T), max(rug.array[,,z], na.rm=T))
    }
    image(rug.array.rescale[,,13], col=rev(grey(1:100/100)))
    range(rug.array.rescale, na.rm=T)
    
    # rotate matrix -90 & renumber from oldest to youngest
    rug.rot.array <- array(data=NA, c(dim(rug.array.rescale)[2], dim(rug.array.rescale)[1], lz))
    for (z in 1:lz) {
      rug.rot.array[,,z] <- apply(t(rug.array.rescale[,,142-z]),2,rev)
    }
    range(rug.rot.array[,,1], na.rm=T)
    image(rot.mat(rug.rot.array[,,20]), col=rev(grey(1:100/100)))
    
    ## interpolate between 1000-yr intervals
    rugentry.array <- rug.rot.array[,,subtentry:(dim(rug.rot.array)[3])]
    rugentry.interp.array <- array(data=0, dim=c(dim(rugentry.array[,,1])[1],dim(rugentry.array[,,1])[2],length(yr.proj.vec)))
    
    for (i in 1:dim(rugentry.array)[1]) {
      for (j in 1:dim(rugentry.array)[2]) {
        if (length(which(is.na(rugentry.array[i,j,]))==T) <  (dim(rugentry.array)[3]-1))
          rugentry.interp.array[i,j,] <- approx(mill.vec, rugentry.array[i,j,], xout=yr.proj.vec, method="linear")$y
        else {
          rugentry.interp.array[i,j,] <- rep(NA, length(yr.proj.vec))
        }
      }
    }
    image(rot.mat(rugentry.interp.array[,,400]), col=rev(grey(1:100/100)))
    
    
    # spatial clustering of catastraphe events controlling parameters
    kappa.mult <- 0.9
    cellslo <- 1
    cellshi <- 3772
    kappa.mult.up <- 1.2*kappa.mod
    kappa.mult.lo <- 0.3*kappa.mod
    kappa.rep <- seq(kappa.mult.up,kappa.mult.lo,-(kappa.mult.up-kappa.mult.lo)/cellshi)
    cells.rep <- seq(cellslo,cellshi, (cellshi-cellslo)/cellshi)
    plot(cells.rep, kappa.rep, type="l")
    kappa.fit <- lm(kappa.rep ~ cells.rep)
    kappaP.func <- function(cells.occ) {
      kappa.pred <- (as.numeric(coef(kappa.fit)[1])) + as.numeric(coef(kappa.fit)[2])*cells.occ
      return(kappa.pred)
    }
    abline(h=kappaP.func(1500),lty=2,col="red")
    abline(v=(1500),lty=2,col="red")
    rpp.mu.mult <- 0.6 # this, with the fluctuating kappa multiplier, keeps overall mean proportion of cells experiencing catastrophic mortality ~ 0.14
    
    proc.sim.start <- proc.time()
    reps <- 100 # number of times to redo spread model to find average values
    first.arrive.array <- array(data=NA, dim=c(dim(rugentry.interp.array[,,1]), reps))
    sub97.vec <- rep(0,reps)
    N.mat <- matrix(data=NA,nrow=reps,ncol=(gen.run+1))
    
    ##########################
    ## start of m reps loop ##
    ##########################
    
    for (m in 1:reps) {
      
      ## parallel set-up
      cores <- 12
      getDoParWorkers()
      getDoParName()
      cl <- makeCluster(cores, type = "SOCK")
      registerDoSNOW(cl)
      getDoParWorkers()
      getDoParName()
      getDoParVersion()
      
      # set up NA array
      NA.array <- Kentry.interp.array * 0
      NA.array <- NA.array[,,1:(gen.run+1)]
      image(rot.mat(K.rot.array[,,subtentry]), col=rev(grey(1:100/100)))
      
      # set direction codes
      dir.vec <- c("NW", "N", "NE", "W", "E", "SW", "S", "SE")
      
      # generations to run
      gen.run * gen.l # numbers years to run
      
      # set up N array
      array.N <- Kentry.interp.array * 0
      array.N <- array.N[,,1:(gen.run+1)]
      
      # set up direction array (dominant direction of influx)
      dir.array <- array.N
      
      array.N[start.row1, start.col1, 1] <- round(runif(1, N.found.lo, N.found.up), 0)
      image((rot.mat(array.N[,,1])), col=rev(grey(1:100/100)))
      
      # log10 relative K array
      Kentry.interp.lrel.array <- log10(Kentry.interp.array / min(Kentry.interp.array, na.rm=T)) # log10 relative K
      
      ## assume same slope between relative NPP and radius movement
      disp.npp.slope <- as.numeric(coefficients((fit.dispkm.rain))[2])
      disp.npp.int <- as.numeric(coefficients((fit.dispkm.rain))[1])
      disp.npp.max.int <- 1.6646371
      
      i.rows <- dim(array.N[,,1])[1]
      j.cols <- dim(array.N[,,1])[2]
      z.layers <- dim(array.N)[3]
      
      # storage vectors
      N.vec <- poparea.vec <- pc.complete <- cat.pr.est <- rep(0,z.layers)
      N.vec[1] <- array.N[start.row1,start.col1,1]
      if (second.col == "yes") {
        N.vec[start.time2] <- array.N[start.row1,start.col1,start.time2]
      }
      poparea.vec[1] <- cell.area/1000
      proc.start <- proc.time() 
      
      #############################
      ## project
      for (t in 1:(z.layers-1)) {
        
        # Poisson-resampled K matrix
        Kentry.interp.poiss <- outer(1:nrow(round(Kentry.interp.array[,,t],0)), 1:ncol(round(Kentry.interp.array[,,t],0)), rpois.vec.fun, round(Kentry.interp.array[,,t],0))
        
        # reduce Ks if lake present
        Kentry.interp.pois <- Kentry.interp.poiss
        for (i in 1:i.rows) { # i rows
          for (j in 1:j.cols) { # j columns
            Kentry.interp.pois[i,j] <- ifelse(sllentry.copy.array[i,j,t] > 1, lake.red * Kentry.interp.poiss[i,j], Kentry.interp.poiss[i,j])
          }
        }
        
        # step through array in t for immigration & emigration
        for (i in 1:i.rows) { # i rows
          for (j in 1:j.cols) { # j columns
            
            if (second.col == "yes") {
              if (t == start.time2) {
                array.N[start.row2, start.col2, start.time2] <- round(runif(1, N.found.lo, N.found.up), 0)
              }
            }
            
            # set cell-cell gradients relative to focal cell
            NW.RK <- ifelse(i > 1 & j > 1, (Kentry.interp.pois[i,j] / ifelse(length(Kentry.interp.pois[i-1,j-1]) > 0, Kentry.interp.pois[i-1,j-1], NA)), NA)  # NW
            N.RK <- ifelse(i > 1, (Kentry.interp.pois[i,j] / ifelse(length(Kentry.interp.pois[i-1,j]) > 0, Kentry.interp.pois[i-1,j], NA)), NA) # N
            NE.RK <- ifelse(i > 1 & j < j.cols, (Kentry.interp.pois[i,j] / ifelse(length(Kentry.interp.pois[i-1,j+1]) > 0, Kentry.interp.pois[i-1,j+1], NA)), NA)  # NE
            W.RK <- ifelse(j > 1, (Kentry.interp.pois[i,j] / ifelse(length(Kentry.interp.pois[i,j-1]) > 0, Kentry.interp.pois[i,j-1], NA)), NA)    # W
            E.RK <- ifelse(j < j.cols, (Kentry.interp.pois[i,j] / ifelse(length(Kentry.interp.pois[i,j+1]) > 0, Kentry.interp.pois[i,j+1], NA)), NA)      # E
            SW.RK <- ifelse(i < i.rows & j > 1, (Kentry.interp.pois[i,j] / ifelse(length(Kentry.interp.pois[i+1,j-1]) > 0, Kentry.interp.pois[i+1,j-1], NA)), NA)  # SW
            S.RK <- ifelse(i < i.rows, (Kentry.interp.pois[i,j] / ifelse(length(Kentry.interp.pois[i+1,j]) > 0, Kentry.interp.pois[i+1,j], NA)), NA)      # S
            SE.RK <- ifelse(i < i.rows & j < j.cols, (Kentry.interp.pois[i,j] / ifelse(length(Kentry.interp.pois[i+1,j+1]) > 0, Kentry.interp.pois[i+1,j+1], NA)), NA)  # SE
            
            ## current population distances from K in each cell
            focal.dK1 <- array.N[i,j,t] / Kentry.interp.pois[i,j]
            NW.dK1 <- ifelse(i > 1 & j > 1, array.N[i-1,j-1,t] / Kentry.interp.pois[i-1,j-1], NA)
            N.dK1 <- ifelse(i > 1, array.N[i-1,j,t] / Kentry.interp.pois[i-1,j], NA)
            NE.dK1 <- ifelse(i > 1 & j < j.cols, array.N[i-1,j+1,t] / Kentry.interp.pois[i-1,j+1], NA)
            W.dK1 <- ifelse(j > 1, array.N[i,j-1,t] / Kentry.interp.pois[i,j-1], NA)
            E.dK1 <- ifelse(j < j.cols, array.N[i,j+1,t] / Kentry.interp.pois[i,j+1], NA)
            SW.dK1 <- ifelse(i < i.rows & j > 1, array.N[i+1,j-1,t] / Kentry.interp.pois[i+1,j-1], NA)
            S.dK1 <- ifelse(i < i.rows, array.N[i+1,j,t] / Kentry.interp.pois[i+1,j], NA)
            SE.dK1 <- ifelse(i < i.rows & j < j.cols, array.N[i+1,j+1,t] / Kentry.interp.pois[i+1,j+1], NA)
            
            focal.dK <- ifelse(focal.dK1 == 0 | is.na(focal.dK1) == T, 0, focal.dK1)
            NW.dK <- ifelse(length(NW.dK1) == 0 | is.na(NW.dK1) == T, 0, NW.dK1)
            N.dK <- ifelse(length(N.dK1) == 0 | is.na(N.dK1) == T, 0, N.dK1)
            NE.dK <- ifelse(length(NE.dK1) == 0 | is.na(NE.dK1) == T, 0, NE.dK1)
            W.dK <- ifelse(length(W.dK1) == 0 | is.na(W.dK1) == T, 0, W.dK1)
            E.dK <- ifelse(length(E.dK1) == 0 | is.na(E.dK1) == T, 0, E.dK1)
            SW.dK <- ifelse(length(SW.dK1) == 0 | is.na(SW.dK1) == T, 0, SW.dK1)
            S.dK <- ifelse(length(S.dK1) == 0 | is.na(S.dK1) == T, 0, S.dK1)
            SE.dK <- ifelse(length(SE.dK1) == 0 | is.na(SE.dK1) == T, 0, SE.dK1)
            
            # direction indices initialised as NA
            NWdir <- Ndir <- NEdir <- Wdir <- Edir <- SWdir <- Sdir <- SEdir <- NA
            
            # immigration into focal cell
            if (is.na(NW.RK) == F & NW.RK > 1 & NW.dK >= runif(1, min=NK.emig.min, max=NK.emig.max)) {
              NW.I <- (pI.func(NW.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i-1,j-1,t] * (rugmovmod.func(rugentry.interp.array[i-1,j-1,t])))
              array.N[i,j,t] <- sum(c(array.N[i,j,t], NW.I), na.rm=T)
              array.N[i-1,j-1,t] <- sum(c(array.N[i-1,j-1,t], -NW.I), na.rm=T)
              NWdir <- ifelse(is.na(NW.RK) == F & NW.RK > 1, NW.I, NA)}
            if (is.na(N.RK) == F & N.RK > 1 & N.dK >= runif(1, min=NK.emig.min, max=NK.emig.max)) {
              N.I <- (pI.func(N.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i-1,j,t] * (rugmovmod.func(rugentry.interp.array[i-1,j,t])))
              array.N[i,j,t] <- sum(c(array.N[i,j,t], N.I), na.rm=T)
              array.N[i-1,j,t] <- sum(c(array.N[i-1,j,t], -N.I), na.rm=T)
              Ndir <- ifelse(is.na(N.RK) == F & N.RK > 1, N.I, NA)}
            if (is.na(NE.RK) == F & NE.RK > 1 & NE.dK >= runif(1, min=NK.emig.min, max=NK.emig.max)) {
              NE.I <- (pI.func(NE.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i-1,j+1,t] * (rugmovmod.func(rugentry.interp.array[i-1,j+1,t])))
              array.N[i,j,t] <- sum(c(array.N[i,j,t], NE.I), na.rm=T)
              array.N[i-1,j+1,t] <- sum(c(array.N[i-1,j+1,t], -NE.I), na.rm=T)
              NEdir <- ifelse(is.na(NE.RK) == F & NE.RK > 1, NE.I, NA)}
            if (is.na(W.RK) == F & W.RK > 1 & W.dK >= runif(1, min=NK.emig.min, max=NK.emig.max)) {
              W.I <- (pI.func(W.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i,j-1,t] * (rugmovmod.func(rugentry.interp.array[i,j-1,t])))
              array.N[i,j,t] <- sum(c(array.N[i,j,t], W.I), na.rm=T)
              array.N[i,j-1,t] <- sum(c(array.N[i,j-1,t], -W.I), na.rm=T)
              Wdir <- ifelse(is.na(W.RK) == F & W.RK > 1, W.I, NA)}
            if (is.na(E.RK) == F & E.RK > 1 & E.dK >= runif(1, min=NK.emig.min, max=NK.emig.max)) {
              E.I <- (pI.func(E.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i,j+1,t] * (rugmovmod.func(rugentry.interp.array[i,j+1,t])))
              array.N[i,j,t] <- sum(c(array.N[i,j,t], E.I), na.rm=T)
              array.N[i,j+1,t] <- sum(c(array.N[i,j+1,t], -E.I), na.rm=T)
              Edir <- ifelse(is.na(E.RK) == F & E.RK > 1, E.I, NA)}
            if (is.na(SW.RK) == F & SW.RK > 1 & SW.dK >= runif(1, min=NK.emig.min, max=NK.emig.max)) {
              SW.I <- (pI.func(SW.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i+1,j-1,t] * (rugmovmod.func(rugentry.interp.array[i+1,j-1,t])))
              array.N[i,j,t] <- sum(c(array.N[i,j,t], SW.I), na.rm=T)
              array.N[i+1,j-1,t] <- sum(c(array.N[i+1,j-1,t], -SW.I), na.rm=T)
              SWdir <- ifelse(is.na(SW.RK) == F & SW.RK > 1, SW.I, NA)}
            if (is.na(S.RK) == F & S.RK > 1 & S.dK >= runif(1, min=NK.emig.min, max=NK.emig.max)) {
              S.I <- (pI.func(S.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i+1,j,t] * (rugmovmod.func(rugentry.interp.array[i+1,j,t])))
              array.N[i,j,t] <- sum(c(array.N[i,j,t], S.I), na.rm=T)
              array.N[i+1,j,t] <- sum(c(array.N[i+1,j,t], -S.I), na.rm=T)
              Sdir <- ifelse(is.na(S.RK) == F & S.RK > 1, S.I, NA)}
            if (is.na(SE.RK) == F & SE.RK > 1 & SE.dK >= runif(1, min=NK.emig.min, max=NK.emig.max)) {
              SE.I <- (pI.func(SE.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i+1,j+1,t] * (rugmovmod.func(rugentry.interp.array[i+1,j+1,t])))
              array.N[i,j,t] <- sum(c(array.N[i,j,t], SE.I), na.rm=T)
              array.N[i+1,j+1,t] <- sum(c(array.N[i+1,j+1,t], -SE.I), na.rm=T)
              SEdir <- ifelse(is.na(SE.RK) == F & SE.RK > 1, SE.I, NA)}
            
            # direction of dominant influx per time layer
            I.vec <- c(NWdir, Ndir, NEdir, Wdir, Edir, SWdir, Sdir, SEdir)
            dir.array[i,j,t] <- ifelse((length(which(I.vec > 1))) > 0, (dir.vec[max(which(I.vec > 1))]), NA)
            
            # emigration out of focal cell
            if ((ifelse(focal.dK >= runif(1, min=NK.emig.min, max=NK.emig.max), 1, 0)) == 1) {
              if (is.na(NW.RK) == F & NW.RK <= 1) {
                NW.E <- (pE.func(NW.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i,j,t] * (rugmovmod.func(rugentry.interp.array[i,j,t])))
                array.N[i,j,t] <- sum(c(array.N[i,j,t], ifelse(NW.E < 0, 0, -NW.E)), na.rm=T)
                array.N[i-1,j-1,t] <- sum(c(array.N[i-1,j-1,t], ifelse(NW.E < 0, 0, NW.E)), na.rm=T)}
              if (is.na(N.RK) == F & N.RK <= 1) {
                N.E <- (pE.func(N.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i,j,t] * (rugmovmod.func(rugentry.interp.array[i,j,t]))) - NW.E
                array.N[i,j,t] <- sum(c(array.N[i,j,t], ifelse(N.E < 0, 0, -N.E)), na.rm=T)
                array.N[i-1,j,t] <- sum(c(array.N[i-1,j,t], ifelse(N.E < 0, 0, N.E)), na.rm=T)}
              if (is.na(NE.RK) == F & NE.RK <= 1) {
                NE.E <- (pE.func(NE.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i,j,t] * (rugmovmod.func(rugentry.interp.array[i,j,t]))) - NW.E - N.E
                array.N[i,j,t] <- sum(c(array.N[i,j,t], ifelse(NE.E < 0, 0, -NE.E)), na.rm=T)
                array.N[i-1,j+1,t] <- sum(c(array.N[i-1,j+1,t], ifelse(NE.E < 0, 0, NE.E)), na.rm=T)}
              if (is.na(W.RK) == F & W.RK <= 1) {
                W.E <- (pE.func(W.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i,j,t] * (rugmovmod.func(rugentry.interp.array[i,j,t]))) - NW.E - N.E - NE.E
                array.N[i,j,t] <- sum(c(array.N[i,j,t], ifelse(W.E < 0, 0, -W.E)), na.rm=T)
                array.N[i,j-1,t] <- sum(c(array.N[i,j-1,t], ifelse(W.E < 0, 0, W.E)), na.rm=T)}
              if (is.na(E.RK) == F & E.RK <= 1) {
                E.E <- (pE.func(E.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i,j,t] * (rugmovmod.func(rugentry.interp.array[i,j,t]))) - NW.E - N.E - NE.E - W.E
                array.N[i,j,t] <- sum(c(array.N[i,j,t], ifelse(E.E < 0, 0, -E.E)), na.rm=T)
                array.N[i,j+1,t] <- sum(c(array.N[i,j+1,t], ifelse(E.E < 0, 0, E.E)), na.rm=T)}
              if (is.na(SW.RK) == F & SW.RK <= 1) {
                SW.E <- (pE.func(SW.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i,j,t] * (rugmovmod.func(rugentry.interp.array[i,j,t]))) - NW.E - N.E - NE.E - W.E - E.E
                array.N[i,j,t] <- sum(c(array.N[i,j,t], ifelse(SW.E < 0, 0, -SW.E)), na.rm=T)
                array.N[i+1,j-1,t] <- sum(c(array.N[i+1,j-1,t], ifelse(SW.E < 0, 0, SW.E)), na.rm=T)}
              if (is.na(S.RK) == F & S.RK <= 1) {
                S.E <- (pE.func(S.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i,j,t] * (rugmovmod.func(rugentry.interp.array[i,j,t]))) - NW.E - N.E - NE.E - W.E - E.E - SW.E
                array.N[i,j,t] <- sum(c(array.N[i,j,t], ifelse(S.E < 0, 0, -S.E)), na.rm=T)
                array.N[i+1,j,t] <- sum(c(array.N[i+1,j,t], ifelse(S.E < 0, 0, S.E)), na.rm=T)}
              if (is.na(SE.RK) == F & SE.RK <= 1) {
                SE.E <- (pE.func(SE.RK) * stoch.beta.func(pmov.mean, pmov.sd) * array.N[i,j,t] * (rugmovmod.func(rugentry.interp.array[i,j,t]))) - NW.E - N.E - NE.E - W.E - E.E - SW.E - S.E
                array.N[i,j,t] <- sum(c(array.N[i,j,t], ifelse(SE.E < 0, 0, -SE.E)), na.rm=T)
                array.N[i+1,j+1,t] <- sum(c(array.N[i+1,j+1,t], ifelse(SE.E < 0, 0, SE.E)), na.rm=T)}
            }
            
            # reset emigration values to zero for next run
            NW.E <- N.E <- NE.E <- W.E <- E.E <- SW.E <- S.E <- SE.E <- 0
            
            # long-range dispersal
            disp.prob <- (exp(-D.vec/(ifelse(Kentry.interp.lrel.array[i, j, t] >= log10(6.706985), D.max.lo, 10^(disp.npp.max.int + as.numeric(disp.npp.slope) * Kentry.interp.lrel.array[i, j, t])))))
            if (is.na(disp.prob[1])==F) {
              disp.prob.ran <- disp.prob} 
            if (is.na(disp.prob[1])==T) {
              disp.prob.ran <- ((Pr.Dmx.up+Pr.Dmx.lo)/2)}
            
            cellD.move <- ldp.mult * ((sample(cellD, size=1, replace=F, prob=disp.prob.ran)) * (2*focal.dK)) # multiply probability upwards if closer to focal cell K
            
            # condition on distance 2 water, where
            # P(reach) = 1 – ((D/max(D))^J); D = distance to water; higher J = more difficult reach gridcell; max(D) = max distance to water
            if (is.na(d2w.mat[i,j]) == F & cellD.move < d2w.mat[i,j]) {
              P.reach <- 1 - ((cellD.move / d2w.mat[i,j])^(watmovresist))
              reach.trial <- rbinom(1, 1, prob=P.reach)
            }
            reach.succ <- ifelse(is.na(d2w.mat[i,j]) == F, reach.trial, 1)
            
            if ((round(cellD.move)) > 0 & reach.succ == 1) {
              dx <- sample(c(-1,1), 1) * rpois(1,(round(cellD.move)))
              dy <- sample(c(-1,1), 1) * rpois(1,(round(cellD.move)))
              if ((i+dy) > 0 & (i+dy) <= i.rows & (j+dx) > 0 & (j+dx) <= j.cols) {
                if (length(which(is.na(array.N[(i+sign(dy)):(i+dy), (j+sign(dx)):(j+dx), t]) == T)) == 0) { # if there is an NA cell anywhere in block between [i,j] & [i+dy,j+dx], don't let dispersal happen
                  N.disp <- round(array.N[i, j, t] * stoch.beta.func(pmov.mean/10, pmov.sd/10) * (rugmovmod.func(rugentry.interp.array[i,j,t])), 0) # number dispersing to new cell
                  array.N[(i+dy), (j+dx), t] <- array.N[(i+dy), (j+dx), t] + N.disp
                  array.N[i, j, t] <- ifelse((array.N[i, j, t] - N.disp) < 0, 0, (array.N[i, j, t] - N.disp))
                } # end if
              } # end if
            } # end if
            
          } # end i loop
        } # end j loop
        
        # remove negative values
        array.N[,,t] <- ifelse(array.N[,,t] < 0, 0, round(array.N[,,t], 0))
        
        # apply dynamical model after movements from previous step
        if (rmax.sel == "low") {
          array.N[,,t+1] <- Nproj.func(Nt=array.N[,,t], rm=rnorm(1, rm.max.NEE, rm.max.NEE.sd), K=Kentry.interp.pois)
        }
        if (rmax.sel == "high") {
          array.N[,,t+1] <- Nproj.func(Nt=array.N[,,t], rm=log((exp(10^(rnorm(1, a.Cole, a.sd.Cole) + rnorm(1, b.Cole, b.sd.Cole) * log10(alpha.Ab))))^gen.l), K=Kentry.interp.pois)
        }
        
        # apply catastrophe; resampling cat.pr (binomial) & M.cat mortality (beta)
        if (catastrophe == "yes") {
          mat.noNA.g0 <- which(is.na(array.N[,,t+1]) != T & array.N[,,t+1] > 0, arr.ind=F)
          #if (spatial.cluster == "no") {
          mat.noNA.g0.sel <- rbinom(length(mat.noNA.g0), 1, cat.pr) # spatially random
          
          if (spatial.cluster == "yes" & length(mat.noNA.g0) > 1) {
            mat.noNA.g0.arr <- which(is.na(array.N[,,t+1]) != T & array.N[,,t+1] > 0, arr.ind=T)
            rows.ng0 <- sort(unique(mat.noNA.g0.arr[,1]))
            cols.ng0 <- sort(unique(mat.noNA.g0.arr[,2]))
            lrowcol.mn <- round(mean(length(rows.ng0),length(cols.ng0)), 0)
            mat.N.it <- array.N[,,t+1]
            Y <- rThomas(round(kappaP.func(length(which(array.N[,,t+1] > 0, arr.ind=F)))*lrowcol.mn, 0), rpp.scale, round(rpp.mu.mult*lrowcol.mn,0), drop=T, nsim = 1)
            #Y <- rThomas(round(kappa.mult*lrowcol.mn, 0), rpp.scale, round(rpp.mu.mult*lrowcol.mn,0), drop=T, nsim = 1)
            Yrows.calc <- round(Y$y * length(rows.ng0), 0)
            Yrows <- Yrows.calc + min(rows.ng0) - 1
            Yrows <- ifelse((Yrows.calc + min(rows.ng0) - 1) == 0, 1, (Yrows.calc + min(rows.ng0) - 1)) 
            Ycols.calc <- round(Y$x * length(cols.ng0), 0)
            Ycols <- Yrows.calc + min(cols.ng0) - 1
            Ycols <- ifelse((Ycols.calc + min(cols.ng0) - 1) == 0, 1, (Ycols.calc + min(cols.ng0) - 1)) 
            Yrowscols <- unique(cbind(Yrows,Ycols))
            mat.N.it[Yrowscols] <- round((stoch.n.beta.func(length(mat.N.it[Yrowscols]), M.cat, M.cat.sd)) * (mat.N.it[Yrowscols]), 0)
            mat.noNA.g0.spat <- which(mat.N.it %in% na.omit(mat.N.it[Yrowscols])[na.omit(mat.N.it[Yrowscols]) > 0])
            cat.pr.est[t] <- length(mat.noNA.g0.spat)/length(mat.noNA.g0)
            #cat.pr.mx <- stoch.beta.func(cat.pr, cat.pr*SD.prop.xbar)
            cat.pr.mx <- rnorm(1, cat.pr, cat.pr*SD.prop.xbar)
            if (cat.pr.est[t] > cat.pr.mx & length(mat.noNA.g0.spat) > 0 & length(mat.noNA.g0) > 1 & spatial.cluster=="yes") {
              mat.noNA.g0.spat <- sample(mat.noNA.g0.spat, round(cat.pr.mx * length(mat.noNA.g0)), replace=F)
              cat.pr.est[t] <- length(mat.noNA.g0.spat)/dim(mat.noNA.g0.arr)[1]}
            if (length(mat.noNA.g0.spat) > 0 & length(mat.noNA.g0) > 1 & spatial.cluster=="yes") {
              array.N[,,t+1][mat.noNA.g0.spat] <- round((stoch.n.beta.func(length(mat.noNA.g0.spat), M.cat, M.cat.sd)) * ((array.N[,,t+1])[mat.noNA.g0.spat]), 0)}
            if (length(mat.noNA.g0.spat) == 0 & length(mat.noNA.g0) > 1 & spatial.cluster=="yes") {
              mat.noNA.g0.sel <- rbinom(length(mat.noNA.g0), 1, cat.pr) # spatially random
              array.N[,,t+1][mat.noNA.g0[which(mat.noNA.g0.sel > 0)]] <- round((stoch.n.beta.func(length(mat.noNA.g0[which(mat.noNA.g0.sel > 0)]), M.cat, M.cat.sd)) * ((array.N[,,t+1])[mat.noNA.g0[which(mat.noNA.g0.sel > 0)]]), 0)
              cat.pr.est[t] <- length(mat.noNA.g0[which(mat.noNA.g0.sel > 0)])/dim(mat.noNA.g0.arr)[1]}
          }
          if (length(mat.noNA.g0[-mat.noNA.g0.sel]) > 0 & length(mat.noNA.g0) > 1 & spatial.cluster=="no") {
            array.N[,,t+1][mat.noNA.g0[which(mat.noNA.g0.sel > 0)]] <- round((stoch.n.beta.func(length(mat.noNA.g0[which(mat.noNA.g0.sel > 0)]), M.cat, M.cat.sd)) * ((array.N[,,t+1])[mat.noNA.g0[which(mat.noNA.g0.sel > 0)]]), 0)}
        }
        
        # apply NA layer post-growth
        array.N[,,t+1] <- ifelse(is.na(NA.array[,,t+1]) == T, NA, array.N[,,t+1])
        
        # apply additional mortality for cells with N < MVP.thresh
        if (small.pop.ext == "yes") {
          ltMVPthresh.sub <- which(array.N[,,t+1] < MVP.thresh) # which cells in array.N[,,t] < MVP.thresh
          ltMVPred.vec <- stoch.n.beta.func(length(ltMVPthresh.sub), (1-ltMVP.red), 0.05*ltMVP.red) # assume 5% SD
          array.N[,,t+1][ltMVPthresh.sub] <- round(array.N[,,t+1][ltMVPthresh.sub] * ltMVPred.vec, 0)
        }
        
        # if extinct, restart colonisation at origin
        if ((length(which(is.na(array.N[,,t+1]) == T | array.N[,,t+1] < 1) == T)) == i.rows*j.cols) {
          array.N[start.row1,start.col1,t+1] <- array.N[start.row1,start.col1,1]
        }
        
         if (print.map == "yes" & save.map == "no") {
          image(rot.mat(array.N[,,t+1]), col=rev(grey(1:100/100)), main=paste(format(entry.date-round(t*gen.l,0), big.mark = ","), " YBP",sep=""))
        }
        
        N.vec[t+1] <- sum(array.N[,,t+1], na.rm=T)
        poparea.vec[t+1] <- length(which(array.N[,,t+1] > 0)) * cell.area/1000
        pc.complete[t+1] <- (round(length(which(array.N[,,t+1] > 0, arr.ind=F)) / length(which(is.na(array.N[,,t+1]) != T, arr.ind=F)) * 99, 2))
        #print(paste("gen = ", t, "; ", round((t/gen.run*100), 1), "% compl; ", pc.complete[t+1], "% cells occ; ", round(as.numeric((proc.time() - proc.start)[3])/60, 1), " min elapsed", "; yrs proj = ", round(t*gen.l,0), "; ", "YBP = ", entry.date-round(t*gen.l,0), "; ", "Ntot = ", format(N.vec[t+1], big.mark=","), "; Pr(cat) = ", round(cat.pr.est[t], 3), sep=""))
        
      } # end t loop
      
      sub97.vec[m] <- (min(which(pc.complete >= 97), na.rm=T)) - 1
      N.mat[m,] <- N.vec

      #############################################################
      ## calculate date of first arrival per cell
      first.arrive <- array.N[,,1] * 0
      for (i in 1:dim(array.N)[1]) {
        for (j in 1:dim(array.N)[2]) {
          first.arrive[i,j] <- entry.date-round((which(array.N[i,j,] >= min.pop.size)[1])*gen.l,0)
        }
      }
      
      first.arrive.array[,,m] <- first.arrive
      
      print("*********************************************************************")
      print(paste("run = ", m, "; ", "time elapsed = ", round(as.numeric((proc.time() - proc.sim.start)[3])/60/60, 1), " hours", sep=""))
      print("*********************************************************************")
      
    } # end m reps loop

    
    ## average over first.arrive.array
    first.arrive.mean <- apply(first.arrive.array, c(1,2), mean, na.rm=T)    
    first.arrive.lo <-  apply(first.arrive.array, c(1,2), quantile, probs=0.025, na.rm=T)
    first.arrive.up <-  apply(first.arrive.array, c(1,2), quantile, probs=0.975, na.rm=T)
    
    image(rot.mat(first.arrive.array[,,1]), col=rev(grey(1:200/200)))
    
    first.arrive.list <- list()
    first.arrive.list$x <- unique(npp.sah$Lat)
    first.arrive.list$y <- unique(npp.sah$Lon)
    first.arrive.list$z <- first.arrive.mean
    first.arrive.rast.mn <- raster(first.arrive.list$z, xmn=range(first.arrive.list$x, na.rm=T)[1], xmx=range(first.arrive.list$x, na.rm=T)[2],
                                ymn=range(first.arrive.list$y, na.rm=T)[1], ymx=range(first.arrive.list$y, na.rm=T)[2],
                                crs=CRS("+proj=longlat +datum=WGS84"))
    par(xaxt="n", yaxt="n", bty="n")
    plot(first.arrive.rast.mn, col=rev(terrain.colors(18, alpha=1)), add=F, legend=F)
    cont.levels <- seq(entry.date, round(0.9*min(first.arrive.mean, na.rm=T),-3), round((entry.date - round(0.9*min(first.arrive.mean, na.rm=T),-3))/-20,-2))
    first.arrive.cont <- contour(first.arrive.rast.mn, levels=cont.levels, labels=cont.levels/1000, add=T,
                                 labcex=1.0, drawlabels=T, col="white", vfont=c("sans serif", "bold italic"))
    par(xaxt="s", yaxt="s", bty="o")
    
    first.arrive.list$z <- first.arrive.lo
    first.arrive.rast.lo <- raster(first.arrive.list$z, xmn=range(first.arrive.list$x, na.rm=T)[1], xmx=range(first.arrive.list$x, na.rm=T)[2],
                                   ymn=range(first.arrive.list$y, na.rm=T)[1], ymx=range(first.arrive.list$y, na.rm=T)[2],
                                   crs=CRS("+proj=longlat +datum=WGS84"))
    par(xaxt="n", yaxt="n", bty="n")
    plot(first.arrive.rast.lo, col=rev(terrain.colors(18, alpha=1)), add=F, legend=F)
    cont.levels <- seq(entry.date, round(0.9*min(first.arrive.lo, na.rm=T),-3), round((entry.date - round(0.9*min(first.arrive.lo, na.rm=T),-3))/-20,-2))
    first.arrive.cont <- contour(first.arrive.rast.lo, levels=cont.levels, labels=cont.levels/1000, add=T,
                                 labcex=1.0, drawlabels=T, col="white", vfont=c("sans serif", "bold italic"))
    par(xaxt="s", yaxt="s", bty="o")
    
    first.arrive.list$z <- first.arrive.up
    first.arrive.rast.up <- raster(first.arrive.list$z, xmn=range(first.arrive.list$x, na.rm=T)[1], xmx=range(first.arrive.list$x, na.rm=T)[2],
                                   ymn=range(first.arrive.list$y, na.rm=T)[1], ymx=range(first.arrive.list$y, na.rm=T)[2],
                                   crs=CRS("+proj=longlat +datum=WGS84"))
    par(xaxt="n", yaxt="n", bty="n")
    plot(first.arrive.rast.up, col=rev(terrain.colors(18, alpha=1)), add=F, legend=F)
    cont.levels <- seq(entry.date, round(0.9*min(first.arrive.up, na.rm=T),-3), round((entry.date - round(0.9*min(first.arrive.up, na.rm=T),-3))/-20,-2))
    first.arrive.cont <- contour(first.arrive.rast.up, levels=cont.levels, labels=cont.levels/1000, add=T,
                                 labcex=1.0, drawlabels=T, col="white", vfont=c("sans serif", "bold italic"))
    par(xaxt="s", yaxt="s", bty="o")
    
    writeRaster(first.arrive.rast.mn, filename=paste("SahSpread.mean.",name.add, ".1starrive.grd", sep=""), format="raster", overwrite=T)
    writeRaster(first.arrive.rast.lo, filename=paste("SahSpread.lo.",name.add, ".1starrive.grd", sep=""), format="raster", overwrite=T)
    writeRaster(first.arrive.rast.up, filename=paste("SahSpread.up.",name.add, ".1starrive.grd", sep=""), format="raster", overwrite=T)
    
    ## time to continental saturation (generations)
    median(sub97.vec, na.rm=T)
    quantile(sub97.vec, probs=0.025, na.rm=T)
    quantile(sub97.vec, probs=0.975, na.rm=T)
    
    median(sub97.vec, na.rm=T) * gen.l
    quantile(sub97.vec, probs=0.025, na.rm=T) * gen.l
    quantile(sub97.vec, probs=0.975, na.rm=T) * gen.l
    
    ## output mean + CI for N trajectory
    N.mean.vec <- apply(N.mat, 2, median, na.rm=T)
    N.lo.vec <- apply(N.mat, 2, quantile, probs=0.025, na.rm=T)
    N.up.vec <- apply(N.mat, 2, quantile, probs=0.975, na.rm=T)

    dat.N.out <- data.frame(N.mean.vec,N.up.vec,N.lo.vec)
    colnames(dat.N.out) <- c("Nm","Nup","Nlo")
    dat.out.name <- paste("ScenarioMean-Nout-",name.add,".csv",sep="")
    write.table(dat.N.out,file=dat.out.name,sep=",", row.names = F, col.names = T)
    
    save.image(paste("ScenarioMean.",name.add,".RData",sep=""))
    
    
