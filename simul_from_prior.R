source("RanalysisFunctions.R") # at least for set_to
source("functions_migration_simul.R") # for gillespie code

maps.tot<-read.csv("roc_p.i_fromIadjustedwithII.csv")
blocks<-read.csv("maps_hunter_blocks.csv")

# focusing on hunter for now
maps<-maps.tot[which(maps.tot$D==7 & maps.tot$X>226000 & maps.tot$Y>8179000 & maps.tot$Y<8180000),]

par(mfrow=c(1,3))
plot_reel(maps$X,maps$Y,maps$infested,base=0)
plot_reel(maps$X,maps$Y,maps$est_p.i_db,base=0)
plot_reel(maps$X,maps$Y,maps$est_p.i_da,base=0,top=max(maps$est_p.i_da))

# zm() here helps

## import reports
den<-read.csv("DENUNCIAS_2012Jun12.csv",sep=";")
den<-set_to(den)
den<-den[which(den$L!=0),]

den$infestedDen<- (den$ADULTOS>0 | den$NINFAS>0)
den$colonizedDen<- (den$NINFAS>0)
den$unicode_gps<-make_unicode_gps(den$UNICODE)

# get denuncias in maps
maps$den<-as.numeric(maps$unicode %in% den$unicode_gps)

# get positive denuncias in maps
maps$denPos<-as.numeric(maps$unicode %in% den$unicode_gps[which(den$infestedDen)])

# plot denuncias on top of previous
par(mfrow=c(1,3))
plot_reel(maps$X,maps$Y,maps$infested,base=0)
with(maps[which(maps$denPos==1),],lines(X,Y,type="p",pty=16,col="blue"))
plot_reel(maps$X,maps$Y,maps$est_p.i_db,base=0)
with(maps[which(maps$denPos==1),],lines(X,Y,type="p",pty=16,col="blue"))
plot_reel(maps$X,maps$Y,maps$est_p.i_da,base=0,top=max(maps$est_p.i_da))
with(maps[which(maps$denPos==1),],lines(X,Y,type="p",pty=16,col="blue"))

halfDistJ<-50
halfDistH<-10
rateMove<-0.0375
useDelta<-TRUE # if TRUE, weightSkipInMove<-0
delta<- 0.25
set.seed(777)

weightHopInMove<-1 # should not be changed, the ref for the two others
weightSkipInMove<-1 # ignored if useDelta == TRUE
weightJumpInMove<-0.07

# set corresponding rates
if (useDelta)
{
	weightSkipInMove<-0 
}

totalWeight<-(weightHopInMove+weightSkipInMove+weightJumpInMove)
rateHopInMove<-weightHopInMove/totalWeight
rateSkipInMove<-weightSkipInMove/totalWeight
rateJumpInMove<-weightJumpInMove/totalWeight # 1 - rateHopInMove - rateSkipInMove

seed <- runif(1, 1, 2^31-1)

dist_out <- makeDistClasses(as.vector(maps$X), as.vector(maps$Y), c(0, 100))

blockIndex <- blocks$block_num
bad <- which(is.na(blockIndex))
blockIndex[bad] <- max(blockIndex[-bad])+(1:length(bad))
infestedDens <- simul_priors_gillespie(prob_inf_vec = maps$est_p.i_da, blockIndex = blockIndex, endTime = 7*52, Nrep = 100, rateMove, seed, halfDistJ = halfDistJ, halfDistH = halfDistH, useDelta = useDelta, delta = delta, rateHopInMove = rateHopInMove, rateSkipInMove = rateSkipInMove, rateJumpInMove = rateJumpInMove, dist_out = dist_out$dists)

par(mfrow=c(1,2))
plot_reel(maps$X,maps$Y,maps$est_p.i_da,base=0)
plot_reel(maps$X,maps$Y,infestedDens,base=0)

