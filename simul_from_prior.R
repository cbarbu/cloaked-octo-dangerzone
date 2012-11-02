source("RanalysisFunctions.R") # at least for set_to
source("functions_migration_simul.R") # for gillespie code

maps.tot<-read.csv("maps_hunter_blocks.csv")

# focusing on hunter for now
maps<-maps.tot[which(maps.tot$D==7 & maps.tot$X>226000 & maps.tot$Y>8179000 & maps.tot$Y<8180000),]
blockIndex <- maps.tot$block_num

# calculate the distances between all houses
dist_out <- makeDistClasses(as.vector(maps$X), as.vector(maps$Y), c(0, 100))
dist_mat <- dist_out$dist_mat

# there were some NAs in block index (corrects by making NA = new, different blocks)
if(any(is.na(blockIndex))){
	bad <- which(is.na(blockIndex))
	blockIndex[bad] <- max(blockIndex[-bad])+(1:length(bad))
}

# plot the initial infested, probability infestation day before spray, infestation day after spray
par(mfrow=c(1,3))
plot_reel(maps$X,maps$Y,maps$infested,base=0)
plot_reel(maps$X,maps$Y,maps$est_p.i_db,base=0)
plot_reel(maps$X,maps$Y,maps$est_p.i_da,base=0,top=max(maps$est_p.i_da))

# zm() here helps

## import reports
## not doing anything with reports yet
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


## set the parameters for the simulation
## need to be fixed with actual fitted parameters instead of dummy parameters
halfDistJ<-50
halfDistH<-10
rateMove<-0.0375
useDelta<-TRUE # if TRUE, weightSkipInMove<-0
delta<- 1
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

## make the cumulative probability matrix
cumulProbMat <- generate_prob_mat_C(halfDistJ, halfDistH, useDelta, delta, rateHopInMove, rateSkipInMove, rateJumpInMove, dist_out$dists, blockIndex, cumul = TRUE)

## set the seed
## seed <- runif(1, 1, 2^31-1)
seed <- 31415926535

## run the gillespie simulation with the priors (maps$est_p.i_da)
out <- simul_priors_gillespie(prob_inf_vec = maps$est_p.i_da, blockIndex = blockIndex, endTime = 3*52, Nrep = 100, rateMove, seed, cumulProbMat)

infestedDens <- out$infestedDens
ageMat <- out$ageMat
infestedMat <- out$infestedMat

## plot infestation density
dev.new()
par(mfrow=c(1,2))
plot_reel(maps$X,maps$Y,maps$est_p.i_da,base=0)
plot_reel(maps$X,maps$Y,infestedDens,base=0)

stop()

#==================
# Time based plotting
# needs further refinement
#==================

# initInf <- which(ageMat[2, ] == 0)
# infHouse <- infestedMat[2, initInf]
# keep <- (tail(initInf)[6]+1):(head(which(ageMat[2, ]==-1))[1]-1)
# infested <- rep(0, length(infestedDens))
# infested[infHouse] <- 1
# plot_reel(maps$X, maps$Y, infested, base = 0)

# Sys.sleep(2)

# for(n in keep){
# Sys.sleep(0.1)
# infested[infestedMat[2, n]] <- 1
# plot_reel(maps$X, maps$Y, infested, base = 0)

# }

# dev.new()
# Sys.sleep(1)
# for(n in seq(1, 7, 0.5)){

#	out <- simul_priors_gillespie(prob_inf_vec = maps$est_p.i_da, blockIndex = blockIndex, endTime = n*52, Nrep = 20, rateMove, seed, halfDistJ = halfDistJ, halfDistH = halfDistH, useDelta = useDelta, delta = delta, rateHopInMove = rateHopInMove, rateSkipInMove = rateSkipInMove, rateJumpInMove = rateJumpInMove, dist_out = dist_out$dists)

#	infestedDens <- out$infestedDens
	
#	plot_reel(maps$X,maps$Y,infestedDens,base=0)
#	Sys.sleep(0.25)

# }



