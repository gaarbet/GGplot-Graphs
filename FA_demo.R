
library(psych)

# simulate data

n <- 10000

population_structure <- matrix(c(
0.8	,	0	,	0	,	0	,
0.7	,	0	,	0	,	0	,
0.6	,	0	,	0	,	0	,
0.5	,	0	,	0	,	0	,
0.4	,	0	,	0	,	0	,
0	,	0.8	,	0	,	0	,
0	,	0.7	,	0	,	0	,
0	,	0.6	,	0	,	0	,
0	,	0.5	,	0	,	0	,
0	,	0.4	,	0	,	0	,
0	,	0	,	0.8	,	0	,
0	,	0	,	0.7	,	0	,
0	,	0	,	0.6	,	0	,
0	,	0	,	0.5	,	0	,
0	,	0	,	0.4	,	0	,
0	,	0	,	0	,	0.8	,
0	,	0	,	0	,	0.7	,
0	,	0	,	0	,	0.6	,
0	,	0	,	0	,	0.5	,
0	,	0	,	0	,	0.4	),20,4,byrow=TRUE)

population_phi <- matrix(c(
1	  ,	0.7	,	0.7	,	0.7	,
0.7	,	1	  ,	0.7	,	0.7	,
0.7	,	0.7	,	1	  ,	0.7	,
0.7	,	0.7	,	0.7	,	1	),4,4,byrow=TRUE)

dat <- sim.structure(fx=population_structure,Phi=population_phi,n=n,items=FALSE)$observed

# check number of factors

nfactors(dat)
fa.parallel(dat)
VSS.scree(cor(dat))

# run factor analyses

fa(dat)
fa.sort(fa(dat,2,rotate="oblimin"))
fa.sort(fa(dat,3,rotate="oblimin"))
fa.sort(fa(dat,4,rotate="oblimin"))
fa.sort(fa(dat,5,rotate="oblimin"))

omega(dat,2,flip=FALSE)
omega(dat,3,flip=FALSE)
omega(dat,4,flip=FALSE)
omega(dat,5,flip=FALSE)

omegaSem(dat,4)

