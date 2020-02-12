
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

difficulties <- seq(-2,2,length.out=20)

dat <- sim.structure(fx=population_structure,Phi=population_phi,n=n,items=TRUE,cat=1,d=difficulties)$observed

# check number of factors

xcor <- polychoric(dat)$rho

nfactors(xcor,n.obs=n)a
fa.parallel.poly(dat)
VSS.scree(xcor)

# run factor analyses

fa(dat)
fa.sort(fa(xcor,2,rotate="oblimin"))
fa.sort(fa(xcor,3,rotate="oblimin"))
fa.sort(fa(xcor,4,rotate="oblimin"))
fa.sort(fa(xcor,5,rotate="oblimin"))
#fa.sort(irt.fa(dat,4,rotate="oblimin")$fa)

omega(xcor,2,flip=FALSE)
omega(xcor,3,flip=FALSE)
omega(xcor,4,flip=FALSE)
omega(xcor,5,flip=FALSE)

omegaSem(xcor,4)