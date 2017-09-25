# This is the model that's failing to converge on one of two machines I use

library(lmerTest)
sessionInfo()

load('dat.RData')

###### this is all setup ####
dat$DM1 = dummy(dat$GroupID, 'SBS')
dat$DM2 = dummy(dat$GroupID, 'SBS-IF')
dat$obs = factor(seq(nrow(dat)))

c0 = c(1/3,1/3,1/3)
c1 = c(0,-1,1)
c2 = c(-1, -0.5, -0.5)
m = rbind(c1,c2,c0)		# c0 gets dropped
cntr = solve(m)
###### done setup ####

m = lmer(asin(sqrt(Lachnospiraceae)) ~ GroupID +
										(1|PatientID) +
										(0+DM1|obs) +
										(0+DM2|obs),
			data = dat,
			contrasts = list(GroupID = cntr),
			control = lmerControl(check.nobs.vs.nRE = 'ignore',
									check.nobs.vs.nlev = 'ignore')
	)

summary(m)
