# Make the data set

dat = data.frame(Genus = 'Escherichia',
					SP1 = rep(10,10),
					SP2 = rep(5,10),
					SP3 = rep(17,10))

# Make the data frame long. This works
dat %>% gather(Species,Count,2:4) -> dat_long

# Specify columns more programmatically. This used to work but now doesn't.

st = 2
dat %>% gather(Species,Count,st:ncol(dat)) -> dat_long

# This also fails, but more excitingly:

st = 1
dat %>% gather(Species,Count,(st+1):ncol(dat)) -> dat_long
