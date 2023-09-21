###
# Modelling a multi-optima OU evolution of
# niche traits on Phanaeini
# author: 'Viridiana Lizardo'


# Set Up
## Libraries

library(tidyverse)
library(bayou)
library(viridis)
library(sf)
library(tmap)
library(patchwork)
library(ggtree)
library(phytools)
source('C:/Users/USER/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/R projects/Theme_virix/theme_virix.r')

## palette
pal.inferno <- function(n){c('grey70', hcl.colors(n+1, palette = 'inferno')[2:n])}

# Load data
## Tree
tree <- ape::read.tree('phan.tree')


## Traits

trait_data <- read.csv('Appendix S1_ Niche data.csv',row.names = 1) %>% 
  as.matrix()

head(trait_data)


# Simulating a multi-optima OU process on the phylogeny

## Bio 12 log Precipitation

### Prior
priorOU_bio12_log <- make.prior(tree, 
                                dists=list(dalpha="dhalfcauchy", 
                                           dsig2="dhalfcauchy", 
                                           dk="cdpois", 
                                           dtheta="dnorm"),
                                param=list(dalpha=list(scale=0.1),
                                           dsig2=list(scale=0.1),
                                           dk=list(lambda=10, kmax=50), 
                                           dsb=list(bmax=1, prob=1), 
                                           dtheta=list(
                                             mean= log(mean(trait_data[,'mean_bio12'])),
                                             sd= log(1.5*mean(trait_data[,'sd_bio12']))))
)

startpars <- priorSim(priorOU_bio12_log, tree, plot=TRUE)$pars[[1]]
priorOU_bio12_log(startpars)



### Run `Bayou.makeMCMC`
nruns <- 1000000
set.seed(1)

mcmcOU_bio12_log <- bayou.makeMCMC(tree, log(trait_data[,'mean_bio12']), 
                         file.dir = 'output/bayou_bio12_LOG',
                         SE= log(trait_data[,'se_bio12']), 
                         prior=priorOU_bio12, 
                         outname="mean_bio12_log", plot.freq=NULL) # Set up the MCMC
mcmcOU_bio12_log$run(nruns) # Run the MCMC

chainOU_bio12_log <- mcmcOU_bio12_log$load()
ch <- chainOU_bio12_log
mc <- mcmcOU_bio12_log


name_ch <- 'Precipitation LOG'

pp.cutoff <- 0.2

chainOU_bio12_log <- set.burnin(chainOU_bio12_log, 0.3)

par(mfrow = c(3,2))
plot(chainOU_bio12_log, auto.layout=FALSE)

ch_sum<- summary(ch)

shift_sum <- shiftSummaries(ch, mc,
                            pp.cutoff = pp.cutoff,
                            branches = NULL)


## Bio 1. Annual Mean Temperature
### Prior

priorOU_bio1 <- make.prior(tree, 
                           dists=list(dalpha="dhalfcauchy", 
                                      dsig2="dhalfcauchy", 
                                      dk="cdpois", 
                                      dtheta="dnorm"),
                           param=list(dalpha=list(scale=0.1),
                                      dsig2=list(scale=0.1),
                                      dk=list(lambda=10, kmax=50), 
                                      dsb=list(bmax=1, prob=1), 
                                      dtheta=list(
                                        mean= mean(trait_data[,'mean_bio1']),
                                        sd=1.5*mean(trait_data[,'sd_bio1'])))
)


startpars <- priorSim(priorOU_bio1, tree, plot=TRUE)$pars[[1]]
priorOU_bio1(startpars)

### Run `Bayou.makeMCMC`

nruns <- 1000000
set.seed(1)
mcmcOU_bio1 <- bayou.makeMCMC(tree, trait_data[,'mean_bio1'], 
                              file.dir = 'output/bayou_bio1_new',
                              SE= trait_data[,'se_bio1'], 
                              prior=priorOU_bio1, 
                              outname="mean_bio1", plot.freq=NULL) # Set up the MCMC
mcmcOU_bio1$run(nruns) # Run the MCMC

chainOU_bio1 <- mcmcOU_bio1$load()
ch <- chainOU_bio1
mc <- mcmcOU_bio1

name_ch <- 'Temperature'

### Results

pp.cutoff <- 0.2

chainOU_bio1 <- set.burnin(chainOU_bio1, 0.3)

par(mfrow = c(3,2))
plot(chainOU_bio1, auto.layout=FALSE)

ch_sum<- summary(ch)
shift_sum <- shiftSummaries(ch, mc,
                            pp.cutoff = pp.cutoff,
                            branches = NULL)

