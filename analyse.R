#!/usr/bin/env R
# Code to load the results and set up the data frame and create plots.
# Called from 'graphmulti.py'

args <- commandArgs(TRUE)
filename <- args[1]
data <- read.table(filename, sep = ",")
colnames(data) <- c("idnum", "size", "nwanteig", "delta", "forcamp", "omratio",
                    "premax", "postmax", "prepert", "postpert", "preiter",
                    "postiter", "failure", "minL", "maxL", "oldE", "newE",
                    "iter", "totnorm", "oldconc", "newconc")

library(plyr)
library(ggplot2)
library(reshape)
library(grid)
library(scales)

################################################################
# SUMMARY
################################################################

# FA <- forcing amplitude
# DE <- delta
# SI <- size
# NW <- nwanteigs
# OR <- omratio

## Summarize variables for energy difference.
# Summarize over one single variable.
ED_FA <- ddply(data, .(forcamp), function(x) c(prepert = mean(x$prepert), postpert = mean(x$postpert),
                            preiter = mean(x$preiter), positer = mean(x$postiter)))
ED_DE <- ddply(data, .(delta), function(x) c(prepert = mean(x$prepert), postpert = mean(x$postpert),
                            preiter = mean(x$preiter), positer = mean(x$postiter)))
ED_SI <- ddply(data, .(size), function(x) c(prepert = mean(x$prepert), postpert = mean(x$postpert),
                            preiter = mean(x$preiter), positer = mean(x$postiter)))
ED_NW <- ddply(data, .(nwanteig), function(x) c(prepert = mean(x$prepert), postpert = mean(x$postpert),
                            preiter = mean(x$preiter), positer = mean(x$postiter)))
ED_OR <- ddply(data, .(omratio), function(x) c(prepert = mean(x$prepert), postpert = mean(x$postpert),
                            preiter = mean(x$preiter), positer = mean(x$postiter)))

# Summarize over two variables.
ED_FA_DE <- ddply(data, c(.(forcamp), .(delta)),
              function(x) c(prepert = mean(x$prepert), postpert = mean(x$postpert),
                            preiter = mean(x$preiter), positer = mean(x$postiter)))
ED_FA_SI <- ddply(data, c(.(forcamp), .(size)),
              function(x) c(prepert = mean(x$prepert), postpert = mean(x$postpert),
                            preiter = mean(x$preiter), positer = mean(x$postiter)))
ED_DE_SI<- ddply(data, c(.(delta), .(size)),
              function(x) c(prepert = mean(x$prepert), postpert = mean(x$postpert),
                            preiter = mean(x$preiter), positer = mean(x$postiter)))

# Summarize over three variables.
ED_FA_DE_SI <-  ddply(data, c(.(forcamp), .(size), .(delta)),
              function(x) c(prepert = mean(x$prepert), postpert = mean(x$postpert),
                            preiter = mean(x$preiter), positer = mean(x$postiter)))

## Summarize variables for maximum voltage.
# Summarize over one single variable.
V_FA <- ddply(data, .(forcamp), function(x) c(premax = mean(x$premax), postmax = mean(x$postmax)))
V_DE <- ddply(data, .(delta), function(x) c(premax = mean(x$premax), postmax = mean(x$postmax)))
V_SI <- ddply(data, .(size), function(x) c(premax = mean(x$premax), postmax = mean(x$postmax)))
V_NW <- ddply(data, .(nwanteig), function(x) c(premax = mean(x$premax), postmax = mean(x$postmax)))
V_OR <- ddply(data, .(omratio), function(x) c(premax = mean(x$premax), postmax = mean(x$postmax)))

# Summarize over two variables.
V_FA_DE <- ddply(data, c(.(forcamp), .(delta)),
              function(x) c(premax = mean(x$premax), postmax = mean(x$postmax)))
V_FA_SI <- ddply(data, c(.(forcamp), .(size)),
              function(x) c(premax = mean(x$premax), postmax = mean(x$postmax)))
V_DE_SI<- ddply(data, c(.(delta), .(size)),
              function(x) c(premax = mean(x$premax), postmax = mean(x$postmax)))

# Summarize over three variables.
V_FA_DE_SI <-  ddply(data, c(.(forcamp), .(size), .(delta)),
              function(x) c(premax = mean(x$premax), postmax = mean(x$postmax)))

## Summarize variables for inductor values.
# Summarize over one single variable.
L_FA <- ddply(data, .(forcamp), function(x) c(minL = mean(x$minL), maxL = mean(x$maxL)))
L_DE <- ddply(data, .(delta), function(x) c(minL = mean(x$minL), maxL = mean(x$maxL)))
L_SI <- ddply(data, .(size), function(x) c(minL = mean(x$minL), maxL = mean(x$maxL)))
L_NW <- ddply(data, .(nwanteig), function(x) c(minL = mean(x$minL), maxL = mean(x$maxL)))
L_OR <- ddply(data, .(omratio), function(x) c(minL = mean(x$minL), maxL = mean(x$maxL)))

# Summarize over two variables.
L_FA_DE <- ddply(data, c(.(forcamp), .(delta)),
              function(x) c(minL = mean(x$minL), maxL = mean(x$maxL)))
L_FA_SI <- ddply(data, c(.(forcamp), .(size)),
              function(x) c(minL = mean(x$minL), maxL = mean(x$maxL)))
L_DE_SI<- ddply(data, c(.(delta), .(size)),
              function(x) c(minL = mean(x$minL), maxL = mean(x$maxL)))

# Summarize over three variables.
L_FA_DE_SI <-  ddply(data, c(.(forcamp), .(size), .(delta)),
              function(x) c(minL = mean(x$minL), maxL = mean(x$maxL)))

## Summarize variables for mean energy values
# Summarize over one single variable.
MEN_FA <- ddply(data, .(forcamp), function(x) c(oldE = mean(x$oldE), newE = mean(x$newE)))
MEN_DE <- ddply(data, .(delta), function(x) c(oldE = mean(x$oldE), newE = mean(x$newE)))
MEN_SI <- ddply(data, .(size), function(x) c(oldE = mean(x$oldE), newE = mean(x$newE)))
MEN_NW <- ddply(data, .(nwanteig), function(x) c(oldE = mean(x$oldE), newE = mean(x$newE)))
MEN_OR <- ddply(data, .(omratio), function(x) c(oldE = mean(x$oldE), newE = mean(x$newE)))

# Summarize over two variables.
MEN_FA_DE <- ddply(data, c(.(forcamp), .(delta)),
              function(x) c(oldE = mean(x$oldE), newE = mean(x$newE)))
MEN_FA_SI <- ddply(data, c(.(forcamp), .(size)),
              function(x) c(oldE = mean(x$oldE), newE = mean(x$newE)))
MEN_DE_SI<- ddply(data, c(.(delta), .(size)),
              function(x) c(oldE = mean(x$oldE), newE = mean(x$newE)))

# Summarize over three variables.
MEN_FA_DE_SI <-  ddply(data, c(.(forcamp), .(size), .(delta)),
              function(x) c(oldE = mean(x$oldE), newE = mean(x$newE)))

## Summarize variables for failure rates.
# Summarize over one single variable.
FAIL_FA <- ddply(data, .(forcamp), function(x) c(failure = mean(x$failure)))
FAIL_DE <- ddply(data, .(delta), function(x) c(failure = mean(x$failure)))
FAIL_SI <- ddply(data, .(size), function(x) c(failure = mean(x$failure)))
FAIL_NW <- ddply(data, .(nwanteig), function(x) c(failure = mean(x$failure)))
FAIL_OR <- ddply(data, .(omratio), function(x) c(failure = mean(x$failure)))

# Summarize over two variables.
FAIL_FA_DE <- ddply(data, c(.(forcamp), .(delta)),
              function(x) c(failure = mean(x$failure)))
FAIL_FA_SI <- ddply(data, c(.(forcamp), .(size)),
              function(x) c(failure = mean(x$failure)))
FAIL_DE_SI<- ddply(data, c(.(delta), .(size)),
              function(x) c(failure = mean(x$failure)))

# Summarize over three variables.
FAIL_FA_DE_SI <-  ddply(data, c(.(forcamp), .(size), .(delta)),
              function(x) c(failure = mean(x$failure)))

## Summarize variables for number of newton iterations.
# Summarize over one single variable.
ITR_FA <- ddply(data, .(forcamp), function(x) c(iter = mean(x$iter)))
ITR_DE <- ddply(data, .(delta), function(x) c(iter = mean(x$iter)))
ITR_SI <- ddply(data, .(size), function(x) c(iter = mean(x$iter)))
ITR_NW <- ddply(data, .(nwanteig), function(x) c(iter = mean(x$iter)))
ITR_OR <- ddply(data, .(omratio), function(x) c(iter = mean(x$iter)))

# Summarize over two variables.
ITR_FA_DE <- ddply(data, c(.(forcamp), .(delta)),
              function(x) c(iter = mean(x$iter)))
ITR_FA_SI <- ddply(data, c(.(forcamp), .(size)),
              function(x) c(iter = mean(x$iter)))
ITR_DE_SI<- ddply(data, c(.(delta), .(size)),
              function(x) c(iter = mean(x$iter)))

# Summarize over three variables.
ITR_FA_DE_SI <-  ddply(data, c(.(forcamp), .(size), .(delta)),
              function(x) c(iter = mean(x$iter)))

## Summarize variables for failure rates.
# Summarize over one single variable.
TOT_FA <- ddply(data, .(forcamp), function(x) c(totnorm = mean(x$totnorm)))
TOT_DE <- ddply(data, .(delta), function(x) c(totnorm = mean(x$totnorm)))
TOT_SI <- ddply(data, .(size), function(x) c(totnorm = mean(x$totnorm)))
TOT_NW <- ddply(data, .(nwanteig), function(x) c(totnorm = mean(x$totnorm)))
TOT_OR <- ddply(data, .(omratio), function(x) c(totnorm = mean(x$totnorm)))

# Summarize over two variables.
TOT_FA_DE <- ddply(data, c(.(forcamp), .(delta)),
              function(x) c(totnorm = mean(x$totnorm)))
TOT_FA_SI <- ddply(data, c(.(forcamp), .(size)),
              function(x) c(totnorm = mean(x$totnorm)))
TOT_DE_SI<- ddply(data, c(.(delta), .(size)),
              function(x) c(totnorm = mean(x$totnorm)))

# Summarize over three variables.
TOT_FA_DE_SI <-  ddply(data, c(.(forcamp), .(size), .(delta)),
              function(x) c(totnorm = mean(x$totnorm)))

EC_FA <- ddply(data, .(forcamp), function(x) c(oldconc = mean(x$oldconc), newconc = mean(x$newconc)))
EC_DE <- ddply(data, .(delta), function(x) c(oldconc = mean(x$oldconc), newconc = mean(x$newconc)))
EC_SI <- ddply(data, .(size), function(x) c(oldconc = mean(x$oldconc), newconc = mean(x$newconc)))
EC_NW <- ddply(data, .(nwanteig), function(x) c(oldconc = mean(x$oldconc), newconc = mean(x$newconc)))
EC_OR <- ddply(data, .(omratio), function(x) c(oldconc = mean(x$oldconc), newconc = mean(x$newconc)))

# Summarize over two variables.
EC_FA_DE <- ddply(data, c(.(forcamp), .(delta)),
              function(x) c(oldconc = mean(x$oldconc), newconc = mean(x$newconc)))
EC_FA_SI <- ddply(data, c(.(forcamp), .(size)),
              function(x) c(oldconc = mean(x$oldconc), newconc = mean(x$newconc)))
EC_DE_SI<- ddply(data, c(.(delta), .(size)),
              function(x) c(oldconc = mean(x$oldconc), newconc = mean(x$newconc)))

# Summarize over three variables.
EC_FA_DE_SI <-  ddply(data, c(.(forcamp), .(size), .(delta)),
              function(x) c(oldconc = mean(x$oldconc), newconc = mean(x$newconc)))


################################################################
# PLOTTING
################################################################
data.m <- melt(EC_FA_DE_SI, measure.vars = c("oldconc", "newconc"))
Eplot <- ggplot(data.m, aes(x = delta, y = value,
                             color = factor(forcamp),
                             shape = factor(variable)))
Eplot <- Eplot + geom_line(aes(group = interaction(variable, forcamp),
                                 color = factor(forcamp),
                                 shape = factor(variable)),
                             size = 0.5)
Eplot <- Eplot + geom_point(aes(group = interaction(variable, forcamp),
                                  color = factor(forcamp),
                                  shape = factor(variable)),
                              size = 2.0)
Eplot <- Eplot + facet_grid(~ size)

xlab_formatter <- function(x){
  lab <- gsub('^0', '', x)
}
Eplot <-  Eplot + scale_x_continuous(breaks = c(0.25, 0.5, 0.75), label=xlab_formatter)
Eplot <-  Eplot + scale_y_continuous(label=percent)
Eplot <-  Eplot + theme(axis.text.x = element_text(colour = 'black', size = 9, angle = 90))
Eplot <-  Eplot + theme(axis.text.y = element_text(colour = 'black', size = 9))
Eplot <-  Eplot + theme(axis.title.x = element_text(colour = 'black', angle = 0, size = 10))
Eplot <-  Eplot + theme(axis.title.y = element_text(colour = 'black', angle = 90, size = 10))
Eplot <-  Eplot + theme(axis.ticks.margin = unit(0.25, "lines"))
Eplot <-  Eplot + theme(axis.ticks.length = unit(0.25, "lines"))

Eplot <- Eplot + labs(shape = "Simulation")
Eplot <- Eplot + scale_shape_discrete(name = "Simulation", breaks = c("oldconc","newconc"), labels = c(expression(kappa[pre]), expression(kappa[post])))
Eplot <- Eplot + scale_color_discrete(name = "Amplitude")
Eplot <- Eplot + xlab(expression(delta))
Eplot <- Eplot + ylab("% Energy in higher harmonics")
Eplot <- Eplot + theme(legend.text = element_text(colour = "black",
                          size = 9))
Eplot <- Eplot + theme(legend.title = element_text(colour = "black",
                          size = 9))
Eplot <- Eplot + theme(panel.margin=unit(1, "lines"))
filebase <- strsplit(filename,"[.]")[[1]][1]
efile <- paste(filebase,'amplitude.eps')
ggsave(Eplot, file = efile, width = 6.0, height = 2.0, units = "in", dpi = 500)

################################################################
#### Voltage plot.
################################################################
data.m <- melt(V_FA_DE_SI, measure.vars = c("premax", "postmax"))
Vplot <- ggplot(data.m, aes(x = delta, y = value,
                             color = factor(forcamp),
                             shape = factor(variable)))
Vplot <- Vplot + geom_line(aes(group = interaction(variable, forcamp),
                                 color = factor(forcamp),
                                 shape = factor(variable)),
                             size = 0.5)
Vplot <- Vplot + geom_point(aes(group = interaction(variable, forcamp),
                                  color = factor(forcamp),
                                  shape = factor(variable)),
                              size = 2.0)
Vplot <- Vplot + facet_grid(~ size)

xlab_formatter <- function(x){
  lab <- gsub('^0', '', x)
}

Vplot <-  Vplot + scale_x_continuous(breaks = c(0.25, 0.5, 0.75), label=xlab_formatter)
Vplot <-  Vplot + theme(axis.text.x = element_text(colour = 'black', size = 9, angle = 90))
Vplot <-  Vplot + theme(axis.text.y = element_text(colour = 'black', size = 9))
Vplot <-  Vplot + theme(axis.title.x = element_text(colour = 'black', angle = 0, size = 10))
Vplot <-  Vplot + theme(axis.title.y = element_text(colour = 'black', angle = 90, size = 10))
Vplot <-  Vplot + theme(axis.ticks.margin = unit(0.25, "lines"))
Vplot <-  Vplot + theme(axis.ticks.length = unit(0.25, "lines"))

Vplot <- Vplot + labs(shape = "Simulation")
Vplot <- Vplot + scale_shape_discrete(name = "Simulation", breaks = c("premax", "postmax"), labels = c(expression(V[pre]), expression(V[post])))
Vplot <- Vplot + scale_color_discrete(name = "Amplitude")
Vplot <- Vplot + theme(legend.text = element_text(colour = "black",
                          size = 9))
Vplot <- Vplot + theme(legend.title = element_text(colour = "black",
                          size = 9))

Vplot <- Vplot + xlab(expression(delta))
Vplot <- Vplot + ylab("Maximum voltage")
Vplot <- Vplot + theme(panel.margin=unit(1, "lines"))
filebase <- strsplit(filename,"[.]")[[1]][1]
vfile <- paste(filebase,'voltage.eps')
ggsave(Vplot, file = vfile, width = 6.0, height = 2.0, units = "in", dpi = 500)
