#!/usr/bin/env R

library(plyr)
library(reshape)
library(ggplot2)

alldata <- read.table('fixedpoint.txt',sep='\t',header=FALSE,skip=0)
mydata <- alldata[2:(dim(alldata)[1]),]
numerr <- alldata[1,2]
names(mydata) <- c("Perturbative","Iterative")
dims <- dim(mydata)
newdata <- melt(mydata)
newdata$x <- c(1L:dims[1], 1L:dims[1])
newdata$value <- newdata$value
horizontal <- numerr

myplot <- ggplot(newdata, aes(x=x,y=value,colour=variable,shape=variable))
myplot <- myplot + geom_line(size=0.5) + geom_point(size=2.0)

myplot <- myplot + scale_color_discrete(name="Method",
       breaks=c("Perturbative","Iterative"),
       labels=c("Perturbative","Iterative"))
myplot <- myplot + scale_shape_discrete(name="Method",
       breaks=c("Perturbative","Iterative"),
       labels=c("Perturbative","Iterative"))

myplot <- myplot + xlab("Iterations") + ylab("Log-Norm difference")
myplot <- myplot + theme(legend.position="top")

myplot <- myplot + geom_hline(aes(yintercept=horizontal))
myplot <- myplot + theme(panel.grid.major = element_blank())
myplot <- myplot + theme(panel.grid.minor = element_blank())
myplot <- myplot + theme(panel.background = element_blank())
myplot <- myplot + theme(axis.line = element_line(colour = "black"))

ggsave(myplot,file='fixedpointerror.eps',width=3.5,height=4,units="in",dpi=500)
