#!/usr/bin/env R

library(plyr)
library(reshape)
library(ggplot2)

mydata <- read.table('fdecay.txt',sep='\t',header=FALSE)
names(mydata) <- c("Iter","Pert","Numerical")
dims <- dim(mydata)
newdata <- melt(mydata)
newdata$x <- c(1L:dims[1], 1L:dims[1],1L:dims[1])
newdata$value <- newdata$value

myplot <- ggplot(newdata, aes(x=x,y=value,colour=variable,shape=variable))
myplot <- myplot + geom_line(size=0.5) + geom_point(size=2.0)

myplot <- myplot + scale_color_discrete(name="Method",
       breaks=c("Iter","Pert","Numerical"),
       labels=c("Iter","Pert","Numerical"))
myplot <- myplot + scale_shape_discrete(name="Method",
       breaks=c("Iter","Pert","Numerical"),
       labels=c("Iter","Pert","Numerical"))

myplot <- myplot + xlab("Modes") + ylab("Log-Norm")
myplot <- myplot + theme(legend.position="top")
myplot <- myplot + theme(panel.grid.major = element_blank())
myplot <- myplot + theme(panel.grid.minor = element_blank())
myplot <- myplot + theme(panel.background = element_blank())
myplot <- myplot + theme(axis.line = element_line(colour = "black"))

mywidth = 3.5
myheight = 4
ggsave(myplot,file='fourierdecay.eps',width=mywidth,height=myheight,units="in",dpi=500)
