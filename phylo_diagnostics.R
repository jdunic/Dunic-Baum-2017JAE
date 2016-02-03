## This function analyzes standardized independent contrasts for
## a continuous trait on a phylogeny, and conducts a series of
## diagnostic tests to determine if the trait evolves consistent
## with a Brownian motion model. It is built on the ape library
## and ape format for representing phylogenies (object type phylo)
##
## David Ackerly, 2007
## dackerly@berkeley.edu
## 
## Some minor improvements by Nick Matzke, 2009-2011.
##

library(ape)
library(picante)

diagnostics = function(x, phy, scaled=TRUE, var.contrasts=TRUE, quadratic=NULL, calcK=TRUE,
    set.minbl=1, rescale.BL=FALSE, bl1=FALSE,logbl=FALSE, small.plot=FALSE, nrows=1, newplot=TRUE, diag5=TRUE) 
    {
        # x should be a trait matrix 
        # with dim = (#taxa, one trait)
        # and species as row.names
    
        #source("~/Documents/Projects/Toolbox/ape/nodeAge.R")
        #source("~/Documents/Projects/Toolbox/ape/K.R")
        #source("~/Documents/Projects/Toolbox/ape/pic3.R")
        
        if (small.plot) 
            {
            if (newplot) par(mfrow=c(nrows,6)) 
            } else par(mfcol=c(2,3))
    
        options('warn'=-1)
        xpos = 80
        ypos = 95
        
        if (small.plot) 
            {
            par(mar=c(1,1,1,1)) 
            }
        else
            {
            par(cex.lab=1.5,mar=c(5.1,5.1,4.1,2.1),cex.axis=1.5)
            }
        
        #tx = data.frame(t(x))
        minbl = min(phy$edge.length)
        if (rescale.BL) phy$edge.length = round(set.minbl*phy$edge.length/minbl)
        if (bl1) phy$edge.length=rep(1,length(phy$edge.length))
        if (logbl) phy$edge.length = log(phy$edge.length+1)+1
    
        pics = pic(x,phy,
            scaled=TRUE,var.contrasts=TRUE)
        
        pics = cbind(pics,ace(x,phy,method="pic")$ace)
    
        Npics = nrow(pics)  
        abspics = abs(pics[,1])
        
        #Diagnostic 1: are abspics drawn from 1/2 normal
        #Examine histogram and 1/2 normal probability plot
        if (small.plot) hist(abspics,main='',xlab='',ylab='',axes=FALSE)        else hist(abspics,main='')  
        dpics = c(-abspics,abspics)
        qs = qqnorm(dpics,plot.it=FALSE)
        if (small.plot) 
            {
            plot(qs$x[-(1:Npics)],abspics,pch=19,cex=1,         xlab='',
                ylab='',xaxt='n',yaxt='n')  
            }
        else 
            {
            plot(qs$x[-(1:Npics)],abspics,pch=19,cex=2,         xlab='normal proability quantiles',
                ylab='absolute value of contrasts')
            }
        
        #Diagnostic 2: are standardized pics independent of 
        # subtending branch lengths (standardization factor)
        picvsd = lm(abspics~I(sqrt(pics[,2])))
        if (small.plot)
            {
            plot(pics[,2],abspics,pch=19,cex=1,xlab='',ylab='',xaxt='n',yaxt='n')
            }
        else
            {
            plot(sqrt(pics[,2]),abspics,xlab='standard deviation of contrast',          ylab='absolute value of contrast',pch=19,cex=2)
            posx = min(sqrt(pics[,2])) + xpos*diff(range(sqrt(pics[,2])))/100
            posy = min(abspics) + ypos*diff(range(abspics))/100
            pval = paste('p = ',
                signif(summary(picvsd)$coeff[2,4],3),sep='')
            text(posx,posy,labels=pval,cex=1.5)
            }
        abline(picvsd)
        
        #Diagnostic 3: are pics independent of trait
        # values at nodes? Node values are calculated by 
        # contrast algorithm, which is not equal to ML 
        # ancestral states
        
        ancx = pics[,3]
        if (small.plot) 
            {
            plot(ancx,abspics,pch=19,cex=1,xlab='',ylab='',xaxt='n',yaxt='n') 
            }
        else
            {
            plot(ancx,abspics,pch=19,cex=2,xlab='node value',
                ylab='absolute value of contrast')
            }
        if (3 %in% quadratic) 
            {
            picvanc = lm(abspics~I(ancx^2)+ancx)
            qcoeff = summary(picvanc)$coeff
            qx = seq(min(ancx),max(ancx),length.out=100)
            qy = qcoeff[1,1] + qx*qcoeff[3,1] + qx^2*qcoeff[2,1]
            lines(qx,qy)
            }
        else
            {
            picvanc = lm(abspics~ancx)
            abline(picvanc)
            }   
        if (!small.plot)
            {
            pval = paste('p = ',
                signif(summary(picvanc)$coeff[2,4],3),sep='')
            posx = min(ancx) + xpos*diff(range(ancx))/100
            posy = min(abspics) + ypos*diff(range(abspics))/100
            text(posx,posy,labels=pval,cex=1.5)
            }
        
        #Diagnostic 4: are pics independent of node age
        # i.e., is there evidence of increasing or decreasing rates
        if (phy$edge[1,1] > 0)
            {
            Nterm = phy$edge[1,1] - 1
            }
        else
            {
            Nterm = max(phy$edge[,2])
            }
        int.nodes = which(phy$edge[,2] > Nterm)
        ages = node.age(phy)$ages
        int.ages = c(0,sort(ages[int.nodes]))
        picvage = lm(abspics~int.ages)
    
        if (small.plot)
            {
            plot(int.ages,abspics,pch=19,cex=1,xlab='',
                ylab='',xaxt='n',yaxt='n')
            }
        else
            {
            plot(int.ages,abspics,pch=19,cex=2,xlab='node age',
                ylab='absolute value of contrast')
            posx = min(int.ages) + xpos*diff(range(int.ages))/100
            posy = min(abspics) + ypos*diff(range(abspics))/100
            pval = paste('p = ',
                signif(summary(picvage)$coeff[2,4],3),sep='')
            text(posx,posy,labels=pval,cex=1.5)
            }
        abline(picvage)

        
        # Diagnostic 5: Plot of pairwise trait differences vs.
        #  cophenetic.phylo distances, and K value
        if (diag5) 
            {
            x2 = as.matrix(x[match(row.names(x),phy$tip.label)])
            row.names(x2)=row.names(x)[match(row.names(x),phy$tip.label)]
            xd = as.matrix(dist(x,diag=TRUE,upper=TRUE))
            phyd = cophenetic.phylo(phy)
            if (small.plot)
                {
                plot(phyd,xd,pch=19,cex=1,xlab='',
                    ylab='',xaxt='n',yaxt='n')
                } else
                {
                plot(phyd,xd,pch=19,cex=2,xlab='cophenetic distance',
                    ylab='pairwise trait difference')
                if (calcK) 
                    {
                    posx = min(phyd) + 20*diff(range(phyd))/100
                    posy = min(xd) + ypos*diff(range(xd))/100
                    x = as.matrix(x,ncol=1)
                    K = Kcalc(x,phy)    
                    text(posx,y=posy,
                        labels=paste('K =',signif(K,3),sep=' '),cex=ifelse              (small.plot,1.5,1.5))
                    }
                }
            }
        if (!small.plot) par(mfcol=c(1,1))  
        options('warn'=0)
    }
