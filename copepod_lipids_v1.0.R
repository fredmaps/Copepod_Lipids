###################################################
#                                                 #
# Script to compute potential diapause duration   #
# from individual copepod size and lipid content  #
#                                                 #
# Schmid, Maps, Fortier 2017                      #
#                                                 #
###################################################


#--- Read data
#setwd("~/mydirectory")
data <- read.csv("copepod_lipids.csv")


#--- Define working variabes & indices

# Flag to get ALL the diagnostic plots
print( "Do you want ALL diagnostic plots ? [TRUE|FALSE] " ); Pflag <- scan(,what="logical")

# species
SP <- data$spec

Ch <- SP=="C. hyperboreus"
Cg <- SP=="C. glacialis"
Ml <- SP=="M. longa"

# stages
SG <- data$stage

C3 <- SG=="C3"
C4 <- SG=="C4"
C5 <- SG=="C5"
CF <- SG=="F"

# vertical casts
HL <- data$haul

Hd <- HL=="18:45h"
Hn <- HL=="2:40h"

# other variables
Z  <- data$depth

PL <- data$prosomelength

TL <- data$totallipid

LF <- data$fullnessratio


#****************************************************************************************
# Figure 1: prosome length frequency distribution / stage

grey  <- rgb(0.3,0.3,0.3,0.5)
blue  <- rgb(0,0,1,0.5)
green <- rgb(0,1,0,0.5)
red   <- rgb(1,0,0,0.5)

# Optional diagnostic plot
if (Pflag) {
  print( "Which species do you want ? C. hyperboreus: Ch | C. glacialis: Cg | M. longa: Ml " ); Sflag <- scan(,what="character")
  if (Sflag=="Ch") {
    Psp   <- Ch
    Pmain <- "C. hyperboreus"
    hist(PL[Psp&C3], xlim=c(1,8), ylim=c(0,40), col=grey, lwd=2, main="", xlab="")
    par(new=TRUE)
    hist(PL[Psp&C4], xlim=c(1,8), ylim=c(0,40), col=blue,     lwd=2, main="", xlab="")
    par(new=TRUE)
  } else if (Sflag=="Cg") {
    Psp   <- Cg
    Pmain <- "C. glacialis"
    hist(PL[Psp&C4], xlim=c(1,8), ylim=c(0,40), col=blue,     lwd=2, main="", xlab="")
    par(new=TRUE)
  } else if (Sflag=="Ml") {
    Psp   <- Ml
    Pmain <- "M. longa"
  }
  hist(PL[Psp&C5], xlim=c(1,8), ylim=c(0,40), col=green,    lwd=2, main=Pmain, xlab="Prosome length (mm)")
  par(new=TRUE)
  hist(PL[Psp&CF], xlim=c(1,8), ylim=c(0,40), col=red,      lwd=2, main="",    xlab="")
  legend("topright", legend=c("CF","C5","C4","C3"), pch=15, col=c(red,green,blue,grey), pt.cex=3, cex=1.5, bty="n")
}

boxplot(PL[Ch]~SG[Ch], ylim=c(1,7), col=c(grey,blue,green,red), main="Prosome length (mm)", bty="n")
par(new=TRUE)
boxplot(PL[Cg]~SG[Cg], ylim=c(1,7), col=c(grey,blue,green,red), main="")
par(new=TRUE)
boxplot(PL[Ml]~SG[Ml], ylim=c(1,7), col=c(grey,blue,green,red), main="")
legend("topleft", legend=c("CF","C5","C4","C3"), pch=15, col=c(red,green,blue,grey), pt.cex=3, cex=1.5, bty="n")

#****************************************************************************************


#****************************************************************************************
# Figure 2: lipid fullness frequency distribution / stage

# Optional diagnostic plot
if (Pflag) {
  print( "Which species do you want ? C. hyperboreus: Ch | C. glacialis: Cg | M. longa: Ml " ); Sflag <- scan(,what="character")
  if (Sflag=="Ch") {
    Psp   <- Ch
    Pmain <- "C. hyperboreus"
    hist(LF[Psp&C3], xlim=c(0,1), ylim=c(0,50), col=grey, lwd=2, main="", xlab="")
    par(new=TRUE)
    hist(LF[Psp&C4], xlim=c(0,1), ylim=c(0,50), col=blue,     lwd=2, main="", xlab="")
    par(new=TRUE)
  } else if (Sflag=="Cg") {
    Psp   <- Cg
    Pmain <- "C. glacialis"
    hist(LF[Psp&C4], xlim=c(0,1), ylim=c(0,50), col=blue,     lwd=2, main="", xlab="")
    par(new=TRUE)
  } else if (Sflag=="Ml") {
    Psp   <- Ml
    Pmain <- "M. longa"
  }
  hist(LF[Psp&C5], xlim=c(0,1), ylim=c(0,50), col=green, main=Pmain, xlab="Lipid fullness")
  par(new=TRUE)
  hist(LF[Psp&CF], xlim=c(0,1), ylim=c(0,50), col=red,   main="")
  legend("topright", legend=c("CF","C5","C4","C3"), pch=15, col=c(red,green,blue,grey), pt.cex=3, cex=1.5, bty="n")
}

boxplot(LF[Ch]~SG[Ch], ylim=c(0.1,0.9), col=c(grey,blue,green,red), main="Lipid fullness")
par(new=TRUE)
boxplot(LF[Cg]~SG[Cg], ylim=c(0.1,0.9), col=c(grey,blue,green,red), main=" ")
par(new=TRUE)
boxplot(LF[Ml]~SG[Ml], ylim=c(0.1,0.9), col=c(grey,blue,green,red), main=" ")
legend("topleft", legend=c("CF","C5","C4","C3"), pch=15, col=c(red,green,blue,grey), pt.cex=3, cex=1.5, bty="n")

#****************************************************************************************


#--- Conversion from TL (mg) to LC (ug) using wax esters data from Kattner & Graeve 1991 (doi:10.3402/polar.v10i2.6760)
#    Almost all lipids are WE, and the most common WE length chain in both Calanus congeners is C:36
#    Lipid mass = nC + 2nH + 2O with C = carbon M (12), H = hydrogen M (1), O = oxygen M (16)

C2L <- (36*12) / (36*12 + 2*36 + 2*16)

LC  <- 1e3 * TL*C2L


#****************************************************************************************
# Figure 3: distribution of lipid carbon content (ug) / species / stage

boxplot(LC[Ch]~SG[Ch], ylim=c(0,2500), col=c(grey,blue,green,red), main="Lipid Carbon (ug)")
par(new=TRUE)
boxplot(LC[Cg]~SG[Cg], ylim=c(0,2500), col=c(grey,blue,green,red), main=" ")
par(new=TRUE)
boxplot(LC[Ml]~SG[Ml], ylim=c(0,2500), col=c(grey,blue,green,red), main=" ")
legend("topleft", legend=c("CF","C5","C4","C3"), pch=15, col=c(red,green,blue,grey), pt.cex=3, cex=1.5, bty="n")

#****************************************************************************************


#--- Conversion from LC (ug) to structural carbon SC (ug)
#    We can vary alpha to change the carbon content of the structural part
#    of the copepod's body relative to the WE carbon from the oil sac.

alpha <- 1 #0.5

TC    <- LC / LF

SC    <- ( TC - LC ) * alpha # conservaive estimate; likely to OVERestimate the structural carbon content

#*** Alternate more complex scenario:
# =>   use TL & DW relationship from Vogedes et al. (2010) doi:10.1093/plankt/fbq068 
# => + use species-specific DW & C relationships from Forest et al. 2010 doi:10.1093/plankt/fbq103

#DW     <- 1e3 * ( TL - 0.037 ) / 0.453

#TC     <- NA*DW
#TC[Ch] <- 0.634 * DW[Ch] + 141.65
#TC[Cg] <- 0.585 * DW[Cg] + 15.09

#SC <- TC * (1-LF)

# => Does not change the results qualitatively !


#****************************************************************************************
# Figure 4: structural carbon frequency distribution / stage

boxplot(SC[Ch]~SG[Ch], ylim=c(0,1500), col=c(grey,blue,green,red), main="Structural carbon (ug)")
par(new=TRUE)
boxplot(SC[Cg]~SG[Cg], ylim=c(0,1500), col=c(grey,blue,green,red), main=" ")
par(new=TRUE)
boxplot(SC[Ml]~SG[Ml], ylim=c(0,1500), col=c(grey,blue,green,red), main=" ")
legend("topleft", legend=c("CF","C5","C4","C3"), pch=15, col=c(red,green,blue,grey), pt.cex=3, cex=1.5, bty="n")

#****************************************************************************************


#--- Potential diapause duration (& overhead for capital breeding) from Maps et al. 2014 (doi:10.1093/plankt/fbt100)
#    Energetic parameters from mixed-effects linear model of regression of multi-species literature review.
#    Data for M. longa from Hirche 1987 (doi:10.1007/BF00428240)

m0 <- vector(length = length(PL))
E  <- m0

# Metabolic constant at 0 degC (ug C^-1/4 s^-1)
m0[Ch] <- 2.56e-7
m0[Cg] <- 2.25e-7
m0[Ml] <- 6.34e-7

# Activation energy of metabolism (eV)
E[Ch] <- 0.6854
E[Cg] <- 0.6038
E[Ml] <- 0.4791

#k  <- 8.6170e-05 # Boltzmann constant (eV)
#T0 <- 273.15     # reference temperature (K)

# Computes duration (days) until ALL lipids are depleted at 0 degC
D <- 4 / (-m0*86400) * ( SC^0.25 - TC^0.25 ) #* exp(-E*0/(k*(0+T0)*T0)) => temperature term not required since ambient temperature is ~ T0


#****************************************************************************************
# Figure 5: diapause duration frequency distribution / species / stage

z <- D[Ch|Cg]
x <- factor(SP[Ch|Cg])
y <- factor(SG[Ch|Cg])

par(mar=c(4,5,3,3))

boxplot(z~x*y, xlim=c(1,8.5), ylim=c(0,600), xlab=" ", ylab="days", main=" ", col=c(blue,red), names=c(" ","C3","C4","C4","C5","C5","F","F"), main=" ", lwd=1.5, cex.lab=2, cex.axis=1.5, frame.plot=F)

abline(a=182, b=0, lty=3, col="black", lwd=2)
abline(a=365, b=0, lty=2, col="black", lwd=2)

legend("bottomright", legend=c( expression(italic("C. glacialis")), expression(italic("C. hyperboreus")) ), pch=c(15,15), pt.cex=2, cex=1.2, col=c(blue,red), bty="n")

# Optional diagnostic plot
if (Pflag) {
  print( "Which species do you want ? C. hyperboreus: Ch | C. glacialis: Cg | M. longa: Ml " ); Sflag <- scan(,what="character")
  if (Sflag=="Ch") {
    Psp   <- Ch
    Pmain <- "C. hyperboreus"
  } else if (Sflag=="Cg") {
    Psp   <- Cg
    Pmain <- "C. glacialis"
  } else if (Sflag=="Ml") {
    Psp   <- Ml
    Pmain <- "M. longa"
  }
  boxplot(D[Psp]~SG[Psp], ylim=c(0,600), col=c(grey,blue,green,red), ylab="Time to deplete TL (d)", main=Pmain)
  abline(a=30,  b=0, lty=3, col="black", lwd=2)
  abline(a=365, b=0, lty=2, col="black", lwd=2)
}

#****************************************************************************************


#****************************************************************************************
# Figure 6: diapause duration vs depth / species / stage

# C. hyperboreus
plot(   D[Ch&C5&Hd], -Z[Ch&C5&Hd], pch=1, col="black", xlim=c(0,800), xlab="days", ylim=c(-350,0), ylab="depth (m)", main=expression( paste( "Time to deplete TL - ",italic("C. hyperboreus") ) ) )
points( D[Ch&C5&Hn], -Z[Ch&C5&Hn], pch=6,  col="blue")
points( D[Ch&CF],    -Z[Ch&CF],    pch=15, col="red")
legend("bottomleft", c("C5 day","C5 night","females"), pch=c(1,6,15), col=c("black","blue","red"), bty="n")

# C. glacialis
plot(   D[Cg&C5&Hd], -Z[Cg&C5&Hd], pch=1, col="black", xlim=c(0,500), xlab="days", ylim=c(-350,0), ylab="depth (m)", main=expression( paste( "Time to deplete TL - ",italic("C. glacialis") ) ) )
points( D[Cg&C5&Hn], -Z[Cg&C5&Hn], pch=6, col="blue")
points( D[Cg&CF],    -Z[Cg&CF],    pch=15, col="red")
legend("bottomleft", c("C5 day","C5 night","females"), pch=c(1,6,15), col=c("black","blue","red"), bty="n")

# M. longa
plot(   D[Ml&C5&Hd], -Z[Ml&C5&Hd], pch=1, col="black", xlim=c(0,100), xlab="days", ylim=c(-350,0), ylab="depth (m)", main=expression( paste( "Time to deplete TL - ",italic("M. longa") ) ) )
points( D[Ml&C5&Hn], -Z[Ml&C5&Hn], pch=6, col="blue")
points( D[Ml&CF],    -Z[Ml&CF],    pch=15, col="red")
legend("bottomright", c("C5 day","C5 night","females"), pch=c(1,6,15), col=c("black","blue","red"), bty="n")

#****************************************************************************************
