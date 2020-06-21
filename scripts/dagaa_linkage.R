######################################
##### Dagaa SNP Linkage Analysis #####
######################################
## J. Rick, June 2020
## Run for each species

# For both sets of loci, we're interested in whether they are more linked than expected by chance
# this would suggest that we are correct in calling them "sex" and "inversion" loci.

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

###############################
# Data import and cleaning ####
###############################

spp <- "limno" # specify "stolo" or "limno"

## import linkage data from PLINK ##
dagaa_dt <- fread(paste("../data/",spp,"_all_120219.ld",sep=""))

setindex(dagaa_dt, CHR_A, BP_A, CHR_B, BP_B)
indices(dagaa_dt)

## import loci identified from DAPC analyses ##
inv.loci <- read.csv("../data/limno_inversion_loci.csv",header=TRUE)
sex.loci <- read.csv(paste("../data/",spp,"_sex_loci.csv",sep=""),header=TRUE)

inv.chr <- unique(inv.loci$scaffold)
sex.chr <- unique(sex.loci$scaffold)

## plot all linkage values ##
par(mfrow=c(1,2))
#hist(dagaa_dt[,DP])
hist(dagaa_dt[,R2])
#plot(limno_ld$R2,limno_ld$DP)

## assign loci as "sex","inv", or neither ##
if (spp == "limno") {
  dagaa_dt[, sig := c("none", "inv", "sex")[1 +
                             1 * (dagaa_dt$CHR_A %in% inv.loci$scaffold & 
                                    dagaa_dt$BP_A %in% inv.loci$pos & 
                                    dagaa_dt$CHR_B %in% inv.loci$scaffold & 
                                    dagaa_dt$BP_B %in% inv.loci$pos) + 
                             2 * (dagaa_dt$CHR_A %in% sex.loci$scaffold & 
                                    dagaa_dt$BP_A %in% sex.loci$pos & 
                                    dagaa_dt$CHR_B %in% sex.loci$scaffold & 
                                    dagaa_dt$BP_B %in% sex.loci$pos)]
  ]
} elseif (spp == "stolo") {
  dagaa_dt[, sig := c("none", "sex")[1 +
                             1 * (dagaa_dt$CHR_A %in% sex.loci$scaffold & 
                                    dagaa_dt$BP_A %in% sex.loci$pos & 
                                    dagaa_dt$CHR_B %in% sex.loci$scaffold & 
                                    dagaa_dt$BP_B %in% sex.loci$pos)]
  ]
}

dagaa_dt <- dagaa_dt[, sig:=as.factor(sig)]

# ns_dt <- dagaa_dt[sig == "none"]
# sex_dt <- dagaa_dt[sig == "sex"]
# inv_dt <- dagaa_dt[sig == "inv"]

sex.test <- wilcox.test(R2 ~ sig, data=dagaa_dt, subset = sig %in% c("none","sex"))

if (spp == "limno") {
  inv.test <- wilcox.test(R2 ~ sig, data=dagaa_dt, subset = sig %in% c("none","inv"))
  sex.inv.test <- wilcox.test(R2 ~ sig, data=dagaa_dt, subset = sig %in% c("sex","inv"))
}


## density plots ##
## WARNING: very large objects ##

dens1 <- ggdensity(dagaa_dt[sig != "sex"], x="R2", add="mean", color="sig", fill="sig",
          palette = c("gray20", #"gray40", 
                      "gray60", "gray80")) +
  scale_x_continuous(limits=c(0.2,1))
#ggexport(dens1,filename="../results/limno_LD_density_sex_only_R2.png",height=500,width=1200)

dens2 <- ggdensity(dagaa_dt[sig != "sex"], x="R2", add="mean", color="sig", fill="sig",
          palette = c("gray20", #"gray40", 
                      "gray60", "gray80"),
          xlim=c(0.2,1))
#ggexport(dens2,filename="../results/limno_LD_density_sex_only_min20_031120.png",height=500,width=1200)

pdf(file=paste("../results/",spp,"_LD_density.pdf",spp=""))      
#print(dens1)
#print(dens2)
print(ggdensity(dagaa_dt[sig != "sex"], x="R2", add="mean", color="sig", fill="sig",
          palette = c("gray20", #"gray40", 
                      "gray60", "gray80"),
          xlim=c(0.2,1)))
#Rmisc::multiplot(dens1,dens2,cols=2)
dev.off()


##########
