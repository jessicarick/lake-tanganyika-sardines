###############################################################
## Code for calculating FSTs and conducting Mantel tests
## for Limnothrissa and Stolothrissa
##
## J. Rick, 8/3/19
## Updated 6/20/20
###############################################################

source("packages_functions.R")

#########################################################
## prepping data in the same way as in PCA/DAPC analyses
## skip this section if "stologen.both" and "limnogen.both" already exist
#########################################################

## LIMNOTHRISSA
limno_vcfR<-read.vcfR("../../data/limno_all_0.5_maf0.01.recode.vcf")

limnogen <- vcfR2genlight(limno_vcfR)
col.names <- unlist(strsplit(indNames(limnogen),"/project/wagnerlab/jrick/dagaa/rad/bamfiles/aln_"))
col.names.clean <- unlist(strsplit(col.names,".sorted.bam"))
indNames(limnogen)<-col.names.clean

fishinfo<-read.csv('../../data/sardine_all_metadata_clean.csv',
                   header=TRUE,stringsAsFactors = FALSE)
pairedinfolimno <- left_join(data.frame(sample_name=col.names.clean),fishinfo,by="sample_name")
head(pairedinfolimno)

## filtering SNPs found only in one group or the other
pop(limnogen) <- pairedinfolimno$sequencing_center

floragenex <- gl.keep.pop(limnogen,
                          "Floragenex")
floragenex <- gl.recalc.metrics(floragenex)
floragenex.highcalls <- gl.filter.callrate(floragenex, method="loc", threshold=0.95)

lausanne <- gl.keep.pop(limnogen,
                        "Lausanne")
lausanne <- gl.recalc.metrics(lausanne)
lausanne.highcalls <- gl.filter.callrate(lausanne, method="loc", threshold=0.95)

snps.both <- intersect(locNames(floragenex.highcalls),locNames(lausanne.highcalls))
limnogen.both <- limnogen[,locNames(limnogen) %in% snps.both]

## STOLOTHRISSA
stolo_vcfR<-read.vcfR("../../data/stolo_all_0.5_maf0.01.recode.vcf")

stologen <- vcfR2genlight(stolo_vcfR)
col.names <- unlist(strsplit(indNames(stologen),"/project/wagnerlab/jrick/dagaa/rad/bamfiles/stolo/aln_"))
col.names.clean <- unlist(strsplit(col.names,".sorted.bam"))
indNames(stologen)<-col.names.clean

fishinfo<-read.csv('../../data/sardine_all_metadata_clean.csv',
                     header=TRUE,stringsAsFactors = FALSE)
pairedinfostolo <- left_join(data.frame(sample_name=col.names.clean),fishinfo,by="sample_name")

## filtering SNPs found only in one group or the other
pop(stologen) <- pairedinfostolo$sequencing_center

floragenex <- gl.keep.pop(stologen,
                          "Floragenex")
floragenex <- gl.recalc.metrics(floragenex)
floragenex.highcalls <- gl.filter.callrate(floragenex, method="loc", threshold=0.95)

lausanne <- gl.keep.pop(stologen,
                        "Lausanne")
lausanne <- gl.recalc.metrics(lausanne)
lausanne.highcalls <- gl.filter.callrate(lausanne, method="loc", threshold=0.95)

snps.both <- intersect(locNames(floragenex.highcalls),locNames(lausanne.highcalls))
stologen.both <- stologen[,locNames(stologen) %in% snps.both]

#######################################
#######################################
## First, limnothrissa
#######################################
#######################################
pcaAll.limno <- read.csv("pcaAll.limno.csv",header=TRUE) # from PCA script
Dgeo <- read.csv("../../data/dist_matrix.csv",header=TRUE,row.names = 1)

limnogen_nolowcov <- limnogen.both[which(indNames(limnogen.both) %in% 
                                           pcaAll.limno$sample.id),]

pcaAll.limno$site <- factor(pcaAll.limno$site)
limnogen_nolowcov@pop <- pcaAll.limno$site

# calculate reich-patterson fst
limno.reich.fst <- reich.fst(limnogen_nolowcov, plot=TRUE, verbose=TRUE, bootstrap=100)

# subset geographic distance
Dgeo_limno <- as.dist(Dgeo[row.names(Dgeo) %in% unique(limnogen_nolowcov@pop),
                           colnames(Dgeo) %in% unique(limnogen_nolowcov@pop)])

# combine together
limno_fst_dgeo <- melt(limno.reich.fst$fsts,value.name="dist.fst") %>%
  full_join(melt(as.matrix(Dgeo_limno),value.name="dist.km"),by=c("Var1","Var2")) 

limno_fst_dgeo$fst_1fst <- limno_fst_dgeo$dist.fst / (1 - limno_fst_dgeo$dist.fst)

# plotting all sampling sites
plot(limno_fst_dgeo$dist.km,
     limno_fst_dgeo$fst_1fst,
     xlab="geographic distance (km)",ylab="fst / 1-fst")

# now, omitting Kivu from the plots and linear model
plot(limno_fst_dgeo$dist.km[limno_fst_dgeo$dist.km < 1000],
     limno_fst_dgeo$fst_1fst[limno_fst_dgeo$dist.km < 1000],
     xlab="geographic distance (km)",ylab="fst / 1-fst")
abline(lm(limno_fst_dgeo$fst_1fst[limno_fst_dgeo$dist.km < 1000] ~
            limno_fst_dgeo$dist.km[limno_fst_dgeo$dist.km < 1000]),lty=2)
summary(lm(limno_fst_dgeo$fst_1fst[limno_fst_dgeo$dist.km < 1000] ~
             limno_fst_dgeo$dist.km[limno_fst_dgeo$dist.km < 1000])) ## not significant

## testing for IBD -- mantel test ##
Dgen_limno_fst_1fst <- as.dist(spread(limno_fst_dgeo[,c(1,2,5)], Var2, fst_1fst)[,-c(1)])
Dgeo_limno_km <- as.dist(spread(limno_fst_dgeo[,c(1,2,4)], Var2, dist.km)[,-c(1)])
ibd_limno <- mantel.randtest(Dgen_limno_fst_1fst,Dgeo_limno_km)
ibd_limno # not significant

# IBD plot
plot(Dgeo_limno_km,Dgen_limno_fst_1fst,
     xlab="Geographic Distance (km)",
     ylab="Genetic Distance (1/(1-Fst))",
     cex.lab=1.5, pch=19, col=scales::alpha("gray20",0.5),
     cex=2,cex.axis=1.2)
abline(lm(Dgen_limno_fst_1fst ~ Dgeo_limno_km),
       lty = 2, lwd=2, col="gray20")
lm.sum <- summary(lm(Dgen_limno_fst_1fst ~ Dgeo_limno_km))
text(x=0,y=0.013,paste("R2 = ",round(lm.sum$adj.r.squared,3),"\np = ",round(lm.sum$coefficients[2,4],3)),pos=4)

## now, IBD removing small pops ##
smallpops <- c("Kivu", "Mbete", "Kabimba")
largepops <- levels(limnogen_nolowcov@pop)[!(levels(limnogen_nolowcov@pop) %in% smallpops)]

limno_fst_dgeo_nosmall <- limno_fst_dgeo[limno_fst_dgeo$Var1 %in% largepops &
                                           limno_fst_dgeo$Var2 %in% largepops,]
Dgen_limno_fst_1fst_nosmall <- as.dist(spread(limno_fst_dgeo_nosmall[,c(1,2,5)], Var2, fst_1fst)[,-c(1)])
Dgeo_limno_km_nosmall <- as.dist(spread(limno_fst_dgeo_nosmall[,c(1,2,4)], Var2, dist.km)[,-c(1)])
ibd_limno_nosmall <- mantel.randtest(Dgen_limno_fst_1fst_nosmall,Dgeo_limno_km_nosmall)
ibd_limno_nosmall # not significant

plot(ibd_limno_nosmall)
plot(Dgeo_limno_km_nosmall,Dgen_limno_fst_1fst_nosmall)
abline(lm(Dgen_limno_fst_1fst_nosmall~Dgeo_limno_km_nosmall))

plot(Dgeo_limno_km_nosmall,Dgen_limno_fst_1fst_nosmall,
     xlab="Geographic Distance (km)",
     ylab="Genetic Distance (Fst / 1-Fst)",
     cex.lab=1.5, pch=19, col=scales::alpha("gray20",0.5),
     cex=2,cex.axis=1.2,
     xlim=c(0,600),ylim=c(0,0.005))
abline(lm(Dgen_limno_fst_1fst_nosmall~Dgeo_limno_km_nosmall),lty=2,lwd=2)
lm.sum <- summary(lm(Dgen_limno_fst_1fst_nosmall~Dgeo_limno_km_nosmall))

text(x=450,y=0.0045,paste("R2 = ",round(lm.sum$adj.r.squared,3),"\np = ",round(lm.sum$coefficients[2,4],3)),pos=4)


#######################################
## Now, stolothrissa
#######################################
pcaAll.stolo <- read.csv("pcaAll.stolo.csv",header=T) # from PCA script
pcaAll.stolo$site <- factor(pcaAll.stolo$site, levels=c("Kilomoni","Lusenda","Kagunga","Kigoma","Kabimba",
                                                        "NorthMahale","SouthMahale","Ikola",
                                                        "Kipili","Mbete"))

stologen_nolowcov <- stologen.both[which(indNames(stologen.both) %in% 
                                           pcaAll.stolo$sample.id),]
stologen_nolowcov@pop <- pcaAll.stolo$site

# calculate reich-patterson fst
stolo.reich.fst <- reich.fst(stologen_nolowcov, plot=TRUE, verbose=TRUE, bootstrap=100)

# subset geographic distance
Dgeo_stolo <- as.dist(Dgeo[row.names(Dgeo) %in% unique(stologen_nolowcov@pop),
                           colnames(Dgeo) %in% unique(stologen_nolowcov@pop)])

# combine together
stolo_fst_dgeo <- melt(stolo.reich.fst$fsts,value.name="dist.fst") %>%
  full_join(melt(as.matrix(Dgeo_stolo),value.name="dist.km"),by=c("Var1","Var2")) 

stolo_fst_dgeo$fst_1fst <- stolo_fst_dgeo$dist.fst / (1 - stolo_fst_dgeo$dist.fst)

plot(stolo_fst_dgeo$dist.km,
     stolo_fst_dgeo$fst_1fst,
     xlab="geographic distance (km)",ylab="fst / 1-fst")
abline(lm(stolo_fst_dgeo$fst_1fst ~
            stolo_fst_dgeo$dist.km),lty=2)
summary(lm(stolo_fst_dgeo$fst_1fst ~
             stolo_fst_dgeo$dist.km))

## testing for IBD -- mantel test ##
Dgen_stolo_fst_1fst <- as.dist(spread(stolo_fst_dgeo[,c(1,2,5)], Var2, fst_1fst)[,-c(1)])
Dgeo_stolo_km <- as.dist(spread(stolo_fst_dgeo[,c(1,2,4)], Var2, dist.km)[,-c(1)])
ibd_stolo <- mantel.randtest(Dgen_stolo_fst_1fst,Dgeo_stolo_km)
ibd_stolo # not significant


plot(fst_1fst ~ dist.km, data=stolo_fst_dgeo[stolo_fst_dgeo$dist.fst > 0,],
     xlab="Geographic Distance (km)",
     ylab="Genetic Distance (1/(1-Fst))",
     cex.lab=1.5, pch=19, col=scales::alpha("gray20",0.5),
     cex=2,cex.axis=1.2)
abline(lm(fst_1fst ~ dist.km, data=stolo_fst_dgeo[stolo_fst_dgeo$dist.fst > 0,]),
       lty = 2, lwd=2, col="gray20")
lm.sum <- summary(lm(fst_1fst ~ dist.km, data=stolo_fst_dgeo[stolo_fst_dgeo$dist.fst > 0,]))
text(x=450,y=0.0045,paste("R2 = ",round(lm.sum$adj.r.squared,3),"\np = ",round(lm.sum$coefficients[2,4],3)),pos=4)

## now, removing small pops ##
smallpops <- c("Kilomoni","Kabimba","SouthMahale")
largepops <- levels(stologen_nolowcov@pop)[!(levels(stologen_nolowcov@pop) %in% smallpops)]

stolo_fst_dgeo_nosmall <- stolo_fst_dgeo[stolo_fst_dgeo$Var1 %in% largepops &
                                           stolo_fst_dgeo$Var2 %in% largepops,]
Dgen_stolo_fst_1fst_nosmall <- as.dist(spread(stolo_fst_dgeo_nosmall[,c(1,2,5)], Var2, fst_1fst)[,-c(1)]) # Edward's distance
Dgen_stolo_fst_1fst_nosmall[Dgen_stolo_fst_1fst_nosmall < 0] <- 0
Dgeo_stolo_km_nosmall <- as.dist(spread(stolo_fst_dgeo_nosmall[,c(1,2,4)], Var2, dist.km)[,-c(1)])
ibd_stolo_nosmall <- mantel.randtest(Dgen_stolo_fst_1fst_nosmall,Dgeo_stolo_km_nosmall)
ibd_stolo_nosmall # not significant; p=0.354

plot(ibd_stolo_nosmall)
plot(Dgeo_stolo_km_nosmall,Dgen_stolo_fst_1fst_nosmall)
abline(lm(Dgen_stolo_fst_1fst_nosmall~Dgeo_stolo_km_nosmall))

plot(fst_1fst ~ dist.km, data=stolo_fst_dgeo_nosmall[stolo_fst_dgeo_nosmall$dist.fst > 0,],
     xlab="Geographic Distance (km)",
     ylab="Genetic Distance (Fst / 1-Fst)",
     cex.lab=1.5, pch=19, col=scales::alpha("gray20",0.5),
     cex=2,cex.axis=1.2)
abline(lm(fst_1fst ~ dist.km, data=stolo_fst_dgeo_nosmall[stolo_fst_dgeo_nosmall$dist.fst > 0,]),lty=2,lwd=2)
lm.sum <- summary(lm(fst_1fst ~ dist.km, data=stolo_fst_dgeo_nosmall[stolo_fst_dgeo_nosmall$dist.fst > 0,]))

text(x=450,y=0.002,paste("R2 = ",round(lm.sum$adj.r.squared,3),"\np = ",round(lm.sum$coefficients[2,4],3)),pos=4)



#################################################
#### Plotting limno and stolo results together

par(mfrow=c(1,2),mar=c(4,5,1,1),oma=c(1,1,1,1))

plot(fst_1fst ~ dist.km, data=stolo_fst_dgeo_nosmall[stolo_fst_dgeo_nosmall$dist.km > 0,],
     xlab="Geographic Distance (km)",
     ylab="Genetic Distance (Fst / 1-Fst)",
     cex.lab=1.5, pch=19, col=scales::alpha("gray20",0.5),
     cex=2,cex.axis=1.2)
abline(lm(fst_1fst ~ dist.km, data=stolo_fst_dgeo_nosmall[stolo_fst_dgeo_nosmall$dist.km > 0,]),lty=2,lwd=2,col="black")

plot(fst_1fst ~ dist.km, data=limno_fst_dgeo_nosmall[limno_fst_dgeo_nosmall$dist.km > 0,],
     xlab="Geographic Distance (km)",
     ylab="Genetic Distance (Fst / 1-Fst)",
     cex.lab=1.5, pch=19, col=scales::alpha("gray20",0.5),
     cex=2,cex.axis=1.2)
abline(lm(fst_1fst ~ dist.km, data=limno_fst_dgeo_nosmall[limno_fst_dgeo_nosmall$dist.km > 0,]),lty=2,lwd=2,col="black")

