#Nicole E Soltis
#04/06/18

#----------------------------------------------------------------------------
rm(list=ls())
setwd("~/Projects/BcAtGWAS/data")
LsEcc <- read.csv("REF_BcLesionTraits/LesionEccGWA.csv")
LsGrn <- read.csv("REF_BcLesionTraits/LesionGreenGWA.csv")
LsYlw <- read.csv("REF_BcLesionTraits/LesionYellowGWA.csv")
LsArea <- read.csv("REF_BcLesionTraits/LesionAreaGWA.csv")

library(ggplot2); 

#Ylw and Ecc: only plot A517
#plot like fx permut in supp table
myDat <- LsEcc 
#myDat <- LsYlw #done

#Grn and Area: plot all isolates but as individual manhattan plots
#myDat <- LsGrn #done
#myDat <- LsArea #done

#create a custom color scale
#myColors <- c("grey20", "grey60", "grey20", "grey60", "grey20")
myColors <- c("royalblue3", "blue4", "royalblue3", "blue4", "royalblue3")
names(myColors) <- levels(myDat$chr)
colScale <- scale_colour_manual(name = "chr",values = myColors)

#sort dataframe rows in order of Chrom, then Pos
str(myDat)
myDat$chrom <- as.numeric(myDat$chr)
myDat <- myDat[with(myDat, order(chrom, pos)), ]

#Make plotting variables
myDat$Index = NA
lastbase = 0

#Redo the positions to make them sequential		-- accurate position indexing
for (i in unique(myDat$chrom)) {
  print(i)
  if (i==1) {
    #for the subset of HEM.plotdata rows with Chromosome 1, set Index variable for each row to equal Pos.
    myDat[myDat$chrom==i, ]$Index=myDat[myDat$chrom==i, ]$pos
    #for all other chromosomes: 
  }	else {
    #lastbase for chromosome i is the greater of:
    #current lastbase counter plus the maxiumum position of chromosome i-1
    #OR 1
    lastbase=lastbase+max(subset(myDat,myDat$chrom==i-1)$pos, 1)
    #and then for the subset of HEM.plotdata rows with Chromosome i, set Index variable for each row to equal Pos + lastbase
    myDat[myDat$chrom==i, ]$Index=myDat[myDat$chrom==i, ]$pos+lastbase
  }
}

hist(myDat$pos)
hist(myDat$Index)
#positions look fine...

#thresholds
my95Thr <- myDat[c(myDat$Thresholds==0.95), c(9:13)]
#my95Thr <- myDat[c(myDat$Thresholds==95), c(9:13)] #for Lesion Yellow only

#Grn and Area: plot all isolates but as individual manhattan plots
#Ylw and Ecc: only plot A517
names(myDat)

#get Thr
#recode this depending on genotype, range = 4:7
for (myz in (4)){
mytrait <- myz 
varname <- names(myDat)[mytrait]
my95val <- my95Thr[1,paste(varname, "1", sep=".")]

setwd("~/Projects/BcAtGWAS")
jpeg(paste("plots/manhattan/",varname,".jpg", sep=""), width=8, height=3, units='in', res=600)
#y can substitute out
print(
  ggplot(myDat, aes(x=Index, y=get(varname)))+
        theme_bw()+
        colScale+
        geom_point(aes(color = factor(chr),alpha=0.001), size=1)+
        labs(list(y="SNP effect estimate", title="Lesion Eccentricity Apple517"))+
                    #paste("Lesion Yellow",varname, sep=" ")))+
        guides(col = guide_legend(nrow = 8, title="Chromosome"))+
        geom_hline(yintercept=my95val, colour = "red", lty=1) +
        geom_hline(yintercept=-my95val, colour = "red", lty=1)+
        theme(legend.position="none", plot.title = element_text(hjust = 0.5))+
        scale_x_continuous(name="Chromosome", breaks = c(15215061, 40272785, 61855075, 82864894, 105641530), labels = c("1", "2", "3", "4", "5"))+
              expand_limits(y=0)
)
      dev.off()
}     
#---------------------------------------------------------------------------------
      #get chromosome midpoints
      my.chroms <- as.data.frame(myDat[!duplicated(myDat$chrom, fromLast=FALSE), "Index"]) #Lower Bounds
      names(my.chroms)[1] <- "Chr.Start"
      my.chroms$Chr.End <- myDat[!duplicated(myDat$chrom, fromLast=TRUE), "Index"] # Upper Bounds
      my.chroms$Chr.Mid <- (my.chroms$Chr.Start + my.chroms$Chr.End)/2
my.chroms$Chr.Mid
