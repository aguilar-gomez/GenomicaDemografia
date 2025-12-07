setwd("~/Documents/GitHub/GenomicaDemografia/3EstructuraPoblacional")

library(ggplot2)

#read in input file
data <- read.table("~/Documents/GitHub/GenomicaDemografia/data/admixture/finwhale.fam")
sampling<-data[,1]


#sort and read in log files
runs <- list()
cbPalette <- list(
  `1` = c("#58BBDC"),
  `2` = c("#e9c856","#58BBDC"),
  `3` = c("#e9c856","#CA6E28", "#58BBDC"))

for (i in 1:3){
  tempdf <- read.table(paste0("~/Documents/GitHub/GenomicaDemografia/data/admixture/","finwhale.", i, ".Q"))
  temp <- as.data.frame(tempdf)
  
  # Add sample names
  temp$names <- sampling
  
  # Store sorted data frame in list
  runs[[i]] <- temp
}

#plot runs 1:3
png("AdmixtureK1-3_Finwhale.png",
    width = 1000, height = 1200, res = 300)  # 8x4 pulgadas aprox.
par(mfrow=c(4,1),
    mar=c(0, 3.7, 0.5, 0.5),
    mgp=c(1, 0.4, 0),
    xaxs="i")


for (i in 1:3){
  colors_for_K <- cbPalette[[as.character(i)]]
  if(i==3){
    barplot(t(as.matrix(runs[[i]])), col=colors_for_K, border="black",
            names.arg=sampling, las=2, cex.names=1,
            cex.axis=0.8, yaxt="n")
  } else {
    barplot(t(as.matrix(runs[[i]])), col=colors_for_K, border="black",
            names.arg=rep("", nrow(runs[[i]])),
            cex.axis=0.8, yaxt="n")
  }
  mtext(paste("K =", i), side=2, line=0.5, cex=0.9)
}
dev.off()


#read in log error values to determine optimal K
log<-read.table("~/Documents/GitHub/GenomicaDemografia/data/admixture/log.errors.txt")[,c(3:4)]
#use double backslash to interpret the opening parentheses literally in the regular expression
log$V3<-gsub("\\(K=", "", log$V3)
log$V3<-gsub("):", "", log$V3)
#interpret K values as numerical
log$V3<-as.numeric(log$V3)
#rename columns
colnames(log)<-c("Kvalue","cross.validation.error")

#make plot showing the cross validation error across K values 1:10
#lowest CV value is the ideal K value
ggplot(data=log, aes(x=Kvalue, y=cross.validation.error, group=1)) +
  geom_line(linetype = "dashed")+
  geom_point()+
  ylab("cross-validation error")+
  xlab("K")+
  scale_x_continuous(breaks = c(1:10))+
  theme_classic(base_size = 18)
ggsave("CrossValidationV.png",width = 4, height=4)

