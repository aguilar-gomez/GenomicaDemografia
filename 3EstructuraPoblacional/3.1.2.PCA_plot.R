setwd("~/Documents/GitHub/GenomicaDemografia/3EstructuraPoblacional")
library(ggplot2)
library(ggrepel)
library(gplots)

#read files
eigvectors <- read.table("../data/pca/finwhale_pca.eigenvec")
eigenval <- read.table("../data/pca/finwhale_pca.eigenval")


vecs<-eigvectors[c(3:22)]

#Name columns
colnames(eigvectors)<-c("sample","names",paste0("PC",c(1:10)))

#Population name is the first three letter
eigvectors$group <- substr(eigvectors$names, 1, 3)


# Eigenvalues: how much variance each component explains
eigvalue <- eigenval$V1/sum(eigenval$V1);
cat(signif(eigvalue , digits=3)*100,"\n");

# Plot
comp<-c(1,2)
title <- paste("PC",comp[1]," (",signif(eigvalue[comp[1]], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eigvalue[comp[2]], digits=3)*100,"%)",sep="",collapse="")
xlabel = paste("PC",comp[1]," (",signif(eigvalue[comp[1]], digits=3)*100,"%)",sep="")
ylabel = paste("PC",comp[2]," (",signif(eigvalue[comp[2]], digits=3)*100,"%)",sep="")
x_axis = paste("PC",comp[1],sep="")
y_axis = paste("PC",comp[2],sep="")

#colores
cbPalette <- c( "#58BBDC","#e9c856")


ggplot() + geom_point(data=eigvectors, aes_string(x=x_axis, y=y_axis, color="group"),alpha = .8,size=3)+
  scale_colour_manual(values=cbPalette)+ theme_bw(base_size = 16) +xlab(xlabel)+ylab(ylabel)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(filename=paste0("PCA_finwhale","PC",comp[1],"PC",comp[2],".png"),width=6,height=4)


#Añadir etiquetas
ggplot(eigvectors, aes(x = .data[[x_axis]],
                       y = .data[[y_axis]],
                       color = group)) +
  geom_point(alpha = 0.8, size = 3) +
  geom_label_repel(aes(label = sample), size = 3, max.overlaps = Inf) +
  scale_colour_manual(values = cbPalette) +
  theme_bw(base_size = 16) +
  xlab(xlabel) +
  ylab(ylabel) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(filename = paste0("PCA_finwhale_labels_PC", comp[1], "_PC", comp[2], ".png"),
  width = 6, height = 4, dpi = 300)

