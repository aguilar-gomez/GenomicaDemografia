# Snigenda Workshop.
# author: Kaden Winspear 
# Purpose: Plot GONE output for the last 50 generations. Produces one pdf
# with both plots. 

# Set working directory
setwd("/Users/snigenda/Documents/workshop_2025/results/gone_ne")

# open the pdf
pdf("GONE2_GOC_ENP_last50.pdf")
par(mfrow = c(2, 1))

# GOC PLOT
goc <- read.table("GOC_GONE2_Ne", header = TRUE)

# get most recent 50 generations (first 50 rows)
goc_50 <- goc[1:50, ]

plot(goc_50$Generation, goc_50$Ne_diploids,
     type = "l", lwd = 2, ylim = c(0, 600), xaxp = c(0, 50, 10),
     xlab = "Generations Before Present",
     ylab = "Effective Population Size (Ne)")


# ENP PLOT
enp <- read.table("ENP_GONE2_Ne", header = TRUE)

# get most recent 50 generations (first 50 rows)
enp_50 <- enp[1:50, ]

plot(enp_50$Generation, enp_50$Ne_diploids,
     type = "l", lwd = 2, ylim = c(0, 25000), xaxp = c(0, 50, 10),
     xlab = "Generations Before Present",
     ylab = "Effective Population Size (Ne)")
# Close pdf
dev.off()
