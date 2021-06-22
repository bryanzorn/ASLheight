library("imager")
library("ggplot2")
library("plyr")
library("randomForest")
library("gridExtra")

# Extract file names and experiment type, build dataframe for recording results
files <- list.files(pattern = "\\.jpg$", recursive=TRUE)
files <- files[!grepl("Logo", files)]

df <- as.data.frame(c())

for (file in files) {
	fname <- as.data.frame(
			t(
			rev(
			unlist(
			strsplit(file, "/")
			))))
	df <- rbind(df, fname)
}

df <- df[complete.cases(df), ]
colnames(df) <- c("jpg","markfind","timepoint")

df2 <- as.data.frame(c())

for (lvl1 in levels(df$markfind)){
	for (lvl2 in levels(df$timepoint)) {
		tmpdf <- subset(df, markfind == lvl1)
		tmpdf <- subset(tmpdf, timepoint == lvl2)
		lendf <- seq(1:nrow(tmpdf))
		if (dim(tmpdf)[1] != 0) {
			df2 <- rbind(df2,cbind(tmpdf, lendf))
		}
	}
}

cutoff <- 486

Size <- c()
histdata <- as.data.frame(c())

# Load and deconstruct images, find bands, define peaks and draw boundaries for review
pdf("ReviewImages.pdf", width=9, height=3)
exclusionvals <- c()
for (file in files) {
		img <- load.image(file) %>% imsub(y < cutoff)
		red <- as.matrix(R(img))
		grn <- as.matrix(G(img))
		blue <- as.matrix(B(img))
		compr <- imrotate(as.cimg(red + grn + blue),90)

		exclusionvals <- append(exclusionvals, mean((compr > mean(compr)+sd(compr)+sd(compr))))

		comprdf <- as.data.frame(compr > mean(compr)+sd(compr)+sd(compr))
		den <- density(comprdf$x)
		val <- as.data.frame(cbind(den$x,den$y))
		pred <- predict(lm(val$V2 ~ poly(val$V1, 16, raw=TRUE)))
		newframe <- cbind(seq(1:length(pred)), as.data.frame(pred))
		colnames(newframe) <- c("x","y")
	
		valley <- as.data.frame(c())
		peaks <- as.data.frame(c())
		
		for (value in seq(1:nrow(newframe))) {
			x <- newframe$y[value]
			y <- newframe$y[value+1]
			z <- newframe$y[value+2]
			if (is.na(x < y && y < z)){
			}	else if (is.na(x > y && y > z)) {
				}	else if (x < y && y > z) {
					peaks <- rbind(peaks, newframe[value,])
					}	else if (x > y && y < z) {
						valley <- rbind(valley, newframe[value,])
						}
		}
		
		target <- peaks[which(peaks$y == max(peaks$y)),]$x
		for (value in seq(1:nrow(valley))){
			if (!is.na(valley$x[value+1]) && valley$x[value] < target && target < valley$x[value+1]) {
				start <- valley$x[value]
				end <- (target+(target-start))
				start <- (target-(end-target))
			}
		}
		sdevmark <- (sd(seq(start:end))*0.25)

		par(mfrow=c(1,3))
		hist(comprdf$x, breaks=25, main="", xlab="", ylab="", xlim=c(0,cutoff))
		plot(pred, xlab="", ylab="")
			title(main=file, cex.main=0.75)
		points(peaks,pch=24,bg="green",col="black", cex=1.25)
		points(valley,pch=25,bg="red",col="black", cex=1.25)
		abline(v=target, col="red", lwd=1.75)
		abline(v=target-sdevmark, col="red", lwd=1.75, lty=2)
		abline(v=target+sdevmark, col="red", lwd=1.75, lty=2)
		abline(v=start, col="blue", lwd=1.75, lty=2)
		abline(v=end, col="blue", lwd=1.75, lty=2)
		plot(imrotate(img, 90))
		abline(v=median(comprdf$x)-sdevmark, col="white", lwd=1.75, lty=2)
		abline(v=median(comprdf$x)+sdevmark, col="white", lwd=1.75, lty=2)

	distance <- (target+sdevmark)-(target-sdevmark)
	Size <- append(Size, distance)
	histdata <- rbind(histdata, val$V2)
}	

dev.off()

clist <- c()
	for (nc in seq(1:ncol(histdata))){
		x <- paste("HIST", nc, sep = "")
		clist <- append(clist, x)
	}

colnames(histdata) <- clist

df <- cbind(df, as.data.frame(Size), files, as.data.frame(exclusionvals),histdata)

# Classify data
df$Quality <- NA
load("ASL.RData")
df$Quality <- predict(rf, df[,6:ncol(df)])

failure <- subset(df, Quality == "Q1" | Quality == "Q2")
success <- subset(df, Quality == "Q3" | Quality == "Q4")

# Output success/failure bank
write.table(df, "alldata.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(failure, "failures.txt", quote=FALSE, sep="\t", row.names=FALSE)
write.table(success, "success.txt", quote=FALSE, sep="\t", row.names=FALSE)

# Descriptive statistics on success
stats <- ddply(success, c("timepoint", "markfind"), summarise,
               N    = length(Size),
               mean = mean(Size),
               sd   = sd(Size),
               se   = sd / sqrt(N)
		)

write.table(stats, "data.stats.txt", sep="\t", quote=FALSE, row.names=FALSE)

totaln <- aggregate(N ~ timepoint, stats, sum)

# Graphical descriptive statistics
ggplot(stats, aes(timepoint, mean, fill=timepoint)) +
		geom_bar(stat="summary", width=0.5) +
		geom_errorbar(aes(timepoint, ymin=value-se, ymax=value+se),
			stat='summary',
			position=position_dodge(width = 0.5),
			size=0.25, width=0.25) +
		ylab("Size (pixels)") +
		xlab("") +
		labs(fill = "") +
		theme_bw()

# Groupwise comparisons
sink("TukeyHSD.txt")
totaln
cat("\n")
TukeyHSD(aov(df$Size ~ df$timepoint))
sink()
