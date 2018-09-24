# June 2018
# CVID analysis

rm(list=ls(all=T))

library(vegan)
library(ggplot2)
library(reshape2)
library(viridis)
library(ineq)
library(gridExtra)
library(grid)
library(entropy)

# cvid metadata
idmet <- read.table('/Volumes/BF_MI_2/TCRseq-sample-dets/metadata-tcrrepseq-cvid.txt', header=T, sep='\t')
# hv metadata
hvmet <- read.table('/Volumes/BF_MI_2/TCRseq-sample-dets/metadata-tcrrepseq-hv.txt', sep='\t', header=T)
matchagehv <- hvmet[which(hvmet$hvagebetafiles >= 20 & hvmet$hvagebetafiles < 80),]

pathtohv <- '/Volumes/BF_MI_1/cvid-analysis/hv_freq'
pathtochudhv <- '/Volumes/BF_MI_3/chudakov-dataset/data'
fileshv <- c(list.files(pathtohv, pattern='beta', full.names=T), list.files(pathtochudhv, pattern='beta.*freq', full.names=T))
matchhvfiles <- unlist(sapply(matchagehv$refname, function(x) grep(x, fileshv, value=T)))

pathtoid <- '/Volumes/BF_MI_1/cvid-analysis/cvid_freq'
filesid <- list.files(pathtoid, pattern='beta', full.names=T)

# read in files to get sample size, s
combfilenames <- c(matchhvfiles, filesid)
comb_abund <- sapply(combfilenames, function(f){
	dat <- read.table(f, header=F, sep=',')
	freq <- sum(dat[,6])
	freq
})

s <- round((90/100)*min(comb_abund))
rm(comb_abund)

label <- c(paste('HV', seq(1:length(matchhvfiles)), sep=''), paste('ID', seq(1:length(filesid)), sep=''))
group <- c(rep('HV', length(matchhvfiles)), rep('ID', length(filesid)))

comb_div_avg <- NULL

##### (a) #####

for(i in 1:length(combfilenames)){
	dat <- read.table(combfilenames[i], header=F, sep=',')
	dat[,7] <- do.call(paste, c(dat[,1:5], sep=':'))
	dat2 <- rep(dat[,7], dat[,6])
	df <- NULL
	for(i in 1:100){
	subs_dat <- as.vector(table(as.vector(sample(dat2, s, replace=F))))
	H_vals <- renyi(subs_dat, scales=c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf), hill=F)
	df <- rbind(df, H_vals)
	}
subs_avg <- apply(df, 2, mean)
comb_div_avg <- rbind(comb_div_avg, subs_avg)
}

rownames(comb_div_avg) <- label
comb <- as.data.frame(comb_div_avg)
comb$group <- group
comb$label <- label

dtp <- melt(comb)

ggplot(dtp, aes(x=as.factor(variable), y=value, group=label, color=group)) + geom_point(size=2, alpha=.6) + geom_line(alpha=.6) + scale_color_manual(values=c('#999999', '#E69F00')) + ylab(expression(paste('H', alpha))) + xlab(expression(alpha)) + theme_bw() + theme(text=element_text(size=16))

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_RenyiHvCvid-beta.png', sep=''), dpi=600)

##### (b) #####

for(i in 1:length(combfilenames)){
	dat <- read.table(combfilenames[i], header=F, sep=',')
	dat[,7] <- do.call(paste, c(dat[,1:5], sep=':'))
	dat2 <- rep(dat[,7], dat[,6])
	df <- NULL
	for(i in 1:100){
	subs_dat <- as.vector(table(as.vector(sample(dat2, s, replace=F))))
	x <- ineq(subs_dat, type='Gini')
	y <- entropy(subs_dat, unit='log')
	comb_vals <- cbind(x, y)
	df <- rbind(df, comb_vals)
	}
subs_avg <- apply(df, 2, mean)
comb_div_avg <- rbind(comb_div_avg, subs_avg)
}

rownames(comb_div_avg) <- label
comb <- as.data.frame(comb_div_avg)
comb$group <- group
comb$label <- label

p1 <- ggplot(comb, aes(x=as.factor(1), y=x, color=group)) + geom_boxplot() + scale_color_manual(values=c('#999999', '#E69F00')) + xlab('') + ylab('Gini index') + theme_bw() + theme(text=element_text(size=20), axis.text.x=element_blank(), axis.ticks.x=element_blank())

p2 <- ggplot(comb, aes(x=as.factor(1), y=y, color=group)) + geom_boxplot() + scale_color_manual(values=c('#999999', '#E69F00')) + xlab('') + ylab('Shannon entropy') + theme_bw() + theme(text=element_text(size=20), axis.text.x=element_blank(), axis.ticks.x=element_blank())

g <- grid_arrange_shared_legend(p1, p2)
ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_GiniShanHvCvid-beta.png', sep=''), g, dpi=600)




