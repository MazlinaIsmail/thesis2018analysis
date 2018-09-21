# Aug 2018

# effects of sampling in Paris samples
# exclude Jurkat

rm(list=ls(all=T))

library(ineq)
library(entropy)
library(vegan)
library(ggplot2)
library(reshape2)

##### Start #####

# read in all files (miseq/nextseq)
parentPath <- '/Volumes/BF_MI_3/Paris-sample-analysis'
#fnames <- list.files(parentPath, pattern=paste(chain, '.*Paris.*freq', sep=''), full.names=T, recursive=T)
fnames_alpha <- list.files(parentPath, pattern='alpha.*Paris.*freq', full.names=T, recursive=T)
fnames_beta <- list.files(parentPath, pattern='beta.*Paris.*freq', full.names=T, recursive=T)

fnames <- c(fnames_alpha, fnames_beta)
dcr <- lapply(fnames, read.table, header=F, sep=',')
# exclude Jurkat sequences
dcrnoj_alpha <- lapply(1:18, function(x)
{
	# TRA
	dcr[[x]][which(dcr[[x]][,1] != 39 & dcr[[x]][,2] != 18),]
	#TRB
	#dcr[[x]][which(dcr[[x]][,1] != 6 & dcr[[x]][,2] != 1),]
})

dcrnoj_beta <- lapply(19:36, function(x)
{
	# TRA
	#dcr[[x]][which(dcr[[x]][,1] != 39 & dcr[[x]][,2] != 18),]
	#TRB
	dcr[[x]][which(dcr[[x]][,1] != 6 & dcr[[x]][,2] != 1),]
})

dcrnoj <- c(dcrnoj_alpha, dcrnoj_beta)

##### Gini/Shannon #####

sampling_sizes <- c(500, 1000, 1500, 2000)

values <- lapply(sampling_sizes, function(s)
{
comb_div_avg <- NULL
	
for(x in 1:length(dcrnoj))
{
	dat <- data.frame(do.call(paste, c(dcrnoj[[x]][,1:5], sep=':')), dcrnoj[[x]][,6])
	
	dat2 <- rep(dat[,1], dat[,2])
	df <- NULL
	for(i in 1:100)
	{
		subs_dat <- as.vector(table(as.vector(sample(dat2, s, replace=F))))
		gini <- ineq(subs_dat, type='Gini')
		shannon <- entropy(subs_dat, unit='log')
		comb_vals <- cbind(gini, shannon)
		df <- rbind(df, comb_vals)
	}
subs_avg <- apply(df, 2, mean)
comb_div_avg <- rbind(comb_div_avg, subs_avg)
}
df <- data.frame(comb_div_avg)
df$platform <- rep(c('miseq', 'nextseq'), 2, each=9)
df$dilution <- rep(c(1000, 100, 10))
df$samplesize <- rep(s)
df$chain <- rep(c('alpha', 'beta'), each=18)
df
})

dtp <- do.call(rbind, values)

#save(dtp, file=paste('/Volumes/BF_MI_3/Paris-sample-analysis/scripts/', Sys.Date(), '_SamplingEff-GiniShannon.RData', sep=''))

dtp2 <- melt(dtp, id=c('platform', 'dilution', 'samplesize', 'chain'))

ggplot(dtp2, aes(x=samplesize, y=value, colour=platform, group=interaction(samplesize, platform))) + geom_boxplot() + scale_color_manual(values=c('#999999', '#D55E00')) + theme_bw() + theme(text=element_text(size=20)) + facet_grid(variable~chain, scales='free')

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_SamplingEff-GiniShannon.png'), dpi=600, width=10, height=5)

##### Renyi #####

sampling_sizes <- c(500, 1000, 1500, 2000)

values <- lapply(sampling_sizes, function(s)
{
comb_div_avg <- NULL
	
for(x in 1:length(dcrnoj))
{
	dat <- data.frame(do.call(paste, c(dcrnoj[[x]][,1:5], sep=':')), dcrnoj[[x]][,6])
	
	dat2 <- rep(dat[,1], dat[,2])
	df <- NULL
	for(i in 1:100)
	{
		subs_dat <- as.vector(table(as.vector(sample(dat2, s, replace=F))))
		H_vals <- renyi(subs_dat, scales=c(0, 1, 2, Inf), hill=F)
		df <- rbind(df, H_vals)
	}
subs_avg <- apply(df, 2, mean)
comb_div_avg <- rbind(comb_div_avg, subs_avg)
}
df <- data.frame(comb_div_avg)
df$platform <- rep(c('miseq', 'nextseq'), 2, each=9)
df$dilution <- rep(c(1000, 100, 10))
df$samplesize <- rep(s)
df$chain <- rep(c('alpha', 'beta'), each=18)
df
})

dtp <- do.call(rbind, values)

#save(dtp, file=paste('/Volumes/BF_MI_3/Paris-sample-analysis/scripts/', Sys.Date(), '_SamplingEff-RenyiSubset.RData', sep=''))

colnames(dtp)[1:4] <- c(0, 1, 2, Inf)

dtp2 <- melt(dtp, id=c('platform', 'dilution', 'samplesize', 'chain'))

ggplot(dtp2, aes(x=samplesize, y=value, colour=platform, group=interaction(samplesize, platform))) + geom_boxplot() + scale_color_manual(values=c('#999999', '#D55E00')) + theme_bw() + theme(text=element_text(size=20)) + facet_grid(variable~chain, scales='free')

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_SamplingEff-RenyiSubset.png'), dpi=600, width=10, height=10)

##### Plotting #####







