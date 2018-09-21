# May 2018
# compare sequence data from two protocols (dcr and mixcr)
# both processed using dcr pipeline

library(ggplot2)

# pre-collapsed
# DCR
pathToN12_1 <- '/Volumes/BF_MI_1/MiSeq107/run1/Decombined'
pathToN12_2 <- '/Volumes/BF_MI_3/MiSeq109/run1/Decombined'
chain1 <- 'alpha.*Paris'
decFiles_1 <- list.files(pathToN12_1, pattern=chain1, full.names=T)
decFiles_2 <- list.files(pathToN12_2, pattern=chain1, full.names=T)
decFiles <- c(decFiles_1, decFiles_2)
dcr <- lapply(decFiles, read.table, header=F, sep=',')
dcr2 <- lapply(1:length(dcr), function(x) do.call(paste, c(dcr[[x]][,1:5], sep=':')))

tabDcr <- lapply(dcr2, function(x) data.frame(table(data.frame(table(x))[,2])))
# getprop
tabDcr2 <- lapply(1:length(tabDcr), function(x)
{
	tabDcr[[x]]$prop <- NULL
	for(i in tabDcr[[x]])
	{
		tabDcr[[x]]$prop <- tabDcr[[x]]$Freq/sum(tabDcr[[x]]$Freq)
	}
	#tabDcr[[x]]$type <- rep('pre', nrow(tabDcr[[x]]))
	tabDcr[[x]]$protocol <- rep('dcr', nrow(tabDcr[[x]]))
	tabDcr[[x]]
})

# MIXCR
pathton12mixcr <- '/Volumes/BF_MI_3/Paris-sample-analysis/dcr_mixcrdat'
chain2 <- 'dcr_alpha.*n12'
decFilesNs <- list.files(pathton12mixcr, pattern=chain2, full.names=T)
dcrNs <- lapply(decFilesNs, read.table, header=F)
dcrNs2 <- lapply(1:length(dcrNs), function(x) do.call(paste, c(dcrNs[[x]][,1:5], sep=':')))

tabDcrNs <- lapply(dcrNs2, function(x) data.frame(table(data.frame(table(x))[,2])))
# getprop
tabDcrNs2 <- lapply(1:length(tabDcrNs), function(x)
{
	tabDcrNs[[x]]$prop <- NULL
	for(i in tabDcrNs[[x]])
	{
		tabDcrNs[[x]]$prop <- tabDcrNs[[x]]$Freq/sum(tabDcrNs[[x]]$Freq)
	}
	#tabDcrNs[[x]]$type <- rep('pre', nrow(tabDcrNs[[x]]))
	tabDcrNs[[x]]$protocol <- rep('mixcr', nrow(tabDcrNs[[x]]))
	tabDcrNs[[x]]
})

##### Plotting #####

dtp1 <- do.call(rbind, tabDcrNs2)
dtp2 <- do.call(rbind, tabDcr2)
dtp <- rbind(dtp1, dtp2)

p <- ggplot(dtp, aes(x=as.integer(Var1), y=prop, color=protocol))
p + geom_point() + scale_x_log10() + scale_y_log10() + scale_color_manual(values=c('#999999', '#E69F00')) + theme_bw() + xlab('number of copies of each TCR') + ylab('proportion') + labs(title=strsplit(chain1, '[.]')[[1]][1]) + theme(text=element_text(size=16))

##############################################################################################

# post-collapsed

chain1 <- 'alpha.*Paris'
pathToFreq <- '/Volumes/BF_MI_3/Paris-sample-analysis/miseq_dcr'
pattF <- paste(chain1, '.*freq', sep='')
freqFiles <- list.files(pathToFreq, pattern=pattF, full.names=T)
freq <- lapply(freqFiles, read.table, header=F, sep=',')
tabFreq <- lapply(freq, function(x) data.frame(table(x[,6])))
tabFreq2 <- lapply(1:length(tabFreq), function(x)
{
	tabFreq[[x]]$prop <- NULL
	for(i in tabFreq[[x]])
	{
		tabFreq[[x]]$prop <- tabFreq[[x]]$Freq/sum(tabFreq[[x]]$Freq)
	}
	#tabFreq[[x]]$type <- rep('post', nrow(tabFreq[[x]]))
	tabFreq[[x]]$protocol <- rep('dcr', nrow(tabFreq[[x]]))
	tabFreq[[x]]
})

pathToFreqNs <- '/Volumes/BF_MI_3/Paris-sample-analysis/dcr_mixcrdat'
freqFilesNs <- list.files(pathToFreqNs, pattern='alpha.*freq', full.names=T)
freqNs <- lapply(freqFilesNs, read.table, header=F, sep=',')

# freq
tabFreqNs <- lapply(freqNs, function(x) data.frame(table(x[,6])))
tabFreqNs2 <- lapply(1:length(tabFreqNs), function(x)
{
	tabFreqNs[[x]]$prop <- NULL
	for(i in tabFreqNs[[x]])
	{
		tabFreqNs[[x]]$prop <- tabFreqNs[[x]]$Freq/sum(tabFreqNs[[x]]$Freq)
	}
	#tabFreqNs[[x]]$type <- rep('post', nrow(tabFreqNs[[x]]))
	tabFreqNs[[x]]$protocol <- rep('mixcr', nrow(tabFreqNs[[x]]))
	tabFreqNs[[x]]
})

dtp2 <- do.call(rbind, tabFreqNs2)
dtp4 <- do.call(rbind, tabFreq2)
dtp <- rbind(dtp2, dtp4)

p <- ggplot(dtp, aes(x=as.integer(Var1), y=prop, color=protocol))
p + geom_point() + scale_x_log10() + scale_y_log10() + scale_color_manual(values=c('#999999', '#E69F00')) + theme_bw() + xlab('number of copies of each TCR') + ylab('proportion') + labs(title=strsplit(chain1, '[.]')[[1]][1]) + theme(text=element_text(size=16))

ggsave(paste('/Volumes/BF_MI_3/Paris-sample-analysis/', Sys.Date(), '_DcrVsDcrMixcrdat-FreqDistAlpha.png', sep=''), height=5, width=6, dpi=600)

##############################################################################################
# jurkat prop

dpath <- '/Volumes/BF_MI_3/Paris-sample-analysis/miseq_dcr'
mpath <- '/Volumes/BF_MI_3/Paris-sample-analysis/dcr_mixcrdat'

dfilesalpha <- list.files(dpath, full.names=T, pattern='alpha.*freq')
mfilesalpha <- list.files(mpath, full.names=T, pattern='alpha.*freq')
dfilesbeta <- list.files(dpath, full.names=T, pattern='beta.*freq')
mfilesbeta <- list.files(mpath, full.names=T, pattern='beta.*freq')

allfiles <- list(dfilesalpha, mfilesalpha, dfilesbeta, mfilesbeta)
pflabs <- c('dcr', 'mixcr', 'dcr', 'mixcr')
dcrdils <- rep(c(1000, 100, 10), 3)
mixcrdils <- rep(c(10, 100, 1000), each=3)

getvalsalpha <- lapply(1:2, function(iter1)
{
	dat <- lapply(allfiles[[iter1]], read.table, header=F, sep=',')
	proportion <- unlist(lapply(1:length(dat), function(iter2)
	{
		t <- sum(dat[[iter2]][,6])
		sum(dat[[iter2]][which(dat[[iter2]][,1] == 39 & dat[[iter2]][,2] == 18),6])/t
	}))
	df <- data.frame(proportion)
	df$protocol <- rep(pflabs[iter1], nrow(df))
	#df$dilution <- rep(dcrdils, 3)
	df$chain <- rep('TRA', nrow(df))
	df
})

getvalsbeta <- lapply(3:4, function(iter1)
{
	dat <- lapply(allfiles[[iter1]], read.table, header=F, sep=',')
	proportion <- unlist(lapply(1:length(dat), function(iter2)
	{
		t <- sum(dat[[iter2]][,6])
		sum(dat[[iter2]][which(dat[[iter2]][,1] == 6 & dat[[iter2]][,2] == 1),6])/t
	}))
	df <- data.frame(proportion)
	df$protocol <- rep(pflabs[iter1], nrow(df))
	#df$dilution <- rep(dils, 3)
	df$chain <- rep('TRB', nrow(df))
	df
})

valsalpha <- do.call(rbind, getvalsalpha)
valsbeta <- do.call(rbind, getvalsbeta)

vals <- rbind(valsalpha, valsbeta)
vals$dilution <- c(dcrdils, mixcrdils, dcrdils, mixcrdils)

p <- ggplot(vals, aes(x=dilution, y=proportion, color=chain, group=interaction(protocol, chain)))
p + geom_point(aes(shape=protocol), alpha=.7, position=position_jitter(w=.1), size=4, stroke=2) + scale_x_log10() + scale_color_manual(values=c('#999999', '#D55E00')) + scale_y_log10(limits=c(0.0001, 0.5)) + theme_bw() + theme(text=element_text(size=16)) + labs(title='Dcr')

p + geom_point(aes(shape=protocol), alpha=.7, position=position_jitter(w=.1), size=4, stroke=2) + scale_x_log10() + scale_color_manual(values=c('#999999', '#D55E00')) + theme_bw() + theme(text=element_text(size=16)) + labs(title='Dcr')



ggsave(paste('/Volumes/BF_MI_3/Paris-sample-analysis/', Sys.Date(), '_DcrVsDcrMixcrdat-JurkatProp.png', sep=''), height=5, width=5, dpi=600)





