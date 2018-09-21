# Apr 2018
# comparing frequency distribution of Paris samples pre and post collapsing

##### Pseudocode #####
# 1. read in n12 files
# 2. get DCR, collapse to get un-corrected abundance
# 3. tabulate abundance (how many 1s, 2s, 3s, etc)
# 4. plot proportion on y, size on x 

rm(list=ls(all=T))

library(ggplot2)

##### Start #####

chain <- 'beta.*Paris'

# MISEQ
# pre-collapsed
pathToN12_1 <- '/Volumes/BF_MI_1/MiSeq107/run1/Decombined'
pathToN12_2 <- '/Volumes/BF_MI_3/MiSeq109/run1/Decombined'
# post-collapsed
pathToFreq <- '/Volumes/BF_MI_3/Paris-sample-analysis/miseq_dcr'
pattF <- paste(chain, '.*freq', sep='')

# NEXTSEQ
pathToN12Ns <- '/Volumes/BF_MI_1/NextSeq009/Decombined_Paris'
pathToFreqNs <- '/Volumes/BF_MI_3/Paris-sample-analysis/nextseq_dcr'

##### get vals for MiSeq #####

decFiles_1 <- list.files(pathToN12_1, pattern=chain, full.names=T)
decFiles_2 <- list.files(pathToN12_2, pattern=chain, full.names=T)
freqFiles <- list.files(pathToFreq, pattern=pattF, full.names=T)

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
	tabDcr[[x]]$type <- rep('pre', nrow(tabDcr[[x]]))
	tabDcr[[x]]$platform <- rep('miseq', nrow(tabDcr[[x]]))
	tabDcr[[x]]
})

# for the freqs
freq <- lapply(freqFiles, read.table, header=F, sep=',')
tabFreq <- lapply(freq, function(x) data.frame(table(x[,6])))
tabFreq2 <- lapply(1:length(tabFreq), function(x)
{
	tabFreq[[x]]$prop <- NULL
	for(i in tabFreq[[x]])
	{
		tabFreq[[x]]$prop <- tabFreq[[x]]$Freq/sum(tabFreq[[x]]$Freq)
	}
	tabFreq[[x]]$type <- rep('post', nrow(tabFreq[[x]]))
	tabFreq[[x]]$platform <- rep('miseq', nrow(tabFreq[[x]]))
	tabFreq[[x]]
})

##### get vals for NextSeq #####

decFilesNs <- list.files(pathToN12Ns, pattern=chain, full.names=T, recursive=F)
freqFilesNs <- list.files(pathToFreqNs, pattern=pattF, full.names=T, recursive=F)

dcrNs <- lapply(decFilesNs, read.table, header=F, sep=',')
freqNs <- lapply(freqFilesNs, read.table, header=F, sep=',')

# n12
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
	tabDcrNs[[x]]$type <- rep('pre', nrow(tabDcrNs[[x]]))
	tabDcrNs[[x]]$platform <- rep('nextseq', nrow(tabDcrNs[[x]]))
	tabDcrNs[[x]]
})

# freq
tabFreqNs <- lapply(freqNs, function(x) data.frame(table(x[,6])))
tabFreqNs2 <- lapply(1:length(tabFreqNs), function(x)
{
	tabFreqNs[[x]]$prop <- NULL
	for(i in tabFreqNs[[x]])
	{
		tabFreqNs[[x]]$prop <- tabFreqNs[[x]]$Freq/sum(tabFreqNs[[x]]$Freq)
	}
	tabFreqNs[[x]]$type <- rep('post', nrow(tabFreqNs[[x]]))
	tabFreqNs[[x]]$platform <- rep('nextseq', nrow(tabFreqNs[[x]]))
	tabFreqNs[[x]]
})


##### Plotting #####

dtp1 <- do.call(rbind, tabDcrNs2)
dtp2 <- do.call(rbind, tabFreqNs2)
dtp3 <- do.call(rbind, tabDcr2)
dtp4 <- do.call(rbind, tabFreq2)
dtp <- rbind(dtp1, dtp2, dtp3, dtp4)

p <- ggplot(dtp, aes(x=as.integer(Var1), y=prop, color=platform, group=type))
p + geom_point(aes(shape=type), alpha=.6, size=3, stroke=1) + scale_x_log10() + scale_y_log10() + scale_color_manual(values=c('#999999', '#E69F00')) + scale_shape_manual(values=c(1, 4)) + theme_bw() + xlab('number of copies of each TCR') + ylab('proportion') + labs(title=strsplit(chain, '[.]')[[1]][1]) + theme(text=element_text(size=16))

ggsave(paste('/Volumes/BF_MI_3/Paris-sample-analysis/', Sys.Date(), '_FreqDistPrePostParis_', strsplit(chain, '[.]')[[1]][1], '.png', sep=''), width=7, height=5)


