# Apr 2018
# analysis of Paris samples

##### Pseudocode #####
# (1) similarity - Jaccard index
# (2) pearson correlation
#	  - read in every two files
#	  - get common dcrs 
#	  - for every commmon dcr, get corresponding abundance for the two files
#	  - get pearson correlation for the two
# plot as triangle heatmap

rm(list=ls(all=T))

library(ggplot2)
library(viridis)

##### Global variables #####

chain <- 'alpha'

# read in all files (miseq/nextseq)
parentPath <- '/Volumes/BF_MI_3/Paris-sample-analysis'
fnames <- list.files(parentPath, pattern=paste(chain, '.*Paris.*freq', sep=''), full.names=T, recursive=T)
labels <- unlist(lapply(fnames, function(x) paste(strsplit(strsplit(x, '/')[[1]][5], '_')[[1]][1], paste(strsplit(strsplit(x, '/')[[1]][6], '_')[[1]][3:4], collapse='_'), sep='_')))

dcr <- lapply(fnames, read.table, header=F, sep=',')
# exclude Jurkat sequences
dcrnoj <- lapply(1:length(dcr), function(x)
{
	# TRA
	#dcr[[x]][which(dcr[[x]][,1] != 39 & dcr[[x]][,2] != 18),]
	#TRB
	dcr[[x]][which(dcr[[x]][,1] != 6 & dcr[[x]][,2] != 1),]
})

##### (1) #####

dcr2 <- lapply(1:length(dcrnoj), function(x) do.call(paste, c(dcrnoj[[x]][,1:5], sep=':')))

combins <- combn(1:length(dcr2), 2, simplify=F)

# jaccard
jiVals <- sapply(1:length(combins), function(x) length(intersect(dcr2[[combins[[x]][1]]], dcr2[[combins[[x]][2]]]))/length(union(dcr2[[combins[[x]][1]]], dcr2[[combins[[x]][2]]])))

combinsNames <- lapply(1:length(combins), function(x) labels[combins[[x]]])
dtp <- data.frame(do.call(rbind, combinsNames), log(jiVals, base=2))
dtp <- data.frame(do.call(rbind, combinsNames), jiVals)
colnames(dtp)[3] <- c('jiVals')

p <- ggplot(dtp, aes(x=X1, y=X2, fill=jiVals)) 
p + geom_tile(color='white') + theme_bw() + theme(axis.text.x=element_text(angle=90, vjust=.5, hjust=.5), text=element_text(size=16)) + scale_fill_viridis(option='cividis') + labs(title=chain) + xlab('') + ylab('')

ggsave(paste('/Volumes/BF_MI_3/Paris-sample-analysis/', Sys.Date(), '_JaccardExcJurkParis_', chain, '.png', sep=''))

##### (2) #####

dcr2 <- lapply(1:length(dcrnoj), function(x)
{
	dcrnoj[[x]]$dcrnoj <- NULL
	dcrnoj[[x]]$prop <- NULL
	for(i in dcrnoj[[x]])
	{
		d <- do.call(paste, c(dcrnoj[[x]][,1:5], sep=':'))
		dcrnoj[[x]]$dcrnoj <- d
		dcrnoj[[x]]$prop <- dcrnoj[[x]]$V6/sum(dcrnoj[[x]]$V6)
	}
	dcrnoj[[x]]
})

combins <- combn(1:length(dcr2), 2, simplify=F)

pearsonVals <- unlist(lapply(1:length(combins), function(x)
{
	sharedDcr <- intersect(dcr2[[combins[[x]][1]]]$dcrnoj, dcr2[[combins[[x]][2]]]$dcrnoj)
	a <- dcr2[[combins[[x]][1]]][match(sharedDcr, dcr2[[combins[[x]][1]]][,7]),8]
	b <- dcr2[[combins[[x]][2]]][match(sharedDcr, dcr2[[combins[[x]][2]]][,7]),8]
	cor(a, b, method='pearson')
}))

combinsNames <- lapply(1:length(combins), function(x) labels[combins[[x]]])
dtp <- data.frame(do.call(rbind, combinsNames), pearsonVals)

p <- ggplot(dtp, aes(x=X1, y=X2, fill=pearsonVals)) 
p + geom_tile(color='black') + theme_bw() + theme(axis.text.x=element_text(angle=90, vjust=.5, hjust=.5), text=element_text(size=16)) + scale_fill_viridis() + labs(title=chain) + xlab('') + ylab('')

ggsave(paste('/Volumes/BF_MI_3/Paris-sample-analysis/', Sys.Date(), '_PearsonExcJurkParis_', chain, '.png', sep=''))



