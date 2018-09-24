# Aug 2018
# intra-patient pairwise similarity (TvT, TvNTL, TvB, BvNTL)
# using (1) Jaccard, (2) Spearman and (3) dot product

library(ggplot2)
library(gtools)

##### Global variables #####

parentdir <- '/Volumes/BF_MI_3/cancer-analysis_2/TRACERx_TCR_up_to_NextSeq017'
ltxdirs <- grep('LTX', list.dirs(parentdir, full.names=T), value=T)

# files used for manuscript
fileId <- read.table('/Volumes/BF_MI_3/cancer-analysis_2/TracerX_TCR_paper/TRACERx_paper1/contributor-Data/2018-08-13_useforpaper.csv', header=T, sep=',')

fileLst <- NULL
for(i in 1:length(fileId$Freq_File_Name))
{
	allfiles <- list.files(ltxdirs, pattern='freq', full.names=T, recursive=T)
	match <- grep(fileId$Freq_File_Name[i], allfiles, value=T)
	fileLst[i] <- match
}

ltxid_filename <- data.frame(fileId$LTX_ID, fileLst)
colnames(ltxid_filename) <- c('ltxid', 'filename')

ltx <- unique(unlist(lapply(ltxid_filename$ltxid, function(x) strsplit(as.character(x), ':')[[1]][1])))

##### (1) jaccard #####

values <- lapply(ltx, function(y)
{
	subset <- ltxid_filename[which(grepl(y, ltxid_filename$ltxid) == T),]
	dat <- subset[grep('PBMCF', subset$ltxid, invert=T),]
	combins <- combn(1:nrow(dat), 2, simplify=F)
	
	dcr <- lapply(as.character(dat$filename), read.table, header=F, sep=',')
	dcr2 <- lapply(1:length(dcr), function(x) do.call(paste, c(dcr[[x]][,1:5], sep=':')))
	
	# jaccard
	jiVals <- sapply(1:length(combins), function(x) length(intersect(dcr2[[combins[[x]][1]]], dcr2[[combins[[x]][2]]]))/length(union(dcr2[[combins[[x]][1]]], dcr2[[combins[[x]][2]]])))

	df <- lapply(1:length(combins), function(x){data.frame(dat[combins[[x]][1],1], dat[combins[[x]][2],1])})
	df2 <- do.call(rbind, df)
	colnames(df2) <- c('var1', 'var2')
	df2$jaccardVals <- jiVals
	df2$patient <- rep(y)
	df2$chain1 <- unlist(lapply(df2$var1, function(x) strsplit(as.character(x), ':')[[1]][3]))
	df2$chain2 <- unlist(lapply(df2$var2, function(x) strsplit(as.character(x), ':')[[1]][3]))
	df2$type1 <- ifelse(grepl(':R', df2$var1) == T, 'T', ifelse(grepl(':N', df2$var1) == T, 'NTL', 'B'))
	df2$type2 <- ifelse(grepl(':R', df2$var2) == T, 'T', ifelse(grepl(':N', df2$var2) == T, 'NTL', 'B'))
	df2[which(df2$chain1 == df2$chain2),]
})

dtp <- do.call(rbind, values)
p <- combinations(3, 2, c('T', 'NTL', 'B'), repeats=T)

dtp$pairwise <- unlist(lapply(1:nrow(dtp), function(x)
{
	lapply(1:nrow(p), function(y)
	{
		if(as.character(dtp[x,7]) == p[y,1] & as.character(dtp[x,8]) == p[y,2] | as.character(dtp[x,7]) == p[y,2] & as.character(dtp[x,8]) == p[y,1])
		paste(p[y,], collapse='.')
	})
}))

ggplot(dtp, aes(x=patient, y=jaccardVals, colour=pairwise)) + geom_point(size=2) + scale_colour_manual(values=c("#999999", "#E69F00", "#009E73", "#CC79A7")) + guides(colour=guide_legend(title="interaction")) + ylab('Similarity (Jaccard index)') + theme_bw() + theme(text=element_text(size=20), axis.text.x=element_text(angle=90, size=8)) + facet_wrap(~chain1, nrow=2, scales='free')

ggsave(paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_Jaccard-LTX-TNB.png', sep=''), dpi=600, width=10, height=8)

ggplot(dtp, aes(x=pairwise, y=jaccardVals)) + geom_boxplot() + ylab('Similarity (Jaccard index)') + theme_bw() + theme(text=element_text(size=20)) + facet_wrap(~chain1, scales='free')

ggsave(paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_JaccardBoxPlot-LTX-TNB.png', sep=''), dpi=600, width=8, height=4)

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# significance
labs <- combn(c('B.NTL', 'B.T', 'NTL.T', 'T.T'), 2, simplify=F)

getSig <- function(chain, dtp)
{
	lapply(1:length(labs), function(x)
	{
		z <- dtp[which(dtp$chain1 == chain & dtp$pairwise == labs[[x]][1]),3]
		y <- dtp[which(dtp$chain1 == chain & dtp$pairwise == labs[[x]][2]),3]
		
		wilcox.test(z, y)
	})
}

getSig(chain='alpha', dtp)
getSig(chain='beta', dtp)

##### (2) Spearman #####

values <- lapply(ltx, function(y)
{
	subset <- ltxid_filename[which(grepl(y, ltxid_filename$ltxid) == T),]
	dat <- subset[grep('PBMCF', subset$ltxid, invert=T),]
	combins <- combn(1:nrow(dat), 2, simplify=F)
	
	dcr <- lapply(as.character(dat$filename), read.table, header=F, sep=',')
	
	dcr2 <- lapply(1:length(dcr), function(x)
	{
		d <- do.call(paste, c(dcr[[x]][,1:5], sep=':'))
		p <- dcr[[x]][,6]/sum(dcr[[x]][,6])
		new_df <- cbind(dcr[[x]], d, p)
		new_df
	})
	
	spVals <- unlist(lapply(1:length(combins), function(x)
	{
		sharedDcr <- intersect(dcr2[[combins[[x]][1]]]$d, dcr2[[combins[[x]][2]]]$d)
		a <- dcr2[[combins[[x]][1]]][match(sharedDcr, dcr2[[combins[[x]][1]]][,7]),8]
		b <- dcr2[[combins[[x]][2]]][match(sharedDcr, dcr2[[combins[[x]][2]]][,7]),8]
		cor(a, b, method='spearman')
	}))
	
	df <- lapply(1:length(combins), function(x){data.frame(dat[combins[[x]][1],1], dat[combins[[x]][2],1])})
	df2 <- do.call(rbind, df)
	colnames(df2) <- c('var1', 'var2')
	df2$spearmanVals <- spVals
	df2$patient <- rep(y)
	df2$chain1 <- unlist(lapply(df2$var1, function(x) strsplit(as.character(x), ':')[[1]][3]))
	df2$chain2 <- unlist(lapply(df2$var2, function(x) strsplit(as.character(x), ':')[[1]][3]))
	df2$type1 <- ifelse(grepl(':R', df2$var1) == T, 'T', ifelse(grepl(':N', df2$var1) == T, 'NTL', 'B'))
	df2$type2 <- ifelse(grepl(':R', df2$var2) == T, 'T', ifelse(grepl(':N', df2$var2) == T, 'NTL', 'B'))
	df2[which(df2$chain1 == df2$chain2),]
})

dtp <- do.call(rbind, values)
p <- combinations(3, 2, c('T', 'NTL', 'B'), repeats=T)

dtp$pairwise <- unlist(lapply(1:nrow(dtp), function(x)
{
	lapply(1:nrow(p), function(y)
	{
		if(as.character(dtp[x,7]) == p[y,1] & as.character(dtp[x,8]) == p[y,2] | as.character(dtp[x,7]) == p[y,2] & as.character(dtp[x,8]) == p[y,1])
		paste(p[y,], collapse='.')
	})
}))

ggplot(dtp, aes(x=patient, y=spearmanVals, colour=pairwise)) + geom_point(size=2) + scale_colour_manual(values=c("#999999", "#E69F00", "#009E73", "#CC79A7")) + guides(colour=guide_legend(title="interaction")) + ylab('Similarity (Spearman correlation)') + theme_bw() + theme(text=element_text(size=20), axis.text.x=element_text(angle=90, size=8)) + facet_wrap(~chain1, nrow=2, scales='free')

ggsave(paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_Spearman-LTX-TNB.png', sep=''), dpi=600, width=10, height=8)

ggplot(dtp, aes(x=pairwise, y=spearmanVals)) + geom_boxplot() + ylab('Similarity (Spearman correlation)') + theme_bw() + theme(text=element_text(size=20)) + facet_wrap(~chain1, scales='free')

ggsave(paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_SpearmanBoxPlot-LTX-TNB.png', sep=''), dpi=600, width=8, height=4)

# significance
labs <- combn(c('B.NTL', 'B.T', 'NTL.T', 'T.T'), 2, simplify=F)

getSig <- function(chain, dtp)
{
	lapply(1:length(labs), function(x)
	{
		z <- dtp[which(dtp$chain1 == chain & dtp$pairwise == labs[[x]][1]),3]
		y <- dtp[which(dtp$chain1 == chain & dtp$pairwise == labs[[x]][2]),3]
		
		wilcox.test(z, y)
	})
}

getSig(chain='alpha', dtp)
getSig(chain='beta', dtp)

##### (3) dot product #####

values <- lapply(ltx, function(y)
{
	subset <- ltxid_filename[which(grepl(y, ltxid_filename$ltxid) == T),]
	dat <- subset[grep('PBMCF', subset$ltxid, invert=T),]
	combins <- combn(1:nrow(dat), 2, simplify=F)
	
	dcr <- lapply(as.character(dat$filename), read.table, header=F, sep=',')
	
	dcr2 <- lapply(1:length(dcr), function(x)
	{
		d <- do.call(paste, c(dcr[[x]][,1:5], sep=':'))
		p <- dcr[[x]][,6]/sum(dcr[[x]][,6])
		new_df <- cbind(dcr[[x]], d, p)
		new_df
	})
	
	dpVals <- unlist(lapply(1:length(combins), function(x)
	{
		sharedDcr <- intersect(dcr2[[combins[[x]][1]]]$d, dcr2[[combins[[x]][2]]]$d)
		a <- dcr2[[combins[[x]][1]]][match(sharedDcr, dcr2[[combins[[x]][1]]][,7]),8]
		b <- dcr2[[combins[[x]][2]]][match(sharedDcr, dcr2[[combins[[x]][2]]][,7]),8]
		
		(a %*% b)/(sqrt(a %*% a)*sqrt(b %*% b))
	}))

	df <- lapply(1:length(combins), function(x){data.frame(dat[combins[[x]][1],1], dat[combins[[x]][2],1])})
	df2 <- do.call(rbind, df)
	colnames(df2) <- c('var1', 'var2')
	df2$dotprodVals <- dpVals
	df2$patient <- rep(y)
	df2$chain1 <- unlist(lapply(df2$var1, function(x) strsplit(as.character(x), ':')[[1]][3]))
	df2$chain2 <- unlist(lapply(df2$var2, function(x) strsplit(as.character(x), ':')[[1]][3]))
	df2$type1 <- ifelse(grepl(':R', df2$var1) == T, 'T', ifelse(grepl(':N', df2$var1) == T, 'NTL', 'B'))
	df2$type2 <- ifelse(grepl(':R', df2$var2) == T, 'T', ifelse(grepl(':N', df2$var2) == T, 'NTL', 'B'))
	df2[which(df2$chain1 == df2$chain2),]
})

dtp <- do.call(rbind, values)
p <- combinations(3, 2, c('T', 'NTL', 'B'), repeats=T)

dtp$pairwise <- unlist(lapply(1:nrow(dtp), function(x)
{
	lapply(1:nrow(p), function(y)
	{
		if(as.character(dtp[x,7]) == p[y,1] & as.character(dtp[x,8]) == p[y,2] | as.character(dtp[x,7]) == p[y,2] & as.character(dtp[x,8]) == p[y,1])
		paste(p[y,], collapse='.')
	})
}))

ggplot(dtp, aes(x=patient, y=dotprodVals, colour=pairwise)) + geom_point(size=2) + scale_colour_manual(values=c("#999999", "#E69F00", "#009E73", "#CC79A7")) + guides(colour=guide_legend(title="interaction")) + ylab('Similarity (dot product)') + theme_bw() + theme(text=element_text(size=20), axis.text.x=element_text(angle=90, size=8)) + facet_wrap(~chain1, nrow=2, scales='free')

ggsave(paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_Dotprod-LTX-TNB.png', sep=''), dpi=600, width=10, height=8)

ggplot(dtp, aes(x=pairwise, y=dotprodVals)) + geom_boxplot() + ylab('Similarity (dot product)') + theme_bw() + theme(text=element_text(size=20)) + facet_wrap(~chain1, scales='free')

ggsave(paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_DotprodBoxPlot-LTX-TNB.png', sep=''), dpi=600, width=8, height=4)

# significance
labs <- combn(c('B.NTL', 'B.T', 'NTL.T', 'T.T'), 2, simplify=F)

getSig <- function(chain, dtp)
{
	lapply(1:length(labs), function(x)
	{
		z <- dtp[which(dtp$chain1 == chain & dtp$pairwise == labs[[x]][1]),3]
		y <- dtp[which(dtp$chain1 == chain & dtp$pairwise == labs[[x]][2]),3]
		
		wilcox.test(z, y)
	})
}

getSig(chain='alpha', dtp)
getSig(chain='beta', dtp)


##### Scribbles #####

               var1            var2 jaccardVals patient chain1 chain2 type1 type2
161 LTX073:R2:alpha LTX073:R3:alpha   0.8107164  LTX073  alpha  alpha     T     T



	







