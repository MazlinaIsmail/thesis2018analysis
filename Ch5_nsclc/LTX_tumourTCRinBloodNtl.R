# Aug 2018
# count proportion of ubiquitous and region specific TCRs in blood and normal

##### Description #####
# (1) - number of ubiq / total num of ubiq, in blood and ntl
#     - number of region-specific / total num of region-specific, in b and ntl
# (2) defining what is enriched tumour

library(ggplot2)
library(reshape2)

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

##### (1) #####

getvals <- function(tcrchain){
vals <- lapply(ltx, function(x)
{
	subset <- ltxid_filename[which(grepl(x, ltxid_filename$ltxid) == T & grepl(tcrchain, ltxid_filename$ltxid) == T),]
	subset

	dat <- subset[grep('PBMCF', subset$ltxid, invert=T),]

	b <- dat[which(grepl('PBMC', dat$ltxid) == T),]
	n <- dat[which(grepl(':N:', dat$ltxid) == T),]
	tr <- dat[which(grepl(':R', dat$ltxid) == T),]
	
	b_counts <- NA
	n_counts <- NA

	if(nrow(b) == 1 | nrow(n) == 1)
	{
		trdcr <- lapply(as.character(tr$filename), read.table, header=F, sep=',')
		trdcr2 <- lapply(1:length(trdcr), function(y)
		{
			d <- do.call(paste, c(trdcr[[y]][,1:5], sep=':'))
			p <- trdcr[[y]][,6]/sum(trdcr[[y]][,6])
			new_df <- cbind(trdcr[[y]], d, p)
			new_df
		})
		
		sharedtrdcr <- Reduce(union, lapply(1:length(trdcr2), function(x){trdcr2[[x]][,7]}))

		sharedtrdcr_df <- data.frame(sharedtrdcr)
		
		if(length(sharedtrdcr) > 0)
		{
		num_cols_to_add <- length(trdcr2)
		i <- 1
		while(i <= num_cols_to_add){
			data <- trdcr2[[i]]
			sharedtrdcr_df[,i+1] <- data[match(sharedtrdcr, data[,7]),6]
			i <- i+1
		}
		sharedtrdcr_df2 <- data.frame(sharedtrdcr_df[,2:length(names(sharedtrdcr_df))])
		rownames(sharedtrdcr_df2) <- sharedtrdcr_df[,1]
		sharedtrdcr_df2[sharedtrdcr_df2 >= 1] <- 1
		
		sharedtrdcr_df2$blood <- rep(NA)
		sharedtrdcr_df2$ntl <- rep(NA)
		
		bdcr_total <- NULL
		ndcr_total <- NULL
		
		if(nrow(b) == 1)
		{
			bdcr <- read.table(as.character(b$filename), header=F, sep=',')
			bdcr$d <- do.call(paste, c(bdcr[,1:5], sep=':'))
			bdcr_total <- sum(bdcr[,6])
			sharedtrdcr_df2$blood <- bdcr[match(sharedtrdcr_df[,1], bdcr[,7]),6]
		}
		if(nrow(n) == 1)
		{
			ndcr <- read.table(as.character(n$filename), header=F, sep=',')
			ndcr$d <- do.call(paste, c(ndcr[,1:5], sep=':'))
			ndcr_total <- sum(ndcr[,6])
			sharedtrdcr_df2$ntl <- ndcr[match(sharedtrdcr_df[,1], ndcr[,7]),6]
		}
		}
		
		if(nrow(sharedtrdcr_df2) > 0 & nrow(tr) > 1)
		{
		sharedtrdcr_df2[is.na(sharedtrdcr_df2)] <- 0
		sharedtrdcr_df2$tr_common <- rowSums(sharedtrdcr_df2[,1:nrow(tr)])
		# edit here for ubiquitous or not
		#ubiq <- sharedtrdcr_df2[which(sharedtrdcr_df2$tr_common == nrow(tr)),]
		ubiq <- sharedtrdcr_df2[which(sharedtrdcr_df2$tr_common < nrow(tr)),]

		b_counts <- ifelse(is.integer(bdcr_total), nrow(ubiq[which(ubiq$blood > 0),])/nrow(ubiq), NA)
		n_counts <- ifelse(is.integer(ndcr_total), nrow(ubiq[which(ubiq$ntl > 0),])/nrow(ubiq), NA)
		}
	}
	dat_df <- data.frame(b_counts, n_counts)
	dat_df$patientid <- x
	dat_df$chain <- tcrchain
	dat_df
})
}
	
valsalpha <- getvals(tcrchain='alpha')
valsbeta <- getvals(tcrchain='beta')

dtp <- do.call(rbind, c(valsalpha, valsbeta))
dtp2 <- melt(dtp, id=c('patientid', 'chain'))
	
ggplot(dtp2, aes(x=variable, y=value)) + geom_boxplot() + geom_line(aes(group=patientid), alpha=.7) + geom_point() + ylab('proportion') + xlab('') + theme_bw() + theme(text=element_text(size=20)) + facet_wrap(~chain)	
	
ggsave(paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_LTX-TNB-regionSpecificBloodNtl.png', sep=''), dpi=600, width=6, height=4)
	
##### (2) #####

getvals <- function(tcrchain){
vals <- lapply(ltx, function(x)
{
	subset <- ltxid_filename[which(grepl(x, ltxid_filename$ltxid) == T & grepl(tcrchain, ltxid_filename$ltxid) == T),]
	subset

	dat <- subset[grep('PBMCF', subset$ltxid, invert=T),]

	b <- dat[which(grepl('PBMC', dat$ltxid) == T),]
	n <- dat[which(grepl(':N:', dat$ltxid) == T),]
	tr <- dat[which(grepl(':R', dat$ltxid) == T),]
	
	times <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)
	b_counts <- rep(NA, length(times))
	n_counts <- rep(NA, length(times))
	
	if(nrow(b) == 1 | nrow(n) == 1)
	{
		trdcr <- lapply(as.character(tr$filename), read.table, header=F, sep=',')
		trdcr2 <- lapply(1:length(trdcr), function(y)
		{
			d <- do.call(paste, c(trdcr[[y]][,1:5], sep=':'))
			p <- trdcr[[y]][,6]/sum(trdcr[[y]][,6])
			new_df <- cbind(trdcr[[y]], d, p)
			new_df
		})
		
		sharedtrdcr <- Reduce(union, lapply(1:length(trdcr2), function(x){trdcr2[[x]][,7]}))

		sharedtrdcr_df <- data.frame(sharedtrdcr)
		
		if(length(sharedtrdcr) > 0)
		{
		num_cols_to_add <- length(trdcr2)
		i <- 1
		while(i <= num_cols_to_add){
			data <- trdcr2[[i]]
			sharedtrdcr_df[,i+1] <- data[match(sharedtrdcr, data[,7]),6]
			i <- i+1
		}
		sharedtrdcr_df2 <- data.frame(sharedtrdcr_df[,2:length(names(sharedtrdcr_df))])
		rownames(sharedtrdcr_df2) <- sharedtrdcr_df[,1]
		#sharedtrdcr_df2[sharedtrdcr_df2 >= 1] <- 1
		
		sharedtrdcr_df2$blood <- rep(NA)
		sharedtrdcr_df2$ntl <- rep(NA)
		
		bdcr_total <- NULL
		ndcr_total <- NULL
		
		if(nrow(b) == 1)
		{
			bdcr <- read.table(as.character(b$filename), header=F, sep=',')
			bdcr$d <- do.call(paste, c(bdcr[,1:5], sep=':'))
			bdcr_total <- sum(bdcr[,6])
			sharedtrdcr_df2$blood <- bdcr[match(sharedtrdcr_df[,1], bdcr[,7]),6]
		}
		if(nrow(n) == 1)
		{
			ndcr <- read.table(as.character(n$filename), header=F, sep=',')
			ndcr$d <- do.call(paste, c(ndcr[,1:5], sep=':'))
			ndcr_total <- sum(ndcr[,6])
			sharedtrdcr_df2$ntl <- ndcr[match(sharedtrdcr_df[,1], ndcr[,7]),6]
		}
		}
		
		if(nrow(sharedtrdcr_df2) > 0 & nrow(tr) > 1)
		{
		sharedtrdcr_df2[is.na(sharedtrdcr_df2)] <- 0
		tr_len <- nrow(tr)
		
		c <- unlist(lapply(1:length(times), function(t)
		{
			cond <- sharedtrdcr_df2[,1:length(tail(colnames(sharedtrdcr_df2), -2))]/sharedtrdcr_df2[,length(colnames(sharedtrdcr_df2))-1] > times[t] & sharedtrdcr_df2[,1:length(tail(colnames(sharedtrdcr_df2), -2))] > times[t]
			shared <- sharedtrdcr_df2[which(apply(cond, 1, any)),]
			nrow(shared)
		}))
		b_counts[1:length(times)] <- c
		}
		
		if(nrow(sharedtrdcr_df2) > 0 & nrow(tr) > 1)
		{
		sharedtrdcr_df2[is.na(sharedtrdcr_df2)] <- 0
		tr_len <- nrow(tr)
		
		c <- unlist(lapply(1:length(times), function(t)
		{
			cond <- sharedtrdcr_df2[,1:length(tail(colnames(sharedtrdcr_df2), -2))]/sharedtrdcr_df2[,length(colnames(sharedtrdcr_df2))] > times[t] & sharedtrdcr_df2[,1:length(tail(colnames(sharedtrdcr_df2), -2))] > times[t]
			shared <- sharedtrdcr_df2[which(apply(cond, 1, any)),]
			nrow(shared)
		}))
		n_counts[1:length(times)] <- c
		}
	}
	dat_df <- data.frame(b_counts, n_counts)
	dat_df$patientid <- rep(x)
	dat_df$threshold <- times
	dat_df$chain <- rep(tcrchain)
	dat_df
})
}

valsalpha <- getvals(tcrchain='alpha')
valsbeta <- getvals(tcrchain='beta')

dtp <- do.call(rbind, c(valsalpha, valsbeta))
dtp2 <- melt(dtp, id=c('patientid', 'threshold', 'chain'))

ggplot(dtp2, aes(x=as.factor(threshold), y=value, group=interaction(threshold, variable), colour=variable)) + geom_boxplot() + scale_colour_manual(values=c('#000000', '#0072B2')) + ylab(expression('number of tumour TCRs')) + xlab(expression(paste(italic('n'), '-times in blood/non-tumour lung', sep=''))) + scale_y_log10(breaks=c(10, 100, 1000)) + theme_bw() + theme(text=element_text(size=20)) + facet_wrap(~chain, scales='free', nrow=2)

ggsave(paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_AbundantThreshold-LTX.png', sep=''), dpi=600, width=8, height=8)


#### test #####

a <- c(2, 4, 6, 8, 2, 1)
b <- c(1, 3, 2, 10, 0, 1)
c <- c(1, 1, 1, 1, 1, 0)
d <- c(1, 2, 3, 1, 20, 10)

df <- data.frame(rbind(a, b, c, d))

cond <- df[,1:length(tail(colnames(df), -2))]/df[,length(colnames(df))-1] > 2 & df[,1:length(tail(colnames(df), -2))] > 2

which(apply(cond, 1, any))






	