# Aug 2018
# most abundant TCRs in tumour compared to presence in blood
# defining 'abundant'

##### Pseudocode #####
# (1) ubiquitous, present > n times in at least 1 tr
# count number of TCRs present 2, 4, 6 etc times in tr compared to blood/ntl
# (2) proportion of ubiquitous TCRs in normal/blood
# 2x and 10x
# (3) proportion of non ubiquitous TCRs in normal/blood

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
		
		sharedtrdcr <- Reduce(intersect, lapply(1:length(trdcr2), function(x){trdcr2[[x]][,7]}))
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
		
		sharedtrdcr_df$blood <- rep(NA)
		sharedtrdcr_df$ntl <- rep(NA)
		
		if(nrow(b) == 1)
		{
			bdcr <- read.table(as.character(b$filename), header=F, sep=',')
			bdcr$d <- do.call(paste, c(bdcr[,1:5], sep=':'))
			sharedtrdcr_df$blood <- bdcr[match(sharedtrdcr_df[,1], bdcr[,7]),6]
		}
		if(nrow(n) == 1)
		{
			ndcr <- read.table(as.character(n$filename), header=F, sep=',')
			ndcr$d <- do.call(paste, c(ndcr[,1:5], sep=':'))
			sharedtrdcr_df$ntl <- ndcr[match(sharedtrdcr_df[,1], ndcr[,7]),6]
		}
		}
		
		# get number of times present in TR compared to blood/ntl
		
		sharedtrdcr_df[is.na(sharedtrdcr_df)] <- 0

		if(nrow(sharedtrdcr_df) > 0)
		{
		tr_len <- nrow(tr)
		
		c <- unlist(lapply(1:length(times), function(t)
		{
		unlist(lapply(2:tr_len+1, function(r)
		{
			shared_df <- sharedtrdcr_df[sharedtrdcr_df[,r] / sharedtrdcr_df$blood > times[t] & sharedtrdcr_df[,r] > times[t],]
			nrow(shared_df)
		}))
		}))
		b_counts[1:length(times)] <- c
		}
		
		if(nrow(sharedtrdcr_df) > 0)
		{
		tr_len <- nrow(tr)
		
		c <- unlist(lapply(1:length(times), function(t)
		{
		unlist(lapply(2:tr_len+1, function(r)
		{
			shared_df <- sharedtrdcr_df[sharedtrdcr_df[,r] / sharedtrdcr_df$ntl > times[t] & sharedtrdcr_df[,r] > times[t],]
			nrow(shared_df)
		}))
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

ggplot(dtp2, aes(x=as.factor(threshold), y=value, group=interaction(threshold, variable), colour=variable)) + geom_boxplot() + scale_colour_manual(values=c('#000000', '#0072B2')) + ylab('number of tumour ubiquitous TCRs') + xlab('threshold') + theme_bw() + theme(text=element_text(size=20)) + facet_wrap(~chain, scales='free', nrow=2)

ggsave(paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_AbundantThreshold-LTX.png', sep=''), dpi=600, width=8, height=8)

##### (2) #####

getvals <- function(tcrchain){
vals <- lapply(ltx, function(x)
{
	subset <- ltxid_filename[which(grepl(x, ltxid_filename$ltxid) == T & grepl(tcrchain, ltxid_filename$ltxid) == T),]
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

		sharedtrdcr <- Reduce(intersect, lapply(1:length(trdcr2), function(x){trdcr2[[x]][,7]}))
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
		
		sharedtrdcr_df$blood <- rep(NA)
		sharedtrdcr_df$ntl <- rep(NA)
		
		bdcr_total <- NULL
		ndcr_total <- NULL
		
		if(nrow(b) == 1)
		{
			bdcr <- read.table(as.character(b$filename), header=F, sep=',')
			bdcr$d <- do.call(paste, c(bdcr[,1:5], sep=':'))
			bdcr_total <- sum(bdcr[,6])
			sharedtrdcr_df$blood <- bdcr[match(sharedtrdcr_df[,1], bdcr[,7]),6]
		}
		if(nrow(n) == 1)
		{
			ndcr <- read.table(as.character(n$filename), header=F, sep=',')
			ndcr$d <- do.call(paste, c(ndcr[,1:5], sep=':'))
			ndcr_total <- sum(ndcr[,6])
			sharedtrdcr_df$ntl <- ndcr[match(sharedtrdcr_df[,1], ndcr[,7]),6]
		}
		}
		
		# get number of times present in TR compared to blood/ntl

		sharedtrdcr_df[is.na(sharedtrdcr_df)] <- 0
		
		if(nrow(sharedtrdcr_df) > 0)
		{
		tr_len <- length(trdcr)
		
		c <- unlist(lapply(1:length(times), function(t)
		{
		unlist(lapply(2:tr_len+1, function(r)
		{
			shared_df <- sharedtrdcr_df[sharedtrdcr_df[,r] / sharedtrdcr_df$blood > times[t] & sharedtrdcr_df[,r] > times[t],]
			ifelse(is.integer(bdcr_total), sum(shared_df$blood)/bdcr_total, NA)
		}))
		}))
		b_counts[1:length(times)] <- c
		}
		
		if(nrow(sharedtrdcr_df) > 0)
		{
		tr_len <-length(trdcr)
		
		c <- unlist(lapply(1:length(times), function(t)
		{
		unlist(lapply(2:tr_len+1, function(r)
		{
			shared_df <- sharedtrdcr_df[sharedtrdcr_df[,r] / sharedtrdcr_df$ntl > times[t] & sharedtrdcr_df[,r] > times[t],]
			ifelse(is.integer(ndcr_total), sum(shared_df$ntl)/ndcr_total, NA)
		}))
		}))
		n_counts[1:length(times)] <- c
		}
	}
	dat_df <- data.frame(b_counts, n_counts)
	dat_df$patientid <- paste(x, times, sep='_')
	dat_df$threshold <- as.factor(times)
	dat_df$chain <- rep(tcrchain)
	dat_df
})
}

valsalpha <- getvals(tcrchain='alpha')
valsbeta <- getvals(tcrchain='beta')

dtp <- do.call(rbind, c(valsalpha, valsbeta))
dtp2 <- melt(dtp, id=c('patientid', 'threshold', 'chain'))

dtp3 <- dtp2[which(dtp2$threshold == 2 | dtp2$threshold == 10),]

ggplot(dtp3, aes(x=variable, y=value)) + geom_line(aes(group=interaction(patientid, threshold), colour=threshold)) + geom_point(aes(colour=threshold)) + scale_colour_manual(values=c("#999999", "#D55E00")) + scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) + annotation_logticks(base=10) + xlab('') + ylab('proportion of tumour ubiquitous TCRs') + theme_bw() + theme(text=element_text(size=18)) + facet_wrap(~chain, scales='free', ncol=2)

ggsave(paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_LTX-tumourUbiqBloodNtl.png', sep=''), dpi=600, width=10, height=5)

##### (3) #####

getvals <- function(tcrchain){
vals <- 

lapply(ltx[2], function(x, tcrchain='alpha')
{
	subset <- ltxid_filename[which(grepl(x, ltxid_filename$ltxid) == T & grepl(tcrchain, ltxid_filename$ltxid) == T),]
	dat <- subset[grep('PBMCF', subset$ltxid, invert=T),]

	b <- dat[which(grepl('PBMC', dat$ltxid) == T),]
	n <- dat[which(grepl(':N:', dat$ltxid) == T),]
	tr <- dat[which(grepl(':R', dat$ltxid) == T),]
	
	tr })
	
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
		
		sharedtrdcr_df2 <- sharedtrdcr_df[,2:length(names(sharedtrdcr_df))]
		rownames(sharedtrdcr_df2) <- sharedtrdcr_df[,1]
		sharedtrdcr_df2[sharedtrdcr_df2 >= 1] <- 1
		sharedtrdcr_df2
		
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
		
		sharedtrdcr_df2[is.na(sharedtrdcr_df2)] <- 0
		
		tr_len <- length(trdcr)
		tr_len
		#sharedtrdcr_df2[,1:tr_len]
		# if(nrow(sharedtrdcr_df2) > 0)
		# {
			# tr_len <- length(tr)
			# sharedtrdcr_df2[which(sum(sharedtrdcr_df2[,1:tr_len] == tr_len)),]
		# }

		}})
		
		
		
		
		
		
		
		
		
		
		if(nrow(sharedtrdcr_df) > 0)
		{
		tr_len <- nrow(tr)
		
		c <- unlist(lapply(1:length(times), function(t)
		{
		unlist(lapply(2:tr_len+1, function(r)
		{
			shared_df <- sharedtrdcr_df[sharedtrdcr_df[,r] / sharedtrdcr_df$blood > times[t] & sharedtrdcr_df[,r] > times[t],]
			ifelse(is.integer(bdcr_total), sum(shared_df$blood)/bdcr_total, NA)
		}))
		}))
		b_counts[1:length(times)] <- c
		}
		
		if(nrow(sharedtrdcr_df) > 0)
		{
		tr_len <- nrow(tr)
		
		c <- unlist(lapply(1:length(times), function(t)
		{
		unlist(lapply(2:tr_len+1, function(r)
		{
			shared_df <- sharedtrdcr_df[sharedtrdcr_df[,r] / sharedtrdcr_df$ntl > times[t] & sharedtrdcr_df[,r] > times[t],]
			ifelse(is.integer(ndcr_total), sum(shared_df$ntl)/ndcr_total, NA)
		}))
		}))
		n_counts[1:length(times)] <- c
		}
	}
	dat_df <- data.frame(b_counts, n_counts)
	dat_df$patientid <- paste(x, times, sep='_')
	dat_df$threshold <- as.factor(times)
	dat_df$chain <- rep(tcrchain)
	dat_df
})
}

valsalpha <- getvals(tcrchain='alpha')
valsbeta <- getvals(tcrchain='beta')

dtp <- do.call(rbind, c(valsalpha, valsbeta))
dtp2 <- melt(dtp, id=c('patientid', 'threshold', 'chain'))

dtp3 <- dtp2[which(dtp2$threshold == 2 | dtp2$threshold == 10),]

ggplot(dtp3, aes(x=variable, y=value)) + geom_line(aes(group=interaction(patientid, threshold), colour=threshold)) + geom_point(aes(colour=threshold)) + scale_colour_manual(values=c("#999999", "#D55E00")) + scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) + annotation_logticks(base=10) + xlab('') + ylab('proportion of tumour ubiquitous TCRs') + theme_bw() + theme(text=element_text(size=18)) + facet_wrap(~chain, scales='free', ncol=2)

ggsave(paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_LTX-tumourEnrichedBloodNtl.png', sep=''), dpi=600, width=10, height=5)


##### test area #####


