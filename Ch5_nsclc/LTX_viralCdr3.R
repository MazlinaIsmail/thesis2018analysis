# Aug 2018
# look for viral associated CDR3s in LTX tumour, normal, blood

library(ggplot2)

##### Global variables #####

parentdir <- '/Volumes/BF_MI_3/cancer-analysis_2/TRACERx_TCR_up_to_NextSeq017'
ltxdirs <- grep('LTX', list.dirs(parentdir, full.names=T), value=T)

# files used for manuscript
fileId <- read.table('/Volumes/BF_MI_3/cancer-analysis_2/TracerX_TCR_paper/TRACERx_paper1/contributor-Data/2018-08-13_useforpaper.csv', header=T, sep=',')

fileLst <- NULL
for(i in 1:length(fileId$Freq_File_Name))
{
	allfiles <- list.files(ltxdirs, pattern='\\.cdr3', full.names=T, recursive=T)
	match <- grep(strsplit(as.character(fileId$Freq_File_Name[i]), '\\.')[[1]][1], allfiles, value=T)
	fileLst[i] <- match
}

ltxid_filename <- data.frame(fileId$LTX_ID, fileLst)
colnames(ltxid_filename) <- c('ltxid', 'filename')

ltx <- unique(unlist(lapply(ltxid_filename$ltxid, function(x) strsplit(as.character(x), ':')[[1]][1])))

##### (1) #####

# read in viral sets
cmvfile1 <- read.table('/Volumes/BF_MI_3/cancer-analysis_2/add_datafiles/viral-CDR3/McPAS-TCR_search-CMV.csv', header=T, sep=',')
cmvfile2 <- read.table('/Volumes/BF_MI_3/cancer-analysis_2/add_datafiles/viral-CDR3/CMV_TcRs_fromYS.csv', header=T, sep=',')
cmvfile3 <- read.table('/Volumes/BF_MI_3/cancer-analysis_2/add_datafiles/viral-CDR3/VDJdb_2018-07-04_CMV.tsv', header=T, sep='\t')
cmvcdr3 <- unique(unlist(list(cmvfile1$CDR3.beta.aa, cmvfile2$CDR3, cmvfile3[which(cmvfile3[,2] == 'TRB'),3])))
#viralcdr3 <- cmvcdr3

ebvfile1 <- read.table('/Volumes/BF_MI_3/cancer-analysis_2/add_datafiles/viral-CDR3/McPAS-TCR_search-EBV.csv', header=T, sep=',')
ebvfile2 <- read.table('/Volumes/BF_MI_3/cancer-analysis_2/add_datafiles/viral-CDR3/VDJdb_2018-07-04_EBV.tsv', header=T, sep='\t')
#viralcdr3 <- unique(unlist(list(ebvfile1$CDR3.beta.aa, ebvfile2[which(ebvfile2[,2] == 'TRB'),3])))
ebvcdr3 <- unique(unlist(list(ebvfile1$CDR3.beta.aa, ebvfile2[which(ebvfile2[,2] == 'TRB'),3])))

flufile1 <- read.table('/Volumes/BF_MI_3/cancer-analysis_2/add_datafiles/viral-CDR3/McPAS-TCR_search-Flu.csv', header=T, sep=',')
flufile2 <- read.table('/Volumes/BF_MI_3/cancer-analysis_2/add_datafiles/viral-CDR3/VDJdb_2018-07-04_Inf.tsv', header=T, sep='\t')
#viralcdr3 <- unique(unlist(list(flufile1$CDR3.beta.aa, flufile2[which(flufile2[,2] == 'TRB'),3])))
flucdr3 <- unique(unlist(list(flufile1$CDR3.beta.aa, flufile2[which(flufile2[,2] == 'TRB'),3])))

hivfile1 <- read.table('/Volumes/BF_MI_3/cancer-analysis_2/add_datafiles/viral-CDR3/McPAS-TCR_search-HIV.csv', header=T, sep=',')
hivfile2 <- read.table('/Volumes/BF_MI_3/cancer-analysis_2/add_datafiles/viral-CDR3/VDJdb-2018-08-21_HIV1.tsv', header=T, sep='\t')
#viralcdr3 <- unique(unlist(list(flufile1$CDR3.beta.aa, flufile2[which(flufile2[,2] == 'TRB'),3])))
hivcdr3 <- unique(unlist(list(hivfile1$CDR3.beta.aa, hivfile2[which(hivfile2[,2] == 'TRB'),3])))

getvals <- function(tcrchain, viralseq, v){
vals <- lapply(ltx, function(x)
{
	subset <- ltxid_filename[which(grepl(x, ltxid_filename$ltxid) == T & grepl(tcrchain, ltxid_filename$ltxid) == T),]
	dat <- subset[grep('PBMCF', subset$ltxid, invert=T),]

	b <- dat[which(grepl('PBMC', dat$ltxid) == T),]
	n <- dat[which(grepl(':N:', dat$ltxid) == T),]
	tr <- dat[which(grepl(':R', dat$ltxid) == T),]

	b_prop <- NA
	n_prop <- NA
	tr_prop <- rep(NA, nrow(tr))
	
	for(y in 1:nrow(tr))
	{
		cdr <- read.table(as.character(tr[y,2]), header=F, sep=',')
		total <- sum(cdr[,2])
		common <- sum(cdr[cdr[,1] %in% viralseq,2])
		tr_prop[y] <- common/total
	}

	if(nrow(b) == 1 | nrow(n) == 1)
	{
		if(nrow(b) == 1)
		{
			bcdr <- read.table(as.character(b$filename), header=F, sep=',')
			total <- sum(bcdr[,2])
			common <- sum(bcdr[bcdr[,1] %in% viralseq,2])
			b_prop <- common/total
		}
		if(nrow(n) == 1)
		{
			ncdr <- read.table(as.character(n$filename), header=F, sep=',')
			total <- sum(ncdr[,2])
			common <- sum(ncdr[ncdr[,1] %in% viralseq,2])
			n_prop <- common/total
		}
	}
	prop <- c(b_prop, n_prop, tr_prop)
	type <- c('blood', 'non-tumour lung', rep('tumour', length(tr_prop)))
	df <- data.frame(type, prop)
	df$pathogen <- rep(v)
	df$patientid <- rep(x)
	df
})}	

# f(x) for hv dataset

getvalshv <- function(pathtofiles, viralseq, v){
	files <- list.files(pathtofiles, pattern='beta.*cdr3', full.names=T) 
	vals <- unlist(lapply(1:length(files), function(y){
		dat <- read.table(files[y], header=F, sep=',')
		total <- sum(dat[,2])
		common <- sum(dat[dat[,1] %in% viralseq,2])
		common/total
	}))
	type <- rep('hv blood', length(vals))
	prop <- vals
	pathogen <- rep(v, length(vals))
	patientid <- paste('HV', seq(1:length(vals)), sep='')
	df <- data.frame(type, prop, pathogen, patientid)
	df
}
		
cmv_beta <- getvals(tcrchain='beta', viralseq=cmvcdr3, v='cmv')		
ebv_beta <- getvals(tcrchain='beta', viralseq=ebvcdr3, v='ebv')		
flu_beta <- getvals(tcrchain='beta', viralseq=flucdr3, v='flu')
hiv_beta <- getvals(tcrchain='beta', viralseq=hivcdr3, v='hiv')

pathtohv <- '/Volumes/BF_MI_1/cvid-analysis/hv_cdr3_chainchud'

dtp <- do.call(rbind, c(cmv_beta, ebv_beta, flu_beta, hiv_beta))

ggplot(dtp, aes(x=patientid, y=prop, colour=type)) + geom_point(aes(shape=type), size=2, alpha=.7) + geom_hline(yintercept=0.008, color='red') + scale_colour_manual(values=c("#0072B2", "#000000", "#CC79A7")) + facet_wrap(~pathogen, nrow=4, scales='free_y') + guides(colour = guide_legend(override.aes = list(alpha=1))) + xlab('') + ylab('proportion of pathogen-associated CDR3s') + theme_bw() + theme(text=element_text(size=12), axis.text.x=element_text(angle = 90, hjust = 0.5, size=12))

ggsave(paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_Viral-LTX-TNB.png', sep=''), dpi=600, width=10, height=8)

ggplot(dtp, aes(x=type, y=prop)) + geom_boxplot() + xlab('') + ylab('proportion') + theme_bw() + theme(text=element_text(size=18), axis.text.x=element_text(size=10)) + facet_wrap(~pathogen, scales='free', nrow=2)

ggsave(paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_Viral-LTX-TNB-boxplot.png', sep=''), dpi=600, width=8, height=6)

# include hiv

dtp <- do.call(rbind, hiv_beta)

ggplot(dtp, aes(x=patientid, y=prop, colour=type)) + geom_point(aes(shape=type), size=2, alpha=.7) + scale_colour_manual(values=c("#0072B2", "#000000", "#CC79A7")) + facet_wrap(~pathogen, nrow=3, scales='free_y') + guides(colour = guide_legend(override.aes = list(alpha=1))) + theme_bw() + theme(text=element_text(size=18), axis.text.x=element_text(angle = 90, hjust = 0.5, size=12))

# include hv
hv_cmv <- getvalshv(pathtohv, cmvcdr3, 'cmv')
hv_ebv <- getvalshv(pathtohv, ebvcdr3, 'ebv')
hv_flu <- getvalshv(pathtohv, flucdr3, 'flu')
hv_hiv <- getvalshv(pathtohv, hivcdr3, 'hiv')

dtp <- do.call(rbind, c(cmv_beta, ebv_beta, flu_beta, hiv_beta))
dtp2 <- rbind(dtp, hv_cmv, hv_ebv, hv_flu, hv_hiv)

ggplot(dtp2, aes(x=type, y=prop)) + geom_boxplot() + geom_hline(yintercept=0.008, color='red') + xlab('') + ylab('proportion') + theme_bw() + theme(text=element_text(size=18), axis.text.x=element_text(size=10, angle=45, vjust=.8, hjust=.8)) + facet_wrap(~pathogen, nrow=2, scales='free')

save(dtp2, file=paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_Viral-LTX-TNBHV-boxplot.RData', sep=''))

ggsave(paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_Viral-LTX-HV-TNB-boxplot.png', sep=''), dpi=600, width=8, height=6)

# significance

a <- dtp2[which(dtp2$pathogen == 'ebv' & dtp2$type == 'blood'),2]
b <- dtp2[which(dtp2$pathogen == 'ebv' & dtp2$type == 'non-tumour lung'),2]
wilcox.test(a, b)

a <- dtp2[which(dtp2$pathogen == 'flu' & dtp2$type == 'blood'),2]
b <- dtp2[which(dtp2$pathogen == 'flu' & dtp2$type == 'tumour'),2]
wilcox.test(a, b)$p.value

a <- dtp2[which(dtp2$pathogen == 'flu' & dtp2$type == 'non-tumour lung'),2]
b <- dtp2[which(dtp2$pathogen == 'flu' & dtp2$type == 'tumour'),2]
wilcox.test(a, b)$p.value

labs <- combn(c('blood', 'non-tumour lung', 'tumour', 'hv blood'), 2, simplify=F)

getSig <- function(dtp)
{
	unlist(lapply(c('cmv', 'ebv', 'flu', 'hiv'), function(y)
	{
	lapply(1:length(labs), function(x)
	{
		v1 <- dtp[which(dtp$pathogen == y & dtp$type == labs[[x]][1]),2]
		v2 <- dtp[which(dtp$pathogen == y & dtp$type == labs[[x]][2]),2]
		
		wilcox.test(v1, v2)$p.value
	})
	}))
}

pvals <- getSig(dtp2)
labs2 <- unlist(lapply(labs, function(x) paste(x, collapse='.')))
vars <- rep(labs2, 4)
a <- rep(c('cmv', 'ebv', 'flu', 'hiv'), each=6)
sig_df <- data.frame(pvals, vars, a)







