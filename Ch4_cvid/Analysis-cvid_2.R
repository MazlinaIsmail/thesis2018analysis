# June 2018
# CVID analysis
# (a) diversity analysis subsetting by age (HV - cvid)
# (b) proportion of maits
# (c) mean insert length
# (d) mean deletion length
# (e) mean insert length - chain set only
# (f) mean cdr3 length
# (g) prop of expanded TCRs

rm(list=ls(all=T))

library(ggplot2)
library(reshape2)
library(viridis)
library(ineq)
library(entropy)

# cvid metadata
idmet <- read.table('/Volumes/BF_MI_2/TCRseq-sample-dets/metadata-tcrrepseq-cvid.txt', header=T, sep='\t')
# hv metadata
hvmet <- read.table('/Volumes/BF_MI_2/TCRseq-sample-dets/metadata-tcrrepseq-hv.txt', sep='\t', header=T)

# subset hv and cvid by age
range <- seq(min(hvmet$ageatsampling, na.rm=T), max(hvmet$ageatsampling, na.rm=T), 20)
hvmet$agecat <- findInterval(hvmet$ageatsampling, range)
idmet$agecat <- findInterval(idmet$ageatsampling, range)

pathtohv <- '/Volumes/BF_MI_1/cvid-analysis/hv_freq'
pathtochudhv <- '/Volumes/BF_MI_3/chudakov-dataset/data'
fileshv <- c(list.files(pathtohv, pattern='beta', full.names=T), list.files(pathtochudhv, pattern='beta.*freq', full.names=T))

pathtoid <- '/Volumes/BF_MI_1/cvid-analysis/cvid_freq'
filesid <- list.files(pathtoid, pattern='beta', full.names=T)

# read in files to get sample size, s
combfilenames <- c(fileshv, filesid)
comb_abund <- sapply(combfilenames, function(f){
	dat <- read.table(f, header=F, sep=',')
	freq <- sum(dat[,6])
	freq
})

s <- round((90/100)*min(comb_abund))
rm(comb_abund)

label <- c(paste('HV', seq(1:length(fileshv)), sep=''), paste('ID', seq(1:length(filesid)), sep=''))
group <- c(rep('HV', length(fileshv)), rep('ID', length(filesid)))
agecat <- c(hvmet$agecat, idmet$agecat)

##### (a) Gini/Shannon #####

# read in cdr3 and abundance 
pathtohv <- '/Volumes/BF_MI_1/cvid-analysis/hv_cdr3_chainchud'
pathtoid <- '/Volumes/BF_MI_1/cvid-analysis/cvid_cdr3_np'

# work out sample size, s
fileshv <- list.files(pathtohv, pattern=paste('beta.*cdr3', sep=''), full.names=T)
filesid <- list.files(pathtoid, pattern=paste('beta.*cdr3', sep=''), full.names=T)

combfilenames <- c(fileshv, filesid)
comb_abund <- sapply(combfilenames, function(f){
	dat <- read.table(f, header=F, sep=',')
	freq <- sum(dat[,2])
	freq
})

s <- round((90/100)*min(comb_abund))
rm(comb_abund)

# cvid metadata
idmet <- read.table('/Volumes/BF_MI_2/TCRseq-sample-dets/metadata-tcrrepseq-cvid.txt', header=T, sep='\t')
# hv metadata
hvmet <- read.table('/Volumes/BF_MI_2/TCRseq-sample-dets/metadata-tcrrepseq-hv.txt', sep='\t', header=T)

# subset hv and cvid by age
range <- seq(min(hvmet$ageatsampling, na.rm=T), max(hvmet$ageatsampling, na.rm=T), 20)
hvmet$agecat <- findInterval(hvmet$ageatsampling, range)
idmet$agecat <- findInterval(idmet$ageatsampling, range)

metrefname <- unlist(list(hvmet$refname, idmet$refnam))

new_order <- unlist(lapply(1:length(metrefname), function(x){
		file <- grep(metrefname[x], combfilenames, value=T)
}))

comb_div_avg <- NULL

for(i in 1:length(new_order)){
	dat <- read.table(new_order[i], header=F, sep=',')
	dat2 <- rep(dat[,1], dat[,2])
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

label <- c(paste('HV', seq(1:length(fileshv)), sep=''), paste('ID', seq(1:length(filesid)), sep=''))
group <- c(rep('HV', length(fileshv)), rep('ID', length(filesid)))
agecat <- c(hvmet$agecat, idmet$agecat)

rownames(comb_div_avg) <- label
comb <- as.data.frame(comb_div_avg)
comb$group <- group
comb$label <- label
comb$agecat <- as.factor(agecat)
colnames(comb)[1:2] <- c('Gini index', 'Shannon entropy')

dtp <- melt(comb)

n_fun <- function(x)
{
	return(data.frame(y=min(x-.15), label=paste('(', length(x), ')', sep='')))
}

ggplot(dtp, aes(x=agecat, y=value, colour=group)) + geom_boxplot(width=.75, outlier.shape=NA) + stat_summary(fun.data=n_fun, geom='text', position=position_dodge(width=.75)) + geom_point(alpha=.5, position=position_dodge(width=.75)) + scale_color_manual(values=c('#999999', '#E69F00')) + theme_bw() + theme(text=element_text(size=20)) + facet_wrap(~variable, scales='free')

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_GiniShanHvCvid-beta.png', sep=''), dpi=600, width=10, height=5)

##### (b) find MAITs #####

chain <- 'alpha'
fileshv <- list.files(pathtohv, pattern=chain, full.names=T)
pathtoid <- '/Volumes/BF_MI_1/cvid-analysis/cvid_freq'
filesid <- list.files(pathtoid, pattern=chain, full.names=T)

combfilenames <- list(fileshv, filesid)

getvals <- lapply(1:length(combfilenames), function(iter1)
{
	dat <- lapply(combfilenames[[iter1]], read.table, header=F, sep=',')
	pt_names <- unlist(lapply(combfilenames[[iter1]], function(f){strsplit(f, '/')[[1]][6]}))
	proportion <- unlist(lapply(1:length(dat), function(iter2)
	{
		t <- sum(dat[[iter2]][,6])
		# alpha maits
		#sum(dat[[iter2]][which(dat[[iter2]][,1] == 1 & dat[[iter2]][,2] == 22 | dat[[iter2]][,1] == 1 & dat[[iter2]][,2] == 2 | dat[[iter2]][,1] == 1 & dat[[iter2]][,2] == 9),6])/t
		# beta maits 
		sum(dat[[iter2]][which(dat[[iter2]][,1] >= 31 & dat[[iter2]][,1] <= 36 | dat[[iter2]][,1] == 15),6])/t
	}))
	df <- data.frame(proportion)
	df
})

dat <- do.call(rbind, getvals)
group <- c(rep('HV', length(fileshv)), rep('ID', length(filesid)))
dat$group <- group

# mann-whitney test comparing control and test
wilcox.test(dat[dat$group=='HV',1], dat[dat$group=='ID',1], correct=T)

n_fun <- function(x)
{
	return(data.frame(y=max(x+.01), label=paste('(', length(x), ')', sep='')))
}

# alpha
#ggplot(dat, aes(x=group, y=proportion, colour=group)) + geom_boxplot(outlier.shape=NA) + geom_jitter(height=0, size=3, alpha=.6) + stat_summary(fun.data=n_fun, geom='text', size=5) + scale_color_manual(values=c('#999999', '#E69F00')) + xlab('') + annotate('text', x=1.8, y=.085, label='p-value = <0.0005, Mann-Whitney test', size=3) + theme_bw() + theme(text=element_text(size=20))

# beta
ggplot(dat, aes(x=group, y=proportion, colour=group)) + geom_boxplot(outlier.shape=NA) + geom_jitter(height=0, size=3, alpha=.6) + stat_summary(fun.data=n_fun, geom='text', size=5) + scale_color_manual(values=c('#999999', '#E69F00')) + xlab('') + annotate('text', x=1.8, y=.085, label='p-value = 0.0003, Mann-Whitney test', size=3) + theme_bw() + theme(text=element_text(size=16))

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_MaitProp-beta.png', sep=''), dpi=600, width=5, height=5)

##### (c) mean insert length #####

# numbers for chudakov dataset
chudnums <- read.table(paste(pathtochudhv, '/basicstats.basicstats.txt', sep=''), header=T, sep='\t')

chuddat <- data.frame(rep('HV', nrow(chudnums)), rep('chudakov', nrow(chudnums)), chudnums$mean_insert_size)
colnames(chuddat) <- c('group', 'label', 'meaninsert')

pathtohv <- '/Volumes/BF_MI_1/cvid-analysis/hv_freq'
fileshv <- list.files(pathtohv, pattern='beta', full.names=T)
pathtoid <- '/Volumes/BF_MI_1/cvid-analysis/cvid_freq'
filesid <- list.files(pathtoid, pattern='beta', full.names=T)

combfilenames <- c(fileshv, filesid)
label <- c(paste('HV', seq(1:length(fileshv)), sep=''), paste('ID', seq(1:length(filesid)), sep=''))
group <- c(rep('HV', length(fileshv)), rep('ID', length(filesid)))

meaninsert <- unlist(lapply(1:length(combfilenames), function(x)
{
	dat <- read.table(combfilenames[x], header=F, sep=',')
	insertlen <- sapply(dat[5], function(y) nchar(as.character(y)))
	mean(insertlen)
}))

df <- data.frame(group, rep('chain', length(meaninsert)), meaninsert)
colnames(df)[2] <- 'label'
dtp <- rbind(df, chuddat)

ggplot(dtp, aes(x=group, y=meaninsert, color=label)) + geom_boxplot(width=.75) + scale_color_manual(values=c('#999999', '#D55E00')) + ylab('mean insert length') + theme_bw() + theme(text=element_text(size=18))

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_MeanInsertLen-CVIDvsHV-beta.png', sep=''), dpi=600, width=4, height=3)

##### (d) mean deletion length #####

pathtohv <- '/Volumes/BF_MI_1/cvid-analysis/hv_freq'
pathtoid <- '/Volumes/BF_MI_1/cvid-analysis/cvid_freq_notempus'

getvals <- function(field, c, g){
	
	fileshv <- list.files(pathtohv, pattern=c, full.names=T)
	filesid <- list.files(pathtoid, pattern=c, full.names=T)
	
	combfilenames <- c(fileshv, filesid)
	label <- c(paste('HV', seq(1:length(fileshv)), sep=''), paste('ID', seq(1:length(filesid)), sep=''))
	group <- c(rep('HV', length(fileshv)), rep('ID', length(filesid)))
	
	vals <- unlist(lapply(1:length(combfilenames), function(x){
	dat <- read.table(combfilenames[x], header=F, sep=',')
	dellen <- dat[,field]
	mean(dellen)
	}))
	chain <- rep(c, length(vals))
	gene <- rep(g, length(vals))
	dat <- data.frame(vals, label, group, chain, gene)
}

meandel_v_alpha <- getvals(field=3, c='alpha', g='V')
meandel_j_alpha <- getvals(field=4, c='alpha', g='J')
meandel_v_beta <- getvals(field=3, c='beta', g='V')
meandel_j_beta <- getvals(field=4, c='beta', g='J')

df <- rbind(meandel_v_alpha, meandel_j_alpha, meandel_v_beta, meandel_j_beta)
dtp <- melt(df)

n_fun <- function(x)
{
	return(data.frame(y=max(x+.4), label=paste('(', length(x), ')', sep='')))
}

ggplot(dtp, aes(x=gene, y=value, fill=group)) + geom_boxplot(width=.75, alpha=.7) + stat_summary(fun.data=n_fun, geom='text', size=4, position=position_dodge(width=.75)) + scale_fill_manual(values=c('#999999', '#D55E00')) + ylab('number of nucleotides\ndeleted (mean)') + theme_bw() + theme(text=element_text(size=18)) + facet_wrap(~chain)

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_MeanDelLen-CVIDvsHV.png', sep=''), dpi=600, width=6, height=3)

# stats test
alpha_hv_v <- df[which(df$chain == 'alpha' & df$gene == 'V' & df$group == 'HV'),1]
alpha_id_v <- df[which(df$chain == 'alpha' & df$gene == 'V' & df$group == 'ID'),1]

wilcox.test(alpha_hv_v, alpha_id_v)

alpha_hv_j <- df[which(df$chain == 'alpha' & df$gene == 'J' & df$group == 'HV'),1]
alpha_id_j <- df[which(df$chain == 'alpha' & df$gene == 'J' & df$group == 'ID'),1]

wilcox.test(alpha_hv_j, alpha_id_j)

beta_hv_v <- df[which(df$chain == 'beta' & df$gene == 'V' & df$group == 'HV'),1]
beta_id_v <- df[which(df$chain == 'beta' & df$gene == 'V' & df$group == 'ID'),1]

wilcox.test(beta_hv_v, beta_id_v)

beta_hv_j <- df[which(df$chain == 'beta' & df$gene == 'J' & df$group == 'HV'),1]
beta_id_j <- df[which(df$chain == 'beta' & df$gene == 'J' & df$group == 'ID'),1]

wilcox.test(beta_hv_j, beta_id_j)


##### (e) mean insert length - chain set only #####

pathtohv <- '/Volumes/BF_MI_1/cvid-analysis/hv_freq'
pathtoid <- '/Volumes/BF_MI_1/cvid-analysis/cvid_freq_notempus'

getvals <- function(c){
	
	fileshv <- list.files(pathtohv, pattern=c, full.names=T)
	filesid <- list.files(pathtoid, pattern=c, full.names=T)
	
	combfilenames <- c(fileshv, filesid)
	label <- c(paste('HV', seq(1:length(fileshv)), sep=''), paste('ID', seq(1:length(filesid)), sep=''))
	group <- c(rep('HV', length(fileshv)), rep('ID', length(filesid)))
	
	vals <- unlist(lapply(1:length(combfilenames), function(x){
	dat <- read.table(combfilenames[x], header=F, sep=',')
	insertlen <- sapply(dat[5], function(y) nchar(as.character(y)))
	mean(insertlen)
	}))
	chain <- rep(c, length(vals))
	dat <- data.frame(vals, label, group, chain)
}

mean_alpha <- getvals(c='alpha')
mean_beta <- getvals(c='beta')

df <- rbind(mean_alpha, mean_beta)
dtp <- melt(df)

n_fun <- function(x)
{
	return(data.frame(y=max(x+2), label=paste('(', length(x), ')', sep='')))
}

ggplot(dtp, aes(x=chain, y=value, fill=group)) + geom_boxplot(width=.75, alpha=.7) + stat_summary(fun.data=n_fun, geom='text', size=4, position=position_dodge(width=.75)) + scale_fill_manual(values=c('#999999', '#D55E00')) + ylab('number of nucleotides\ninserted (mean)') + theme_bw() + theme(text=element_text(size=18))

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_MeanInsLen-CVIDvsHV.png', sep=''), dpi=600, width=4, height=3)

# stats test
alpha_hv <- df[which(df$group == 'HV' & df$chain == 'alpha'),1]
alpha_id <- df[which(df$group == 'ID' & df$chain == 'alpha'),1]

wilcox.test(alpha_hv, alpha_id)

beta_hv <- df[which(df$group == 'HV' & df$chain == 'beta'),1]
beta_id <- df[which(df$group == 'ID' & df$chain == 'beta'),1]

wilcox.test(beta_hv, beta_id)

##### mean cdr3 length (split by tcr chain) #####

pathtochud <- '/Volumes/BF_MI_3/chudakov-dataset/data'
chudfiles <- list.files(pathtochud, pattern='txt.gz', full.names=T)
pathtohv <- '/Volumes/BF_MI_1/cvid-analysis/hv_cdr3_np'
pathtoid <- '/Volumes/BF_MI_1/cvid-analysis/cvid_cdr3_np_notempus'

# chain set
getvals <- function(c){
	
	fileshv <- list.files(pathtohv, pattern=paste(c, '.*cdr3', sep=''), full.names=T)
	filesid <- list.files(pathtoid, pattern=paste(c, '.*cdr3', sep=''), full.names=T)
	
	g <- c(rep('HV', length(fileshv)), rep('ID', length(filesid)))
	label <- c(paste('HV', seq(1:length(fileshv)), sep=''), paste('ID', seq(1:length(filesid)), sep=''))
	
	combfilenames <- c(fileshv, filesid)
	
	vals <- unlist(lapply(1:length(combfilenames), function(x){
	dat <- read.table(combfilenames[x], header=F, sep=',')
	cdr3len <- sapply(dat[1], function(y) nchar(as.character(y)))
	mean(cdr3len)
	}))
	source <- rep('chain', length(vals))
	chain <- rep(c, length(vals))
	group <- g
	dat <- data.frame(vals, source, chain, group, label)
}

chain_alpha <- getvals(c='alpha')
chain_beta <- getvals(c='beta')

# chudakov set
chud_vals <- unlist(lapply(1:length(chudfiles), function(x){
	dat <- read.table(chudfiles[x], header=T, sep='\t')
	cdr3len <- sapply(dat[4], function(y) nchar(as.character(y)))
	mean(cdr3len)
	}))
source <- rep('chudakov', length(chud_vals))
chain <- rep('beta', length(chud_vals))
group <- rep('HV', length(chud_vals))
label <- paste('HV', seq(30:108), sep='')
chud_dat <- data.frame(chud_vals, source, chain, group, label)

colnames(chud_dat) <- colnames(chain_alpha)

df <- rbind(chud_dat, chain_alpha, chain_beta)
dtp <- melt(df)

n_fun <- function(x)
{
	return(data.frame(y=max(x+.1), label=paste('(', length(x), ')', sep='')))
}

ggplot(dtp, aes(x=chain, y=value, fill=group)) + geom_boxplot(width=.75, alpha=.7) + stat_summary(fun.data=n_fun, geom='text', size=4, position=position_dodge(width=.75)) + scale_fill_manual(values=c('#999999', '#D55E00')) + ylab('mean CDR3 length') + theme_bw() + theme(text=element_text(size=18))

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_MeanCdr3Len-CVIDvsHV.png', sep=''), dpi=600, width=4, height=3)

# stats test
alpha_hv <- dtp[which(dtp$group == 'HV' & dtp$chain == 'alpha'),6]
alpha_id <- dtp[which(dtp$group == 'ID' & dtp$chain == 'alpha'),6]

wilcox.test(alpha_hv, alpha_id)

beta_hv <- dtp[which(dtp$group == 'HV' & dtp$chain == 'beta'),6]
beta_id <- dtp[which(dtp$group == 'ID' & dtp$chain == 'beta'),6]

wilcox.test(beta_hv, beta_id)

##########################################################
##### Global variables for the following code blocks #####
##########################################################

# cdr3
pathtohv <- '/Volumes/BF_MI_1/cvid-analysis/hv_cdr3_chainchud'
pathtoid <- '/Volumes/BF_MI_1/cvid-analysis/cvid_cdr3_np_notempus'

# cvid metadata
idmet <- read.table('/Volumes/BF_MI_2/TCRseq-sample-dets/metadata-tcrrepseq-cvid.txt', header=T, sep='\t')
# hv metadata
hvmet <- read.table('/Volumes/BF_MI_2/TCRseq-sample-dets/metadata-tcrrepseq-hv.txt', sep='\t', header=T)

# subset hv and cvid by age
range <- seq(min(hvmet$ageatsampling, na.rm=T), max(hvmet$ageatsampling, na.rm=T), 20)
hvmet$agecat <- findInterval(hvmet$ageatsampling, range)
idmet$agecat <- findInterval(idmet$ageatsampling, range)

##### (f) mean cdr3 length by age #####

getvals <- function(pathtofiles, c, g, metdata){
	
	files <- list.files(pathtofiles, pattern=paste(c, '.*cdr3', sep=''), full.names=T)
	new_order <- unlist(lapply(1:nrow(metdata), function(x){
		file <- grep(metdata[x,1], files, value=T)
	}))
	vals <- unlist(lapply(1:length(new_order), function(y){
		dat <- read.table(new_order[y], header=F, sep=',')
		cdr3len <- sapply(dat[1], function(y) nchar(as.character(y)))
		mean(cdr3len)
	}))
	chain <- rep(c, length(vals))
	group <- rep(g, length(vals))
	ageCat <- metdata$agecat
	dat <- data.frame(vals, chain, group, ageCat)
}

hv_chainchud <- getvals(pathtohv, c='beta', g='HV', metdata=hvmet)
id_chain <- getvals(pathtoid, c='beta', g='ID', metdata=idmet)

dtp <- rbind(hv_chainchud, id_chain)
dtp$ageCat <- as.factor(dtp$ageCat)

ggplot(dtp, aes(x=ageCat, y=vals, fill=group)) + geom_boxplot(width=.75, alpha=.7) + stat_summary(fun.data=n_fun, geom='text', size=4, position=position_dodge(width=.75)) + scale_fill_manual(values=c('#999999', '#D55E00')) + ylab('Mean CDR3 length') + theme_bw() + theme(text=element_text(size=18))

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_MeanCdr3Len-CVIDvsHV-byAge.png', sep=''), dpi=600, width=6, height=3)

##### (g) proportion of expanded TCRs #####

# read in cdr3 and abundance 
pathtohv <- '/Volumes/BF_MI_1/cvid-analysis/hv_cdr3_chainchud'
pathtoid <- '/Volumes/BF_MI_1/cvid-analysis/cvid_cdr3_np'

getvals <- function(pathtofiles, c, g, metdata){
	
	files <- list.files(pathtofiles, pattern=paste(c, '.*cdr3', sep=''), full.names=T)
	new_order <- unlist(lapply(1:nrow(metdata), function(x){
		file <- grep(metdata[x,1], files, value=T)
	}))
	vals <- unlist(lapply(1:length(new_order), function(y){
		dat <- read.table(new_order[y], header=F, sep=',')
		total <- sum(dat[,2])
		expanded <- sum(dat[which(dat[,2] >=8),2])
		expanded/total
	}))
	chain <- rep(c, length(vals))
	group <- rep(g, length(vals))
	ageCat <- metdata$agecat
	dat <- data.frame(vals, chain, group, ageCat)
}

# cvid metadata
idmet <- read.table('/Volumes/BF_MI_2/TCRseq-sample-dets/metadata-tcrrepseq-cvid.txt', header=T, sep='\t')
# hv metadata
hvmet <- read.table('/Volumes/BF_MI_2/TCRseq-sample-dets/metadata-tcrrepseq-hv.txt', sep='\t', header=T)

# subset hv and cvid by age
range <- seq(min(hvmet$ageatsampling, na.rm=T), max(hvmet$ageatsampling, na.rm=T), 20)
hvmet$agecat <- findInterval(hvmet$ageatsampling, range)
idmet$agecat <- findInterval(idmet$ageatsampling, range)

hv_chainchud <- getvals(pathtohv, c='beta', g='HV', metdata=hvmet)
id_chain <- getvals(pathtoid, c='beta', g='ID', metdata=idmet)

dtp <- rbind(hv_chainchud, id_chain)
dtp$ageCat <- as.factor(dtp$ageCat)

ggplot(dtp, aes(x=ageCat, y=vals, fill=group)) + geom_boxplot() + scale_fill_manual(values=c('#999999', '#D55E00')) + ylab('Fraction expanded (>8)') + theme_bw() + theme(text=element_text(size=18))

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_FractionExpanded-CVIDvsHV-byAge.png', sep=''), dpi=600, width=8, height=4)

x <- dtp[which(dtp$group ==  'HV' & dtp$ageCat == 4),1]
y <- dtp[which(dtp$group ==  'ID' & dtp$ageCat == 4),1]
wilcox.test(x, y)

#####################
##### Test area #####
#####################

fileshv <- list.files(pathtohv, pattern='beta.*cdr3', full.names=T)

# return file that matched order of metadata
lapply(1:nrow(hvmet), function(x){
	file <- grep(hvmet[x,1], fileshv, value=T)
	dat <- read.table(file, sep=',') 
	head(dat)
})

files <- list.files(pathtoid, pattern=paste('beta', '.*cdr3', sep=''), full.names=T)
new_order <- unlist(lapply(1:nrow(idmet), function(x){
	file <- grep(idmet[x,1], files, value=T, fixed=T)
}))



