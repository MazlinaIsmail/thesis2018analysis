# Aug 2018
# (1) total tcr numbers - justify sample size
# (2) renyi diversity profile for T, NTL, B - solid line = median
# (3) save values in text file
# (4) gini/shannon 

library(vegan)
library(ineq)
library(entropy)
library(ggplot2)
library(reshape2)
library(corrplot)
library(pals)

##### Global variables #####

parentdir <- '/Volumes/BF_MI_3/cancer-analysis_2/TRACERx_TCR_up_to_NextSeq017'
ltxdirs <- grep('LTX', list.dirs(parentdir, full.names=T), value=T)

# files used for manuscript
fileId <- read.table('/Volumes/BF_MI_3/cancer-analysis_2/TracerX_TCR_paper/TRACERx_paper1/contributor-Data/2018-08-13_useforpaper.csv', header=T, sep=',')

# subset for n = 37
to_remove <- c('LTX019', 'LTX093', 'LTX012', 'LTX021', 'LTX031', 'LTX033', 'LTX046', 'LTX073', 'LTX175', 'LTX251', 'LTXL11', 'LTX163:R3', 'LTX110:R2')

fileId2 <- fileId[grep(paste(to_remove, collapse='|'), fileId[,1], invert=T),]
write.table(fileId2, paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_useforpaper.csv', sep=''), quote=F, sep=',', row.names=F, col.names=T)

fileLst <- NULL
for(i in 1:length(fileId2$Freq_File_Name))
{
	allfiles <- list.files(ltxdirs, pattern='freq', full.names=T, recursive=T)
	match <- grep(fileId2$Freq_File_Name[i], allfiles, value=T)
	fileLst[i] <- match
}

ltxid_filename <- data.frame(fileId2$LTX_ID, fileLst)
colnames(ltxid_filename) <- c('ltxid', 'filename')

ltx <- unique(unlist(lapply(ltxid_filename$ltxid, function(x) strsplit(as.character(x), ':')[[1]][1])))
# ltx N = 39
#load('/Volumes/BF_MI_3/cancer-analysis_2/TracerX_TCR_paper/TRACERx_paper1/contributor-Data/2018-08-25_39ltxiduseforpaper.RData')
#ltx <- gsub('LTX0', 'LTX', ltx)

##### (1) tcr numbers #####

values <- lapply(ltx, function(y)
{
	subset <- ltxid_filename[which(grepl(y, ltxid_filename$ltxid) == T),]
	dat <- subset[grep('PBMCF', subset$ltxid, invert=T),]
	
	dcr <- lapply(as.character(dat$filename), read.table, header=F, sep=',')
	
	total <- unlist(lapply(1:length(dcr), function(x){sum(dcr[[x]][,6])}))
	
	df <- data.frame(dat$ltxid, total)
	df$type <- ifelse(grepl(':R', df[,1]) == T, 'T', ifelse(grepl(':N', df[,1]) == T, 'NTL', 'B'))
	df$chain <- unlist(lapply(df[,1], function(x) strsplit(as.character(x), ':')[[1]][3]))
	df
})

dtp <- do.call(rbind, values)
ggplot(dtp, aes(x=type, y=total, group=interaction(chain, type), colour=chain)) + geom_boxplot() + scale_y_continuous(limits=c(NA, 10000))

save(dtp, file=paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_LTXtcrrepseq-tcrnums.RData', sep=''))

##### (2-3) renyi ######

# subset based on ltx id
#ltxid_filename <- do.call(rbind, lapply(ltx, function(x){ltxid_filename[grep(x, ltxid_filename$ltxid),]}))

s <- 2500

samplenames <- sapply(1:nrow(ltxid_filename), function(f)
{
	dat <- read.table(as.character(ltxid_filename[f,2]), header=F, sep=',')
	freq <- sum(dat[,6])
	idx <- ifelse(freq > s, f, NA)
	idx
})

row_idx <- samplenames[!is.na(samplenames)]

comb_renyi_avg <- NULL

for(i in row_idx){
	dat <- read.table(as.character(ltxid_filename[i,2]), header=F, sep=',')
	dat[,7] <- do.call(paste, c(dat[,1:5], sep=':'))
	dat2 <- rep(dat[,7], dat[,6])
	df <- NULL
	for(i in 1:100){
	subs_dat <- as.vector(table(as.vector(sample(dat2, s, replace=F))))
	H_vals <- renyi(subs_dat, scales=c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf), hill=F)
	df <- rbind(df, H_vals)
	}
subs_avg <- apply(df, 2, mean) # mean or median?
comb_renyi_avg <- rbind(comb_renyi_avg, subs_avg)
}

rownames(comb_renyi_avg) <- ltxid_filename[row_idx,1]
unique(unlist(lapply(rownames(comb_renyi_avg), function(x) strsplit(as.character(x), ':')[[1]][1])))

#write.table(comb_renyi_avg, paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_LTXtcrrepseq-renyi.csv', sep=''), quote=F, sep=',', row.names=T, col.names=NA)
save(comb_renyi_avg, file=paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_LTX-n37-tcrrepseq-renyi.RData', sep=''))

dtp <- data.frame(comb_renyi_avg)

dtp$chain <- unlist(lapply(rownames(dtp), function(x) strsplit(as.character(x), ':')[[1]][3]))
dtp$type <- ifelse(grepl(':R', rownames(dtp)) == T, 'T', ifelse(grepl(':N', rownames(dtp)) == T, 'NTL', 'B'))
dtp$ltxid <- rownames(dtp)

# another data frame for solid line (median)
solidlinevals <- function(c)
{
combvals <- lapply(c('T', 'NTL', 'B'), function(y)
{
	vals <- unlist(lapply(1:11, function(x)
	{
	median(dtp[which(dtp$chain == c & dtp$type == y),x])
	}))	
	df <- data.frame(vals)
	df$param <- seq(1, 11, 1)
	df$type <- rep(y)
	df$chain <- rep(c)
	df
})
}

solidalpha <- do.call(rbind, solidlinevals(c='alpha'))
solidbeta <- do.call(rbind, solidlinevals(c='beta'))

#dtp_solidline <- rbind(solidalpha, solidbeta)

dtp2 <- melt(dtp)

# reproduce BC's fig
dtp3 <- dtp2[which(dtp2$chain == 'alpha'),]

ggplot() + geom_line(data=dtp3, aes(x=variable, y=value, group=interaction(ltxid, type), colour=type), alpha=.1) + geom_line(data=solidbeta, aes(x=param, y=vals, colour=type), size=1) + geom_point(data=solidbeta, aes(x=param, y=vals, colour=type), size=2) + scale_colour_manual(values=c('blue', 'black', 'red')) + scale_x_discrete(labels=c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf)) + ylab('R\u{E9}nyi diversity') + xlab('R\u{E9}nyi scale') + labs(title='alpha') + theme_bw() + theme(text=element_text(size=20), axis.text.x=element_text(size=14))

# exclude blood
solidbeta2 <- solidbeta[which(grepl('T|NTL', solidbeta$type) == T),]
solidalpha2 <- solidalpha[which(grepl('T|NTL', solidalpha$type) == T),]
dtp4 <- dtp3[which(grepl('T|NTL', dtp3$type) == T),]

ggplot() + geom_line(data=dtp4, aes(x=variable, y=value, group=interaction(ltxid, type), colour=type), alpha=.1) + geom_line(data=solidalpha2, aes(x=param, y=vals, colour=type), size=1) + geom_point(data=solidalpha2, aes(x=param, y=vals, colour=type), size=2) + scale_colour_manual(values=c('darkblue', 'darkred')) + scale_x_discrete(labels=c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf)) + ylab('R\u{E9}nyi diversity') + xlab('R\u{E9}nyi scale') + labs(title='alpha') + theme_bw() + theme(text=element_text(size=20), axis.text.x=element_text(size=14))

ggsave(paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_LTX-n37-TNalpha-renyi.png', sep=''), dpi=600, width=6, height=5)

# colorblind palette
ggplot() + geom_line(data=dtp2, aes(x=variable, y=value, group=interaction(ltxid, type), colour=type), alpha=.1) + geom_line(data=dtp_solidline, aes(x=param, y=vals, colour=type), size=1) + geom_point(data=dtp_solidline, aes(x=param, y=vals, colour=type), size=2) + scale_colour_manual(values=c('#0072B2', 'black', '#CC79A7')) + scale_x_discrete(labels=c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf)) + ylab('Diversity') + xlab('R\u{E9}nyi scale') + theme_bw() + theme(text=element_text(size=20), axis.text.x=element_text(size=10)) + facet_wrap(~chain, scales='free')

ggsave(paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_LTX-TNB-renyi.png', sep=''), dpi=600, width=10, height=5)

# reproduce BC's boxplot fig
dtp5 <- dtp4[which(dtp4$variable == 'X0' | dtp4$variable == 'X1' | dtp4$variable == 'Inf.'),]
dtp5$type <- as.factor(dtp5$type)
dtp5$type <- factor(dtp5$type, levels=c('T', 'NTL'))

n_fun <- function(x)
{
	return(data.frame(y=max(x+.5), label=paste('(', length(x), ')', sep='')))
}

ggplot(dtp5, aes(x=type, y=value)) + geom_boxplot() + stat_summary(fun.data=n_fun, geom='text', position=position_dodge(width=.75)) + theme_bw() + theme(text=element_text(size=20), axis.text.x=element_text(size=20)) + xlab('') + ylab('R\u{E9}nyi diversity') +  facet_wrap(~variable)

ggplot(dtp5, aes(x=type, y=value)) + geom_boxplot() + theme_bw() + theme(text=element_text(size=20), axis.text.x=element_text(size=20)) + xlab('') + ylab('R\u{E9}nyi diversity') +  facet_wrap(~variable)

ggsave(paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_LTX-n37-TNalpha-renyisubset.png', sep=''), dpi=600, width=6, height=4)

# significance
labs <- combn(c('T', 'B', 'NTL'), 2, simplify=F)
labs <- combn(c('T', 'NTL'), 2, simplify=F)

getSig <- function(c, dtp)
{
	unlist(lapply(c('X0', 'X1', 'Inf.'), function(y)
	{
	lapply(1:length(labs), function(x)
	{
		v1 <- dtp[which(dtp$chain == c & dtp$type == labs[[x]][1] & dtp$variable == y),5]
		v2 <- dtp[which(dtp$chain == c & dtp$type == labs[[x]][2] & dtp$variable == y),5]
		
		wilcox.test(v1, v2)$p.value
	})
	}))
}

pvals <- getSig(c='alpha', dtp=dtp4)
labs2 <- unlist(lapply(labs, function(x) paste(x, collapse='.')))
vars <- rep(labs2, 3)
a <- rep(c(0, 1, Inf), 3)
sig_df <- data.frame(pvals, vars, a)

# rnaseq module score
renyisubset <- dtp[grep(':R', dtp$ltxid),c('X0', 'X1', 'Inf.', 'ltxid')]
renyisubset$label <- unlist(lapply(rownames(renyisubset), function(x) paste(strsplit(as.character(x), ':')[[1]][1:2], collapse=':')))
mod_scores <- read.table('/Volumes/BF_MI_3/cancer-analysis_2/add_datafiles/mod_scores.csv', header=T, sep=',')
mod_scores$label <- paste(mod_scores$PatientID, mod_scores$Region, sep=':')

# subset matching label between renyi table and mod scores
colnames(renyisubset)[ncol(renyisubset)] <- 'label' 
merged <- merge(renyisubset[,c('label', 'ltxid', 'X0', 'X1', 'Inf.')], mod_scores[,grep('PatientID|Region', names(mod_scores), invert=T)], by='label')
rownames(merged) <- merged$ltxid
merged2 <- merged[,grep('label|ltxid', names(merged), invert=T)]

corr_merged <- cor(merged2)
corr_test <- psych::corr.test(merged2, adjust='none')$p
# omit IL17
corr_merged2 <- corr_merged[which(rownames(corr_merged) == c('X0', 'X1', 'Inf.')),-which(colnames(corr_merged) == c('IL17'))]
corr_test2 <- corr_test[which(rownames(corr_test) == c('X0', 'X1', 'Inf.')),-which(colnames(corr_test) == c('IL17'))]

#new_order <- c('X0', 'X1', 'Inf.', 'T', 'CD4', 'CD8', 'B', 'mono', 'NK', 'Neut', 'TypeI_IFN', 'IFN_gamma', 'TNF_alpha')
new_order <- c('T', 'CD4', 'CD8', 'B', 'mono', 'NK', 'Neut', 'TypeI_IFN', 'IFN_gamma')
corr_merged3 <- corr_merged2[grep('Inf.', rownames(corr_merged2), invert=T),new_order]
corr_test3 <- corr_test2[grep('Inf.', rownames(corr_test2), invert=T),new_order]

colnames(corr_merged3)[5] <- c('Mono')
colnames(corr_test3)[5] <- c('Mono')
rownames(corr_merged3) <- c('Renyi(0)', 'Renyi(Inf)')
rownames(corr_test3) <- c('Renyi(0)', 'Renyi(Inf)')

png(paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_RenyiSubsetVals-LTX-TNB-modulescore.png', sep=''), units='in', width=5, height=3, res=600)
corrplot(corr_merged3, p.mat=corr_test3, insig='label_sig', tl.col='black', tl.srt=90, sig.level = c(0.001, 0.01, 0.05), pch.cex=1, pch.col='black', tl.cex=1, cl.cex=0.5)
dev.off()


##### (4) Gini/Shannon #####

s <- 2500

samplenames <- sapply(1:nrow(ltxid_filename), function(f)
{
	dat <- read.table(as.character(ltxid_filename[f,2]), header=F, sep=',')
	freq <- sum(dat[,6])
	idx <- ifelse(freq > s, f, NA)
	idx
})

row_idx <- samplenames[!is.na(samplenames)]

comb_div_avg <- NULL

for(i in row_idx){
	dat <- read.table(as.character(ltxid_filename[i,2]), header=F, sep=',')
	dat[,7] <- do.call(paste, c(dat[,1:5], sep=':'))
	dat2 <- rep(dat[,7], dat[,6])
	df <- NULL
	for(i in 1:100){
	subs_dat <- as.vector(table(as.vector(sample(dat2, s, replace=F))))
	gini <- ineq(subs_dat, type='Gini')
	shannon <- entropy(subs_dat, unit='log')
	comb_vals <- cbind(gini, shannon)
	df <- rbind(df, comb_vals)
	}
subs_avg <- apply(df, 2, mean) # mean or median?
comb_div_avg <- rbind(comb_div_avg, subs_avg)
}

rownames(comb_div_avg) <- ltxid_filename[row_idx,1]

write.table(comb_div_avg, paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_LTXtcrrepseq-ginishannon.csv', sep=''), quote=F, sep=',', row.names=T, col.names=NA)
save(comb_div_avg, file=paste('/Volumes/BF_MI_3/cancer-analysis_2/scripts2/', Sys.Date(), '_LTXtcrrepseq-ginishannon.RData', sep=''))

##### test area #####
