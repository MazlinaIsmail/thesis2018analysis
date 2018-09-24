# July 2018
# renyi diversity [0, 1, Inf] for cvid vs hv (age-matched)

##### Pseudocode #####
# 1. read in metadata to get age matched individuals
# 2. calculate renyi values for the two cohorts
# 3. plot 

rm(list=ls(all=T))

library(vegan)
library(ggplot2)
library(reshape2)
library(ineq)
library(entropy)
library(gridExtra)
library(grid)
library(viridis)
library(ggfortify)

##### Start #####

# cvid metadata
idmet <- read.table('/Volumes/BF_MI_2/TCRseq-sample-dets/metadata-tcrrepseq-cvid.txt', header=T, sep='\t')
# hv metadata
hvmet <- read.table('/Volumes/BF_MI_2/TCRseq-sample-dets/metadata-tcrrepseq-hv.txt', sep='\t', header=T)
matchagehv <- hvmet[which(hvmet$ageatsampling >= 20 & hvmet$ageatsampling < 80),]

pathtohv <- '/Volumes/BF_MI_1/cvid-analysis/hv_cdr3_chainchud'
fileshv <- list.files(pathtohv, full.names=T)
matchhvfiles <- unlist(sapply(matchagehv$refname, function(x) grep(x, fileshv, value=T)))

pathtoid <- '/Volumes/BF_MI_1/cvid-analysis/cvid_cdr3_np'
filesid <- list.files(pathtoid, pattern='beta.*cdr3', full.names=T)

# read in files to get sample size, s
combfilenames <- c(matchhvfiles, filesid)
comb_abund <- sapply(combfilenames, function(f){
	dat <- read.table(f, header=F, sep=',')
	freq <- sum(dat[,2])
	freq
})

s <- round((90/100)*min(comb_abund))
s <- 13159 # beta
rm(comb_abund)

metrefname <- unlist(list(matchagehv$refname, idmet$refnam))

new_order <- unlist(lapply(1:length(metrefname), function(x){
		file <- grep(metrefname[x], combfilenames, value=T)
}))


label <- c(paste('HV', seq(1:length(matchhvfiles)), sep=''), paste('ID', seq(1:length(filesid)), sep=''))
group <- c(rep('HV', length(matchhvfiles)), rep('ID', length(filesid)))

comb_div_avg <- NULL

for(i in 1:length(new_order)){
	dat <- read.table(new_order[i], header=F, sep=',')
	dat2 <- rep(dat[,1], dat[,2])
	df <- NULL
	for(i in 1:100){
	subs_dat <- as.vector(table(as.vector(sample(dat2, s, replace=F))))
	H_vals <- renyi(subs_dat, scales=c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf), hill=F)
	df <- rbind(df, H_vals)
	}
subs_avg <- apply(df, 2, mean)
comb_div_avg <- rbind(comb_div_avg, subs_avg)
}

rownames(comb_div_avg) <- label
comb <- as.data.frame(comb_div_avg)
comb$group <- group
comb$label <- label

dtp <- melt(comb)

ggplot(dtp, aes(x=as.factor(variable), y=value, group=label, color=group)) + geom_point(size=2, alpha=.6) + geom_line(alpha=.6) + scale_color_manual(values=c('#999999', '#E69F00')) + ylab(expression(paste('H', alpha))) + xlab(expression(alpha)) + theme_bw() + theme(text=element_text(size=16))

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_RenyiHvCvid-betaCdr3.png', sep=''), dpi=600, width=5, height=4)

##### diversity with locid #####
# 1. calculate renyi 0, 1, Inf for hv and cvid
# 2. colour according to locid

# cvid metadata
idmet <- read.table('/Volumes/BF_MI_2/TCRseq-sample-dets/CVID_TCRrepseq-clinicalandphenotype.txt', header=T, sep='\t')
# hv metadata
hvmet <- read.table('/Volumes/BF_MI_2/TCRseq-sample-dets/metadata-tcrrepseq-hv.txt', sep='\t', header=T)

# subset hv and cvid by age
range <- seq(min(hvmet$ageatsampling, na.rm=T), max(hvmet$ageatsampling, na.rm=T), 20)
hvmet$agecat <- findInterval(hvmet$ageatsampling, range)
idmet$agecat <- findInterval(idmet$ageatsampling, range)

matchagehv <- hvmet[which(hvmet$ageatsampling >= 20 & hvmet$ageatsampling < 80),]

pathtohv <- '/Volumes/BF_MI_1/cvid-analysis/hv_cdr3_chainchud'
fileshv <- list.files(pathtohv, full.names=T)
matchhvfiles <- unlist(sapply(matchagehv$refname, function(x) grep(x, fileshv, value=T)))

pathtoid <- '/Volumes/BF_MI_1/cvid-analysis/cvid_cdr3_np'
filesid <- list.files(pathtoid, pattern='beta.*cdr3', full.names=T)
#filesid <- list.files(pathtoid, pattern='alpha.*cdr3', full.names=T)

combfilenames <- c(matchhvfiles, filesid)
# reorder files according to metadata
metrefname <- unlist(list(matchagehv$refname, idmet$refnam))
new_order <- unlist(lapply(1:length(metrefname), function(x){
		file <- grep(metrefname[x], combfilenames, value=T)
}))

s <- 13519

label <- c(paste('HV', seq(1:length(matchhvfiles)), sep=''), paste('ID', seq(1:length(filesid)), sep=''))
group <- c(rep('HV', length(matchhvfiles)), rep('ID', length(filesid)))
locid <- c(rep(1, length(matchhvfiles)), idmet$locid)
glild <- c(rep(1, length(matchhvfiles)), idmet$glild)
agecat <- c(matchagehv$agecat, idmet$agecat)

comb_div_avg <- NULL

for(i in 1:length(new_order)){
	dat <- read.table(new_order[i], header=F, sep=',')
	dat2 <- rep(dat[,1], dat[,2])
	df <- NULL
	for(i in 1:100){
	subs_dat <- as.vector(table(as.vector(sample(dat2, s, replace=F))))
	H_vals <- renyi(subs_dat, scales=c(0, 1, Inf), hill=F)
	df <- rbind(df, H_vals)
	}
subs_avg <- apply(df, 2, mean)
comb_div_avg <- rbind(comb_div_avg, subs_avg)
}

rownames(comb_div_avg) <- label
comb <- as.data.frame(comb_div_avg)
comb$group <- group
comb$label <- label
comb$locid <- as.factor(locid)
levels(comb$locid) <- c('not LOCID', 'LOCID')

dtp <- melt(comb, id=c('group', 'label', 'locid'))

n_fun <- function(x)
{
	return(data.frame(y=max(x+.5), label=paste('(', length(x), ')', sep='')))
}

ggplot(dtp, aes(x=group, y=value, color=locid)) + geom_boxplot(width=.75, outlier.shape=NA) + geom_point(position=position_jitterdodge(dodge.width=.75), size=3, alpha=.7) + stat_summary(fun.data=n_fun, geom='text', position=position_dodge(width=.75)) + scale_color_manual(values=c('#999999', '#D55E00')) + ylab('Renyi diversity') + theme_bw() + theme(text=element_text(size=20)) + facet_wrap(~variable)

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_renyisubset-cvidLocid-beta.png', sep=''), dpi=600, width=10, height=5)

##### Gini/Shannon #####

comb_div_avg <- NULL

for(i in 1:length(new_order)){
	dat <- read.table(new_order[i], header=F, sep=',')
	dat2 <- rep(dat[,1], dat[,2])
	df <- NULL
	for(i in 1:100){
	subs_dat <- as.vector(table(as.vector(sample(dat2, s, replace=F))))
	gini <- ineq(subs_dat, type='Gini')
	shannon <- entropy(subs_dat, unit='log')
	comb_vals <- cbind(gini, shannon)
	df <- rbind(df, comb_vals)
	}
subs_avg <- apply(df, 2, mean)
comb_div_avg <- rbind(comb_div_avg, subs_avg)
}

rownames(comb_div_avg) <- label
comb <- as.data.frame(comb_div_avg)
comb$group <- group
comb$label <- label
comb$locid <- as.factor(locid)
levels(comb$locid) <- c('not LOCID', 'LOCID')
comb$agecat <- as.factor(agecat)

n_fun <- function(x)
{
	return(data.frame(y=max(x+.2), label=paste('(', length(x), ')', sep='')))
}

p1 <- ggplot(comb, aes(x=agecat, y=gini, colour=group, shape=locid, group=interaction(group, locid, agecat))) + geom_boxplot(width=.75, outlier.shape=NA) + geom_point(position=position_jitterdodge(dodge.width=.75), size=3, alpha=.7) + stat_summary(fun.data=n_fun, geom='text', position=position_dodge(width=.75)) + scale_color_manual(values=c('#999999', '#D55E00')) + xlab('') + ylab('Gini index') + theme_bw() + theme(text=element_text(size=20))

p2 <-  ggplot(comb, aes(x=agecat, y=shannon, colour=group, shape=locid, group=interaction(group, locid, agecat))) + geom_boxplot(width=.75, outlier.shape=NA) + geom_point(position=position_jitterdodge(dodge.width=.75), size=3, alpha=.7) + stat_summary(fun.data=n_fun, geom='text', position=position_dodge(width=.75)) + scale_color_manual(values=c('#999999', '#D55E00')) + xlab('') + ylab('Shannon entropy') + theme_bw() + theme(text=element_text(size=20))

g <- grid_arrange_shared_legend(p1, p2)
ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_GiniShanHvCvidAgeCat-beta.png', sep=''), g, dpi=600, width=10, height=6)

##### diversity with phenotype data

new_order <- unlist(lapply(1:nrow(idmet), function(x){
		file <- grep(idmet$refnam[x], filesid, value=T)
}))

s <- 13519 # beta
#s <- 2514 # alpha

label <- paste('ID', seq(1:length(filesid)), sep='')
locid <- idmet$locid
perc_cd4memory <- idmet$perc_cd4memory
perc_cd4naive <- idmet$perc_cd4naive
perc_swmemorybc <- idmet$perc_swmemorybc
perc_cd21lobc <- idmet$perc_cd21lobc

comb_div_avg <- NULL

for(i in 1:length(new_order)){
	dat <- read.table(new_order[i], header=F, sep=',')
	dat2 <- rep(dat[,1], dat[,2])
	df <- NULL
	for(i in 1:100){
	subs_dat <- as.vector(table(as.vector(sample(dat2, s, replace=F))))
	gini <- ineq(subs_dat, type='Gini')
	shannon <- entropy(subs_dat, unit='log')
	comb_vals <- cbind(gini, shannon)
	df <- rbind(df, comb_vals)
	}
subs_avg <- apply(df, 2, mean)
comb_div_avg <- rbind(comb_div_avg, subs_avg)
}

rownames(comb_div_avg) <- label
comb <- as.data.frame(comb_div_avg)

df <- data.frame(comb, label, locid, perc_cd4memory, perc_cd4naive, perc_swmemorybc, perc_cd21lobc)

write.table(df, paste('/Volumes/BF_MI_1/cvid-analysis/scripts/', Sys.Date(), '_GiniShannonCvidVals-alpha.csv', sep=''), quote=F, sep=',', col.names=T, row.names=F)

# gini
df <- data.frame(comb$gini, label, locid, perc_cd4memory, perc_cd4naive, perc_swmemorybc, perc_cd21lobc)
levels(df$locid) <- c('not LOCID', 'LOCID')
dtp <- melt(df, id=c('label', 'locid', 'comb.gini'))

ggplot(dtp, aes(x=comb.gini, y=value)) + geom_point(aes(shape=locid, colour=locid), size=3, alpha=.7) + geom_smooth(method='lm') + xlab('Gini index') + scale_color_manual(values=c('#999999', '#D55E00')) + scale_y_continuous(limits=c(NA, NA)) + ylab('%') + theme_bw() + theme(text=element_text(size=20)) + facet_wrap(~variable, scales='free_y')

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_GiniCvidScatter-alpha2.png', sep=''), dpi=600, width=10, height=8)

# shannon
df <- data.frame(comb$shannon, label, locid, perc_cd4memory, perc_cd4naive, perc_swmemorybc, perc_cd21lobc)
levels(df$locid) <- c('not LOCID', 'LOCID')
dtp <- melt(df, id=c('label', 'locid', 'comb.shannon'))

ggplot(dtp, aes(x=comb.shannon, y=value)) + geom_point(aes(shape=locid, colour=locid), size=3, alpha=.7) + geom_smooth(method='lm') + xlab('Shannon entropy') + scale_color_manual(values=c('#999999', '#D55E00')) + scale_y_continuous(limits=c(NA, NA)) + ylab('%') + theme_bw() + theme(text=element_text(size=20)) + facet_wrap(~variable, scales='free_y')

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_ShannonCvidScatter-alpha2.png', sep=''), dpi=600, width=10, height=8)

##### pca using counts #####

infile <- '/Volumes/BF_MI_2/TCRseq-sample-dets/CVID_TCRrepseq-clinicalandphenotype.txt'
dat <- read.table(infile, header=T, sep='\t')

df <- dat[seq(8, ncol(dat), 1)]
autoplot(prcomp(na.omit(df), scale=T), data=na.omit(dat[seq(4, ncol(dat), 1)]), colour='locid', frame=T)

##### using glild #####

comb_div_avg <- NULL

for(i in 1:length(new_order)){
	dat <- read.table(new_order[i], header=F, sep=',')
	dat2 <- rep(dat[,1], dat[,2])
	df <- NULL
	for(i in 1:100){
	subs_dat <- as.vector(table(as.vector(sample(dat2, s, replace=F))))
	gini <- ineq(subs_dat, type='Gini')
	shannon <- entropy(subs_dat, unit='log')
	comb_vals <- cbind(gini, shannon)
	df <- rbind(df, comb_vals)
	}
subs_avg <- apply(df, 2, mean)
comb_div_avg <- rbind(comb_div_avg, subs_avg)
}

rownames(comb_div_avg) <- label
comb <- as.data.frame(comb_div_avg)
comb$group <- group
comb$label <- label
comb$glild <- as.factor(glild)
comb$glild <- addNA(comb$glild)
levels(comb$glild) <- c('no GLILD', 'GLILD', 'NA')
comb$agecat <- as.factor(agecat)

n_fun <- function(x)
{
	return(data.frame(y=max(x+.2), label=paste('(', length(x), ')', sep='')))
}

p1 <- ggplot(comb, aes(x=agecat, y=gini, colour=group, shape=glild, group=interaction(group, glild, agecat))) + geom_boxplot(width=.75, outlier.shape=NA) + geom_point(position=position_jitterdodge(dodge.width=.75), size=3, alpha=.7) + stat_summary(fun.data=n_fun, geom='text', position=position_dodge(width=.75)) + scale_color_manual(values=c('#999999', '#D55E00')) + xlab('') + ylab('Gini index') + theme_bw() + theme(text=element_text(size=20))

p2 <-  ggplot(comb, aes(x=agecat, y=shannon, colour=group, shape=glild, group=interaction(group, glild, agecat))) + geom_boxplot(width=.75, outlier.shape=NA) + geom_point(position=position_jitterdodge(dodge.width=.75), size=3, alpha=.7) + stat_summary(fun.data=n_fun, geom='text', position=position_dodge(width=.75)) + scale_color_manual(values=c('#999999', '#D55E00')) + xlab('') + ylab('Shannon entropy') + theme_bw() + theme(text=element_text(size=20))

g <- grid_arrange_shared_legend(p1, p2)
ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_GiniShanHvCvidAgeCat-beta.png', sep=''), g, dpi=600, width=10, height=6)

##### based on phenotype data #####
# switched memory B cell < 0.4
# cd21lo bc < 20

dat <- read.table('/Volumes/BF_MI_1/cvid-analysis/scripts/2018-07-20_GiniShannonCvidVals-beta.csv', header=T, sep=',')

dat$bcCat <- ifelse(dat$perc_swmemorybc < 0.4, ifelse(dat$perc_cd21lobc < 20, 'low', 'not low'), 'not low')

dat2 <- dat[1:36,]

# shannon
df <- subset(dat2, select=-c(gini, locid))
df$bcCat <- factor(df$bcCat)
#levels(df$bcCat) <- c('not LOCID', 'LOCID')
dtp <- melt(df, id=c('label', 'bcCat', 'shannon'))

ggplot(dtp, aes(x=shannon, y=value)) + geom_point(aes(colour=bcCat), size=3, alpha=.7) + geom_smooth(method='lm') + xlab('Shannon entropy') + scale_color_manual(values=c('#0072B2', '#D55E00')) + scale_y_continuous(limits=c(NA, NA)) + ylab('%') + theme_bw() + theme(text=element_text(size=20)) + facet_wrap(~variable, scales='free_y')

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_ShannonCvidScatterLowNotLow-beta.png', sep=''), dpi=600, width=10, height=8)

idmet$bcCat <- ifelse(dat$perc_swmemorybc < 0.4, ifelse(dat$perc_cd21lobc < 20, 'low', 'not low'), 'not low')

rownames(comb_div_avg) <- label
comb <- as.data.frame(comb_div_avg)
comb$group <- group
comb$label <- label
bcCat <- c(rep(1, length(matchhvfiles)), idmet$bcCat)
bcCat[bcCat == '1'] <- NA
bcCat[is.na(bcCat)] <- 'missing'
comb$bcCat <- bcCat

n_fun <- function(x)
{
	return(data.frame(y=max(x+.2), label=paste('(', length(x), ')', sep='')))
}
  
dtp <- comb[1:89,]
  
ggplot(dtp, aes(x=group, y=shannon, colour=bcCat, fill=bcCat)) + geom_point(position=position_jitterdodge(dodge.width=.9), size=3, alpha=.7) + geom_boxplot(outlier.shape=NA, alpha=.4, position=position_dodge(width=.9)) + stat_summary(fun.data=n_fun, geom='text', position=position_dodge(width=.9)) + scale_colour_manual(values=c('#0072B2', '#999999', '#D55E00')) + scale_fill_manual(values=c('#0072B2', '#999999', '#D55E00')) + xlab('') + ylab('Shannon entropy') + theme_bw() + theme(text=element_text(size=20))

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_GiniShanHvCvidBcCat-beta.png', sep=''), dpi=600, width=6, height=6)

##### CVID only #####

# (1) Gini/Shannon plot
# (2) Gini/Shannon against phenotype data

dat <- read.table('/Volumes/BF_MI_1/cvid-analysis/scripts/2018-07-20_GiniShannonCvidVals-beta.csv', header=T, sep=',')

dtp <- comb[1:89,]

n_fun <- function(x)
{
	return(data.frame(y=max(x+.2), label=paste('(', length(x), ')', sep='')))
}

p1 <- ggplot(dtp, aes(x=group, y=gini, colour=locid, group=interaction(group, locid))) + geom_boxplot(width=.75, outlier.shape=NA) + geom_point(position=position_jitterdodge(dodge.width=.75), size=3, alpha=.7) + stat_summary(fun.data=n_fun, geom='text', position=position_dodge(width=.75)) + scale_color_manual(values=c('#999999', '#D55E00')) + xlab('') + ylab('Gini index') + theme_bw() + theme(text=element_text(size=20))

p2 <-  ggplot(dtp, aes(x=group, y=shannon, colour=locid, group=interaction(group, locid))) + geom_boxplot(width=.75, outlier.shape=NA) + geom_point(position=position_jitterdodge(dodge.width=.75), size=3, alpha=.7) + stat_summary(fun.data=n_fun, geom='text', position=position_dodge(width=.75)) + scale_color_manual(values=c('#999999', '#D55E00')) + xlab('') + ylab('Shannon entropy') + theme_bw() + theme(text=element_text(size=20))

g <- grid_arrange_shared_legend(p1, p2)
ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_GiniShanHvCvidLocid-beta.png', sep=''), g, dpi=600, width=8, height=4)

write.table(comb, '/Volumes/BF_MI_1/cvid-analysis/scripts/2018-07-24_GiniShannonHvCvidVals-beta.csv', quote=F, sep=',', col.names=T, row.names=F)

# correlation test
dat <- read.table('/Volumes/BF_MI_1/cvid-analysis/scripts/2018-07-20_GiniShannonCvidVals-beta.csv', header=T, sep=',')

dat2 <- dat[1:36,]

cor(dat2$gini, dat2$perc_cd4memory, method='pearson', use='pairwise.complete.obs')
cor(dat2$gini, dat2$perc_cd4naive, method='pearson', use='pairwise.complete.obs')
cor(dat2$shannon, dat2$perc_cd4memory, method='pearson', use='pairwise.complete.obs')
cor(dat2$shannon, dat2$perc_cd4naive, method='pearson', use='pairwise.complete.obs')

cor(dat2$shannon, dat2$perc_swmemorybc, method='pearson', use='pairwise.complete.obs')
cor(dat2$shannon, dat2$perc_cd21lobc, method='pearson', use='pairwise.complete.obs')

# age match freq files
pathtohv <- '/Volumes/BF_MI_1/cvid-analysis/hv_freq_chainchud'
fileshv <- list.files(pathtohv, full.names=T)
matchhvfiles <- unlist(sapply(matchagehv$refname, function(x) grep(x, fileshv, value=T)))

file.copy(matchhvfiles, '/Volumes/BF_MI_1/cvid-analysis/hv_freq_chainchud_agematch')

