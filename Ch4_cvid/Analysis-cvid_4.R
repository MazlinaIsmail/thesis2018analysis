# July 2018

# cvid metadata
idmet <- read.table('/Volumes/BF_MI_2/TCRseq-sample-dets/CVID_TCRrepseq-clinicalandphenotype.txt', header=T, sep='\t')
# hv metadata
hvmet <- read.table('/Volumes/BF_MI_2/TCRseq-sample-dets/metadata-tcrrepseq-hv.txt', sep='\t', header=T)

range <- c(0, 20, 40, 60)
hvmet$agecat <- findInterval(hvmet$ageatsampling, range)
idmet$agecat <- findInterval(idmet$ageatsampling, range)

#matchagehv <- hvmet[which(hvmet$ageatsampling >= 20 & hvmet$ageatsampling < 80),]

pathtohv <- '/Volumes/BF_MI_1/cvid-analysis/hv_cdr3_chainchud'
fileshv <- list.files(pathtohv, full.names=T)
#matchhvfiles <- unlist(sapply(matchagehv$refname, function(x) grep(x, fileshv, value=T)))

pathtoid <- '/Volumes/BF_MI_1/cvid-analysis/cvid_cdr3_np_notempus'
filesid <- list.files(pathtoid, pattern='beta.*cdr3', full.names=T)

combfilenames <- c(fileshv, filesid)
metrefname <- unlist(list(hvmet$refname, idmet$refnam))
new_order <- unlist(lapply(1:length(metrefname), function(x){
		file <- grep(metrefname[x], combfilenames, value=T)
}))

s <- 13519 # beta

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

label <- c(paste('HV', seq(1:length(fileshv)), sep=''), paste('ID', seq(1:length(filesid)), sep=''))
group <- c(rep('HV', length(fileshv)), rep('ID', length(filesid)))
locid <- c(rep(1, length(fileshv)), idmet$locid)
glild <- c(rep(1, length(fileshv)), idmet$glild)
agecat <- c(hvmet$agecat, idmet$agecat)

rownames(comb_div_avg) <- label
comb <- as.data.frame(comb_div_avg)
comb$group <- group
comb$label <- label
comb$locid <- as.factor(locid)
levels(comb$locid) <- c('not LOCID', 'LOCID')
#comb$agecat <- as.factor(agecat)
comb$agecat <- addNA(as.factor(agecat))

save(comb, file=paste('/Volumes/BF_MI_1/cvid-analysis/scripts/', Sys.Date(), '_GiniShannonHvCvidVals-beta.RData', sep=''))

n_fun <- function(x)
{
	return(data.frame(y=max(x+.2), label=paste('(', length(x), ')', sep='')))
}

p1 <- ggplot(comb, aes(x=agecat, y=gini, colour=group, group=interaction(group, agecat))) + geom_boxplot(width=.75, outlier.shape=NA) + geom_point(position=position_jitterdodge(dodge.width=.75), size=3, alpha=.7) + stat_summary(fun.data=n_fun, geom='text', position=position_dodge(width=.75)) + scale_color_manual(values=c('#999999', '#D55E00')) + xlab('') + ylab('Gini index') + theme_bw() + theme(text=element_text(size=20))

p2 <-  ggplot(comb, aes(x=agecat, y=shannon, colour=group, group=interaction(group, agecat))) + geom_boxplot(width=.75, outlier.shape=NA) + geom_point(position=position_jitterdodge(dodge.width=.75), size=3, alpha=.7) + stat_summary(fun.data=n_fun, geom='text', position=position_dodge(width=.75)) + scale_color_manual(values=c('#999999', '#D55E00')) + xlab('') + ylab('Shannon entropy') + theme_bw() + theme(text=element_text(size=20))

g <- grid_arrange_shared_legend(p1, p2)
ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_GiniShanHvCvidAgeCat-beta.png', sep=''), g, dpi=600, width=10, height=6)

##### Gini Shannon scatter #####

dat <- read.table('/Volumes/BF_MI_1/cvid-analysis/scripts/2018-07-20_GiniShannonCvidVals-beta.csv', header=T, sep=',')

dat2 <- dat[1:36,]

# scatter for memory/naive cd4

# gini
df <- subset(dat2, select=-c(shannon, perc_swmemorybc, perc_cd21lobc))
levels(df$locid) <- c('not LOCID', 'LOCID')
dtp <- melt(df, id=c('label', 'locid', 'gini'))

ggplot(dtp, aes(x=gini, y=value)) + geom_point(aes(shape=locid, colour=locid), size=4, alpha=.8) + xlab('Gini index') + scale_color_manual(values=c('#999999', '#D55E00')) + scale_y_continuous(limits=c(NA, NA)) + ylab('%') + theme_bw() + theme(text=element_text(size=20), panel.spacing=unit(2, 'lines')) + facet_wrap(~variable, scales='fixed')

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_GiniCvidScatter-beta.png', sep=''), dpi=600, width=10, height=4)

# shannon
df <- subset(dat2, select=-c(gini, perc_swmemorybc, perc_cd21lobc))
levels(df$locid) <- c('not LOCID', 'LOCID')
dtp <- melt(df, id=c('label', 'locid', 'shannon'))

ggplot(dtp, aes(x=shannon, y=value)) + geom_point(aes(shape=locid, colour=locid), size=4, alpha=.8) + xlab('Shannon entropy') + scale_color_manual(values=c('#999999', '#D55E00')) + scale_y_continuous(limits=c(NA, NA)) + ylab('%') + theme_bw() + theme(text=element_text(size=20), panel.spacing=unit(2, 'lines')) + facet_wrap(~variable, scales='fixed')

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_ShannonCvidScatter-beta.png', sep=''), dpi=600, width=10, height=4)

# scatter for switched memory bc and cd21 lo

dat2$bcCat <- ifelse(dat2$perc_swmemorybc < 0.4, ifelse(dat2$perc_cd21lobc < 20, 'low', 'not low'), 'not low')
df <- subset(dat2, select=-c(gini, perc_cd4memory, perc_cd4naive, locid))

dtp <- melt(df, id=c('label', 'bcCat', 'shannon'))

ggplot(dtp, aes(x=shannon, y=value)) + geom_point(aes(colour=bcCat), size=4, alpha=.8) + xlab('Shannon entropy') + scale_color_manual(values=c('#0072B2', '#D55E00')) + scale_y_continuous(limits=c(NA, NA)) + ylab('%') + theme_bw() + theme(text=element_text(size=20), panel.spacing=unit(2, 'lines')) + facet_wrap(~variable, scales='fixed')

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_ShannonCvidBcCatScatter-beta.png', sep=''), dpi=600, width=10, height=4)


