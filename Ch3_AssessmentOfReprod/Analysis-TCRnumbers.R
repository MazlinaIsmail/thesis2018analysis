# Feb 2018
# plots showing TCR numbers for Paris samples
# get numbers from summary log files since it's easier to read in

rm(list=ls(all=T))

##### Load library #####
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(grid)

##### Fx blocks #####

fancy_scientific <- function(l){
	# https://stackoverflow.com/questions/11610377/how-do-i-change-the-formatting-of-numbers-on-an-axis-with-ggplot
	
	# turn in to character string in scientific notation
	l <- format(l, scientific=T)
	# quote the part before the exponent to keep all the digits
	l <- gsub('^(.*)e', "'\\1'e", l)
	# turn the 'e+' into platmath format
	l <- gsub('e', '%*%10^', l)
	# return this as an expression
	parse(text=l)
}

##### Start #####

##### (1) #####
# number of decombined reads
# alpha vs beta
miseqLogFile <- '/Volumes/BF_MI_3/Paris-sample-analysis/docs/Summary_MiSeq-Paris.csv'
miseqVals <- read.table(miseqLogFile, header=T, sep=',', stringsAsFactors=F)
miseqVals$chain <- sapply(miseqVals$sample, function(x){
	c <- strsplit(x, '_')[[1]][1]
	c
})
miseqVals$group <- sapply(miseqVals$sample, function(x){
	g <- strsplit(x, '_')[[1]][3]
	g
})
miseqVals$platform <- rep('MiSeq', nrow(miseqVals))

nextseqLogFile <- '/Volumes/BF_MI_1/NextSeq009/scripts/Summary_NextSeq009-Legion.csv'
nextseqVals <- read.table(nextseqLogFile, header=F, sep=',', stringsAsFactors=F)
nextseqValsParis <- nextseqVals[44:nrow(nextseqVals),]
nextseqValsParis$chain <- sapply(nextseqValsParis[,1], function(x){
	c <- strsplit(x, '_')[[1]][1]
	c
})
nextseqValsParis$group <- sapply(nextseqValsParis[,1], function(x){
	g <- strsplit(x, '_')[[1]][3]
	g
})
nextseqValsParis$platform <- rep('NextSeq', nrow(nextseqValsParis))

colnames(nextseqValsParis) <- colnames(miseqVals)

# Number of decombined reads = 3
# Total DCRs post-collapsing = 9
# Average RNA duplication = 14
# Unique DCRs post-collapsing = 8
type <- 9
ylabel <- 'Total DCRs post-collapsing'
dat <- rbind(miseqVals[,c(1,type,15,16,17)], nextseqValsParis[,c(1,type,15,16,17)])

p <- ggplot(dat, aes(x=dat[,3], y=as.numeric(dat[,2]), group=dat[,4], colour=dat[,5])) 
p + geom_point(size=2) + geom_line(alpha=.7, size=1.2) + scale_y_continuous(labels=fancy_scientific) + scale_color_manual(values=c('#000000', '#999999')) + theme(text=element_text(size=18))+ xlab('') + ylab(ylabel) + labs(colour='Platform') + theme_bw()

ggsave(filename=paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '-MiseqNextSeq-Paris-', gsub(' ', '', ylabel), '.png', sep=''), plot=last_plot(), width=4.5, height=4, units='in', dpi=600)

##### modified multiple plotting #####

dat <- rbind(miseqVals, nextseqValsParis)

ylabel <- c('number reads \n decombined', 'unique DCRs \n passing filters', 'total DCRs \n passing filters', 'unique DCRs \n post-collapsing', 'total DCRs \n post-collapsing')

xdat <- dat[,15]
gr <- dat[,16]
p1 <- ggplot(dat, aes(x=xdat, y=as.numeric(dat[,3]), group=gr, color=dat[,17])) + geom_point(size=1) + geom_line(alpha=.7, size=1) + scale_y_continuous(labels=fancy_scientific) + theme_bw() + theme(text=element_text(size=14), legend.position='none')+ xlab('') + ylab(ylabel[1]) + scale_color_manual(values=c('#000000', '#999999')) + labs(colour='Platform')

p2 <- ggplot(dat, aes(x=xdat, y=as.numeric(dat[,5]), group=gr, color=dat[,17])) + geom_point(size=1) + geom_line(alpha=.7, size=1) + scale_y_continuous(labels=fancy_scientific) + theme_bw() + theme(text=element_text(size=14))+ xlab('') + ylab(ylabel[2]) + scale_color_manual(values=c('#000000', '#999999')) + labs(colour='Platform')

p3 <- ggplot(dat, aes(x=xdat, y=as.numeric(dat[,6]), group=gr, color=dat[,17])) + geom_point(size=1) + geom_line(alpha=.7, size=1) + scale_y_continuous(labels=fancy_scientific) + theme_bw() + theme(text=element_text(size=14), legend.position='none')+ xlab('') + ylab(ylabel[3]) + scale_color_manual(values=c('#000000', '#999999')) + labs(colour='Platform')

p4 <- ggplot(dat, aes(x=xdat, y=as.numeric(dat[,8]), group=gr, color=dat[,17])) + geom_point(size=1) + geom_line(alpha=.7, size=1) + scale_y_continuous(labels=fancy_scientific) + theme_bw() + theme(text=element_text(size=14), legend.position='none')+ xlab('') + ylab(ylabel[4]) + scale_color_manual(values=c('#000000', '#999999')) + labs(colour='Platform')

p5 <- ggplot(dat, aes(x=xdat, y=as.numeric(dat[,9]), group=gr, color=dat[,17])) + geom_point(size=1) + geom_line(alpha=.7, size=1) + scale_y_continuous(labels=fancy_scientific) + theme_bw() + theme(text=element_text(size=14), legend.position='none')+ xlab('') + ylab(ylabel[5]) + scale_color_manual(values=c('#000000', '#999999')) + labs(colour='Platform')

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {

  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))

  grid.newpage()
  grid.draw(combined)

  # return gtable invisibly
  invisible(combined)

}

g <- grid_arrange_shared_legend(p1, p2, p3, p4, p5, ncol=2, nrow=3)
g <- arrangeGrob(p1, p2, p3, p4, p5)
ggsave(filename=paste('/Volumes/BF_MI_3/Paris-sample-analysis/', Sys.Date(), '-ParisDecombinatorNums.png', sep=''), g, dpi=600, device='png')

##### (2) #####

# create new column = chain_group_run
# remove leading zero in miseq group values
dat$cgr <- sapply(1:nrow(dat), function(x){
	gsub('0', '', paste(dat[x,3], dat[x,4], sep=''))
})

#dat2 <- rbind(dat[1:12,], dat[which(duplicated(dat[,6])),])
x1 <- dat[1:12,]
x_alpha <- x1[which(x1[,3]=='alpha'),]
x_beta <- x1[which(x1[,3]=='beta'),]

y1 <- dat[which(duplicated(dat[,6])),]
y_alpha <- y1[which(y1[,3]=='alpha'),]
y_beta <- y1[which(y1[,3]=='beta'),]

dtp_alpha <- data.frame(as.numeric(x_alpha[,2]), as.numeric(y_alpha[,2]), x_alpha$chain, x_alpha$group)
dtp_beta <- data.frame(as.numeric(x_beta[,2]), as.numeric(y_beta[,2]), x_beta$chain, x_beta$group)
colnames(dtp_alpha) <- c('x', 'y', 'chain', 'group')
colnames(dtp_beta) <- c('x', 'y', 'chain', 'group')

dtp <- rbind(dtp_alpha, dtp_beta)

ggplot(dtp, aes(x=dtp[,1], y=dtp[,2], shape=dtp[,3])) + geom_point(size=2) + scale_y_continuous(labels=function(x)x/1000) + scale_x_continuous(labels=function(x)x/1000) + geom_text_repel(aes(label=dtp[,4])) + xlab(expression('Total DCRs post-collapsing (MiSeq), x10'^3)) + ylab(expression('Total DCRs post-collapsing (NextSeq), x10'^3)) + theme_bw() + labs(shape='chain')

ggsave(filename=paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '-Paris-MiseqVsNextseq-Coll.png', sep=''), plot=last_plot(), dpi=600)

##############################################################################################

dat$dilution <- rep(c(1000, 100, 10), each=2)
dat$TRchain <- ifelse(dat$chain == 'alpha', 'TRA', 'TRB')

dat$mygroup <- rep('unique', nrow(dat))
dat$mygroup <- rep('total', nrow(dat))

dat1 <- dat
dat2 <- dat

colnames(dat1)[2] <- 'values'
colnames(dat2)[2] <- 'values'

dtp <- rbind(dat1, dat2)

#p <- ggplot(dat, aes(x=as.factor(dat[,7]), y=as.numeric(dat[,2]), color=dat[,8]))
#p + geom_point(aes(shape=dat[,5]), position=position_jitter(w=.2), size=4, stroke=2) + geom_text_repel(aes(label=dat[,4])) + scale_color_manual(values=c('#999999', '#D55E00'), name='chain') + theme(text=element_text(size=16)) + scale_shape_manual(values=c(4,16), name='platform') + ylab(ylabel) + xlab('dilution') + theme_bw() + scale_y_continuous(limits=c(7000, 72000), breaks=c(0, 10000, 20000, 30000, 40000, 50000, 60000, 70000)) + facet_wrap(~mygroup) 

p <- ggplot(dtp, aes(x=as.factor(dilution), y=as.numeric(values), color=TRchain))
p + geom_point(aes(shape=platform), position=position_jitter(w=.2), size=4, stroke=2) + geom_text_repel(aes(label=group)) + scale_color_manual(values=c('#999999', '#D55E00'), name='chain') + theme(text=element_text(size=16)) + scale_shape_manual(values=c(4,16), name='platform') + xlab('dilution') + theme_bw() + scale_y_continuous(limits=c(7000, 72000), breaks=c(0, 10000, 20000, 30000, 40000, 50000, 60000, 70000)) + facet_wrap(~mygroup) + ylab('number of decombined reads (DCRs)')

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_Paris_TotalAndUniq.png', sep=''), width=10, height=5)

##### (3) #####

dat <- rbind(miseqVals, nextseqValsParis)

dat$dilution <- rep(c(1000, 100, 10), each=2)

p <- ggplot(dat, aes(x=as.factor(dilution), y=as.numeric(UniqueDCRsPostCollapsing), color=platform, group=interaction(dilution, chain, platform)))
p + geom_boxplot(width=.75) + geom_point(position=position_jitterdodge()) + scale_color_manual(values=c('#999999', '#D55E00')) + theme(text=element_text(size=16)) + scale_shape_manual(values=c(4,16), name='platform') + xlab('Number of Jurkat mRNA (dilution)') + theme_bw() + theme(text=element_text(size=20)) + scale_y_continuous(breaks=c(0, 10000, 20000, 30000, 40000, 50000, 60000, 70000), labels=fancy_scientific) + ylab('Unique DCRs post-collapsing') + facet_wrap(~chain)

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_Paris_UniqDcrsPostColl.png', sep=''), dpi=600, width=10, height=5)

##### (4) number of UMIs #####

miseq <- list.files('/Volumes/BF_MI_3/Paris-sample-analysis/dcr_barcode_miseq', full.names=T)
nextseq <- list.files('/Volumes/BF_MI_3/Paris-sample-analysis/dcr_barcode_nextseq', full.names=T)

combfiles <- c(miseq, nextseq)
chain <- rep(c('alpha', 'beta'), 2, each=9)
df <- data.frame(chain)
df$dilution <- rep(c(1000, 100, 10))
df$platform <- rep(c('miseq', 'nextseq'), each=18)

df$umi <- unlist(lapply(1:length(combfiles), function(x) 
{
	dat <- read.table(combfiles[x], sep=',', header=F)
	nrow(dat)
}))

p <- ggplot(df, aes(x=as.factor(dilution), y=as.numeric(umi), color=platform, group=interaction(dilution, chain, platform)))
p + geom_boxplot(width=.75) + geom_point(position=position_jitterdodge()) + scale_color_manual(values=c('#999999', '#D55E00')) + scale_shape_manual(values=c(4,16), name='platform') + xlab('Number of Jurkat mRNA (dilution)') + theme_bw() + theme(text=element_text(size=20)) + scale_y_continuous(labels=fancy_scientific) + ylab('UMI') + facet_wrap(~chain)

ggsave(paste('/Volumes/BF_MI_1/figs-thesis-analysis/', Sys.Date(), '_Paris_UMI.png', sep=''), dpi=600, width=10, height=5)

##### (5) distribution of UMI usage #####

miseq <- list.files('/Volumes/BF_MI_3/Paris-sample-analysis/dcr_barcode_miseq', full.names=T)
nextseq <- list.files('/Volumes/BF_MI_3/Paris-sample-analysis/dcr_barcode_nextseq', full.names=T)

combfiles <- c(miseq, nextseq)
c <- rep(c('alpha', 'beta'), 2, each=9)
d <- rep(c('1-1000', '1-100', '1-10', '2-1000', '2-100', '2-10', '3-1000', '3-100', '3-10'), 4)
p <- rep(c('miseq', 'nextseq'), 1, each=18)

values <- lapply(1:length(combfiles), function(x) 
{
	dat <- read.table(combfiles[x], sep=',', header=F)
	df <- as.data.frame(dat[,2])
	df$chain <- rep(c[x])
	df$dilution <- rep(d[x])
	df$platform <- rep(p[x])
	df
})

dtp <- do.call(rbind, values)
colnames(dtp)[1] <- 'Freq'
cols <- c('chain', 'dilution', 'platform')
dtp[cols] <- lapply(dtp[cols], factor)

dtp_mseq_alpha <- dtp[which(dtp$chain == 'alpha' & dtp$platform == 'miseq'),]
dtp_nseq_alpha <- dtp[which(dtp$chain == 'alpha' & dtp$platform == 'nextseq'),]

brx_mseq <- pretty(range(dtp_mseq_alpha$Freq), n = nclass.Sturges(dtp_mseq_alpha$Freq),min.n = 1)
brx_nseq <- pretty(range(dtp_nseq_alpha$Freq), n = nclass.Sturges(dtp_nseq_alpha$Freq),min.n = 1)

ggplot() + geom_histogram(data=dtp_mseq_alpha, aes(x=Freq, colour=platform), fill='transparent', breaks=seq(1, 100, 1)) + geom_histogram(data=dtp_nseq_alpha, aes(x=Freq, colour=platform), fill='transparent', breaks=seq(1, 100, 1)) + scale_color_manual(values=c('#000000', '#D55E00')) + scale_x_continuous() + theme_bw() + theme(text=element_text(size=20)) + facet_wrap(~dilution, scales='free')

##### test area #####

a <- c(1, 2, 1, 1, 1, 1, 2, 3, 4, 5, 6)
df <- data.frame(a)

ggplot(df, aes(a)) + geom_histogram()

dtp <- data.frame(dtp_mseq_alpha[which(dtp_mseq_alpha$platform == 'miseq' & dtp_mseq_alpha$dilution == '1-1000'),1])

brx <- pretty(range(dtp$Freq), n = nclass.Sturges(dtp$Freq),min.n = 1)
brx <- seq(min(dtp$Freq), max(dtp$Freq), 1)

ggplot() + geom_histogram(data=dtp, aes(Freq), fill='grey', breaks=brx, bins=50) + theme_bw() + theme(text=element_text(size=20))

hist(dtp$Freq, breaks=brx, xlim=c(1, 100))

