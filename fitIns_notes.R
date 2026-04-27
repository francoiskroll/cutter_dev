#####################################################
# ~ cutter: fitting insertions, with the goal of simulating them ~
#
#
# Francois Kroll 2025
# francois@kroll.be
#####################################################

# approach would be
# detect newly synthesised sequences from the Cas9DB dataset, with start positions
# ! should not add extendNewlySeq
# two things to fit:
# start position from the nick
# length of newly synthesised
# then we'll draw a position at random & draw a length at random
# and generate the sequence
# hmmm, not sure
# issue is generating reads with both the substitutions & the insertions will be overly complex
# I might not do it... not sure how to it in a nice way

# simulating complete reads like we do for the deletions seem complicated
# as insertions typically are a mix of substitutions and insertions
# rather, we can simply simulate "newly synthesised sequences"
# and for each detect whether they are templated
# this is what detectTemplatedIns does anyways,
# for the detection, it does not make any difference which part of the sequence is substitution or insertion

# 01/08/2025: detectTemplatedIns will also add columns newlyseq & newlyseq_bp
# which will make it easier to fit the lengths

# question arises with extendNewlySeq
# default is currently 3 bp
# meaning that the newlyseq is artificially 6 bp longer
# and we are guaranteed to always find a 3-bp LCS
# if we simulate random sequences, we will miss this aspect: the smallest LCS will be 0, not 3 bp
# trying to add 3 matching bp on either side of the random newlyseq
# would make us go back to deciding which part of the sequence is a substitution or an insertion so we can place it correctly in the read
# > too complicated
# ignoring the difference is a problem though
# in the real sample: an LCS of 3 bp is meaningless
# in the simulated sample: an LCS of 3 bp is meaningful
# but we cannot just add 3 bp to the LCS we find,
# two LCSs of 3 and 2 bp (frequent) is not equivalent (equal probabilities) to one LCS of 5 bp (not frequent)
# > I think solution is to add the 3-bp extensions, do not need to be perfectly realistic with the positions
# all we need for them to be in the search window so that the probabilities match
# (i.e. that probability of a 3-bp match is always 100%)
# in practice, we can just take the 3-bp on either side of the cut

# work here is for function simulateDel
# which simulates only deletions; it does not simulate an entire sample with some reference reads, some insertions, etc.
# we want to give it realistic probabilities for deletions of various lengths
# e.g. if deletions of 55 bp never happen, we should never generate them
# accordingly, what we are interested is frequency of deletion of various lengths out of all reads with deletions
# for example we want to have: "6-bp deletions represent 2% of all deletions"
# so that we know to generate 2% of 6-bp deletions


# packages & functions ----------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

source('~/Dropbox/cutter/allelesToMutations.R')
source('~/Dropbox/cutter/detectTemplatedIns.R')
source('~/Dropbox/cutter/simulateDel.R')


# import Cas9MiSeqDB ------------------------------------------------------
# mutcalls

mut <- read.csv('~/Dropbox/Cas9MiSeqDB/mutcalls.csv')

# just take one sample as test
# mut <-  mut %>%
#   filter(sample=='190805_H05_tbx16-ab')


# detect templated insertions on Cas9MiSeqDB ------------------------------

# 27/01/2026 extendNewlySeq=3 >> =0
# mutins <- detectTemplatedIns(mut=mut,
#                              allowStretchMatches=5,
#                              extendNewlySeq=0,
#                              searchWindowStarts=20,
#                              minLCSbp=5)
# > takes a while!
# write.csv(mutins, '~/Dropbox/cutter/dev/Cas9MiSeqDB_mutins.csv', row.names=FALSE)

### CHECKPOINT
mutins <- read.csv('~/Dropbox/cutter/dev/Cas9MiSeqDB_mutins.csv')

### keep only insertions
# other rows are all NA for the templated insertion info anyways
ins <- mutins %>%
  filter(type=='ins')


# number of insertion reads
length(unique(ins$urid))

# number of samples
length(unique(ins$sample))

# number of loci
length(unique(ins$locus))

# delete reads with multiple insertions -----------------------------------

# which reads have multiple insertions?
ninsperread <- ins %>%
  group_by(urid) %>%
  tally(name='nins') %>%
  filter(nins>1) # directly remove all the 1s (usual case)

# stop; we may not have to delete them if they always have the same templated ins detection
# below = number of unique newlyseq per urid
tmp <- ins %>%
  filter(urid %in% ninsperread$urid) %>%
  group_by(urid) %>%
  summarise(nnewseq=n_distinct(newlyseq))
unique(tmp$nnewseq) # all 1s
# i.e. when we have multiple insertions in a read,
# the templated ins detection always took both insertions into account as one mutation (one newly synthesised sequence)
# therefore, we do not need to delete all these reads
# we just need to delete duplicates

# there are
nrow(ninsperread) # 2466 unique reads with multiple insertions
sum(duplicated(ins$urid)) # 2471 duplicates; more than above because a few reads have 3 insertions

# check makes sense:
all( ninsperread$urid %in% ins[duplicated(ins$urid),'urid'] ) # all the reads with multiple insertions are duplicates in ins, makes sense
all( ins[duplicated(ins$urid),'urid'] %in% ninsperread$urid ) # all the reads that are duplicates in ins have multiple insertions, makes sense

# or in other words, they are the same reads:
identical( sort(unique(ins[duplicated(ins$urid),'urid'])) , sort(ninsperread$urid))

# simply removing duplicates from ins will keep one instance of each
nrow(ins) # 36,960
nrow(ins[!duplicated(ins$urid),]) # 34,489
# so removing 2,471, makes sense
ins <- ins[!duplicated(ins$urid),]

# now, each row of ins is one read / one templated ins detection

# number of insertion reads now
length(unique(ins$urid)) # did not change, as expected

### reset nreads
# we need to count again number of reads with each mutation
# before, nreads was simply number of reads with this mutation
# total does not necessarily give the total number of reads with any insertion
# e.g. 100 reads with two insertions each
# we count first insertion: 100 reads have it; second insertion: 100 reads have it
# total 200 reads
# in any case, here we want number of reads with each newlyseq
# ! by chance, it happens that we have the same newlyseq for two different insertions
# e.g. 190131_F10_cd2ap-ab CCAACAAACCAAATCCAAAG
# we will pool them at this stage
# not perfectly ideal because they do look like different mutations
# but I think detectTemplatedIns is guaranteed to reach the same conclusion anyways, as the newlyseq is the main input it takes
# so will just say now that they are the same mutation
# summary: we count, for each sample and each newlyseq, number of reads which have it
nreads_new <- ins %>%
  group_by(sample, newlyseq) %>%
  tally(name='nreads')
# delete old nreads
ins$nreads <- NULL
# add the new counts
ins <- left_join(ins, nreads_new, by=c('sample', 'newlyseq'))


# calculate templated ins frequencies -------------------------------------
# more precisely, newlyseq frequencies
# within each sample, we want to calculate the frequency of each newlyseq out of all reads with ins (i.e. with a newlyseq) in that sample

# to obtain this, we calculate number of reads with an ins / with a newlyseq for each sample
instal <- ins %>%
  group_by(sample) %>%
  summarise(spl_nreadsins=n_distinct(urid))

# we add those counts to ins
ins <- left_join(ins, instal, by='sample')
# check nreads_ins is always lower than splcov
ins[ins$spl_nreadsins > ins$splcov,]
# yes, always the case

# spl_nreadsins is now, for each sample, total number of reads with an insertion
# (may have been more than insertion, precise wording is: with a newlyseq)

# would be a good idea to delete samples which have a low number of reads with insertions
# so we do not extrapolate probabilities of newlyseq lengths from low sample size (number of reads with insertions)
# say min. 100
length( unique(ins[which(ins$spl_nreadsins < 50), 'sample']) ) 
# Cas9DB v0:
# 68 samples to exclude if min. 100 reads, 44%
# 41 samples to exclude if min. 50 reads, 27%
# will do min. 50 reads, makes sense that we cannot be as stringent as for deletions as insertions are less frequent
# could make it a bit more stringent if add samples later

ins <- ins %>%
  filter(spl_nreadsins >= 50)

length(unique(ins$sample))
# Cas9DB v0: 113 samples left

length(unique(ins$locus))

# above should have automatically deleted all ni samples
# check that
unique(ins$grp) # yes, only inj left

# have a look at current distribution
ggplot(ins, aes(x=newlyseq_bp)) + 
  geom_density()
# looks like a similar distribution than deletions
# note: min. is 7 bp
# this actually represents a 1-bp insertion,
# but the newlyseq is extended by 3 bp on either side (so + 6 bp)
# >> 27/01/2026 now extendNewlySeq=0 so min. is 1 bp

### now to get the frequency of each deletion, out of all deletions in this sample
# we can simply divide nreads / spl_nreadsdel
# nreads is number of reads (in that sample) with that deletion
ins <- ins %>%
  mutate(insfreq=nreads/spl_nreadsins)


# simplify dataset --------------------------------------------------------
# we do not need to keep unique reads anymore
# we can simply have one row for each sample, each newlyseq

# summarised ins: inssum
inssum <- ins %>%
  distinct(sample, newlyseq, insfreq, .keep_all=TRUE)
# i.e. we just keep one row for each sample & newlyseq & frequency (to differentiate rare cases where the newlyseq is the same for two different insertions)

# delete useless columns
inssum$edit_type <- NULL
inssum$expedit <- NULL
inssum$pestrand <- NULL
inssum$rhapos <- NULL
inssum$rhadist <- NULL
inssum$ref <- NULL
inssum$ali <- NULL
inssum$rid <- NULL
inssum$urid <- NULL

### check: within each sample, sum of frequencies should be 1.0
check <- inssum %>%
  group_by(sample) %>%
  summarise_at(vars(insfreq),
               list(
                 sum=~sum(.)
               ))
unique( as.numeric(as.character(check$sum)) ) # very weird glitch where it gives multiple 1s, even after conversion
# but anyways, all 1s
# I think it is something about storage of numbers because dput() gives two values as 0.999999999999999


# average newlyseq frequencies per locus ----------------------------------

# number of samples for each locus vary
# so there would be some biases if we look now at frequency of each insertion length
# e.g. locus X has 20 samples and often has a CGGGATGAGC newlyseq (10 bp)
# we would conclude that 10-bp newlyseqs are very frequent just because we have many locus X samples
# solution is to use the multiple samples for each locus as replicates
# and calculate, for each unique newlyseq, its average frequency across loci
# e.g. locus X has 3 samples; and we find  a given newlyseq at frequencies 0.2, 0.1, 0.5
# then we say: at locus X, frequency of this newlyseq is is 0.27
# if a newlyseq is not found in a sample, we should assign it frequency 0.0

# simple group-by / summarise_at works but will not assign frequency 0.0 to samples of that locus which have frequency 0

# first have to expand the dataset so that, within locus, each sample mentions all possible newlyseqs
# including those not observed with frequencies = 0.0

# (with help from GPT)
# first create every option
# i.e. lists, for each locus, all possible newlyseqs
allunewseq <- inssum %>%
  distinct(locus, newlyseq) %>%
  left_join(inssum %>% distinct(locus, sample), by = 'locus')

# then add those deletions, with frequencies = 0 when needed
insexp <- allunewseq %>%
  left_join(inssum, by = c('sample', 'locus', 'newlyseq')) %>%
  mutate(insfreq = replace_na(insfreq, 0))

# looks correct, e.g.
# slc45a2-aa; newlyseq TAGGTGG was observed in all 3 samples, and that is the same in insexp (frequencies = 0.51, 0.14, 0.16)
# slc45a2-aa; newlyseq TCTGGCCACTCAGGG was observed in only 1 sample (190805_A03_slc24a5-aa, frequency = 0.14), now we have 0.14, 0.0, 0.0
# too many columns are NA (e.g. sid could give the real ones) but may not be an issue

# now ready to calculate average frequency of each newlyseq at each locus
insav <- insexp %>%
  group_by(locus, newlyseq) %>%
  summarise_at(vars(insfreq),
               list(
                 avgfreq=~mean(.),
                 nspl= ~length(.)
               ),
               .groups='drop')

# for a given locus & newlyseq, number of samples should always be the same
# as we added rows with frequencies = 0.0 above
check <- insav %>%
  group_by(locus, newlyseq) %>%
  summarise(n_nspl=n_distinct(nspl), .groups='drop')
unique(check$n_nspl) # 1, correct

# re-fill newlyseq_bp
# new rows we added were filled with NA
# ! there are some newlyseq which have hyphens in them
# we should not count those as they are not synthesised sequence
sum(str_detect(insav$newlyseq, '-')) # 100 of them, out of 1345
# count lengths of newlyseq, without hyphens
insav <- insav %>%
  mutate(newlyseq_bp=nchar(gsub('-', '', newlyseq)))

# e.g. cd2ap-ab, newlyseq ATGGAGTATGAGTACG was observed in only one sample frequency = 0.23
# there are 9 cd2ap-ab samples, so that is (0.23 + 0.0 + 0.0 + ...) / 9 = 0.025
# correct!


# average frequency by deletion length ------------------------------------

# we now repeat similar steps as above
# to obtain frequency of newlyseq of various lengths
# we should add frequencies = 0.0 for loci where a particular newlyseq length was never observed
# e.g. if a 74-bp newlyseq is observed at only one locus, we should add frequency of 74-bp newlyseq at every other locus = 0.0
# so that it impacts the overall frequency of 74-bp newlyseqs (i.e. should be very rare if observed at only one locus)

# there are
length(unique(insav$newlyseq_bp)) # unique newlyseq lengths

# first create every option of locus & newlyseq length
# i.e. lists, for each locus, all possible newlyseq lengths
allinslen <- insav %>%
  distinct(locus, nspl) %>%
  crossing(newlyseq_bp=unique(insav$newlyseq_bp))

# then add those lengths, with frequencies = 0 when not observed for a given locus
inslen <- allinslen %>%
  left_join(insav, by = c('locus', 'newlyseq_bp', 'nspl')) %>%
  mutate(avgfreq=replace_na(avgfreq, 0))
# looks correct, e.g.
# cd2ap-ab, we observed many insertion lengths
# but never a 6-bp one, so we added these as 0

# another way to check,
# now each locus should have the same number of unique deletion lengths
tmp <- inslen %>%
  group_by(locus) %>%
  summarise(ninslen=n_distinct(newlyseq_bp))
unique(tmp$ninslen)
# yes, all the number of unique deletion lengths


# sum frequencies of each deletion length within locus --------------------

# e.g. cd2ap-ab, we observe multiple 7-bp newlyseqs
# we should sum to have one overall frequency of 7-bp newlyseqs at cd2ap-ab
inslens <- inslen %>%
  group_by(locus, newlyseq_bp) %>%
  summarise_at(vars(avgfreq),
               list(
                 sumfreq=~sum(.)
               ))
# checked for cd2ap-ab, 7-bp newlyseqs, correct


# average frequencies of deletion lengths across loci ---------------------

# now drop locus labels
# and calculate average frequencies of each newlyseq length
inslenf <- inslens %>%
  group_by(newlyseq_bp) %>%
  summarise_at(vars(sumfreq),
               list(
                 freq=~mean(.),
                 nloci= ~length(.)
               ),
               .groups='drop')
# nloci is always the same, which makes sense

inslenf$newlyseq_bp <- as.integer(inslenf$newlyseq_bp)
inslenf <- inslenf[order(inslenf$newlyseq_bp),]

# complete with missing lengths
# in case we have some lengths between min and max which were never observed
# example
missingLens <- which(! min(inslenf$newlyseq_bp):max(inslenf$newlyseq_bp) %in% inslenf$newlyseq_bp ) # gives missing lengths, e.g. 58-bp deletion was never observed

# add those as frequency = 0.0
missingDf <- data.frame(newlyseq_bp=missingLens,
                        freq=rep(0, length(missingLens)),
                        nloci=rep(unique(inslenf$nloci), length(missingLens)))

inslenf <- rbind(inslenf, missingDf)

# order again the dataframe
inslenf <- inslenf[order(inslenf$newlyseq_bp),]


# barplot -----------------------------------------------------------------

# to fit below, I think should remove 6-bp lengths
# bit dodgy as not actually possible
# they are from newlyseq which are longer but have some hyphens -
# otherwise, looks good
# inslenf <- inslenf[-which(inslenf$newlyseq_bp==6),]
# >> 26/01/2026 no need to do this anymore
# min. newlyseq_bp is 1 bp which is logical

ggNewlyBp <- ggplot(inslenf, aes(x=newlyseq_bp, y=freq)) +
  geom_col()
ggNewlyBp

ggsave('~/Dropbox/cutter/dev_plots/newlyseqfreq_Cas9DB.pdf', ggNewlyBp, width=70, height=70, units='mm')


# fit distribution --------------------------------------------------------

# can try all options offered by fitdistr (package MASS)
# uses MLE fitting
# https://www.rdocumentation.org/packages/MASS/versions/7.3-61/topics/fitdistr

# below, we need to fit counts, not frequencies directly
# so we create a sort of sample that represents all loci
# which we will say has 1M reads
# so that lowest frequencies (~ 1e-5) still give a positive read counts
inslenf$counts1M <- round(inslenf$freq * 1000000)
# because of the rounding we get total count 1000005, ignore

# ! some distributions do not want any counts = 0
# put them extremely low
# the smallest observed counts currently is n = 11 for 38-bp newlyseq (Cas9DB v0)
# will simply put 1
inslenf_nozero <- inslenf
inslenf_nozero[which(inslenf_nozero$counts1M==0), 'counts1M'] <- 1

# from this, we create a vector of 1M values
# e.g. 97,627 observations of 7-bp insertions
# so we write "1" 97,627 times
lens <- rep(inslenf$newlyseq_bp, times=inslenf$counts1M)
lens_nozero <- rep(inslenf_nozero$newlyseq_bp, times=inslenf_nozero$counts1M)

# check
freqBins(lens_nozero, every=1)
# yes correct

# "start" does not mean lower bound but some initial values for the parameters
# so that fitdistr can then iterate from them
# some distributions do accept a "lower" argument to set lower bound
# which we know here is 7 bp

### beta
# did not work for deletion length, will skip

### cauchy
fcau <- fitdistr(lens, 'cauchy', lower=7)

### chi-squared
# similar than for beta, does not want any 0s
# using mean as start was suggested by GPT
start_df <- mean(lens) # I also tried median (as I have integers), very close result
fchi <- fitdistr(lens, 'chi-squared', start=list(df = start_df))

### exponential
fexp <- fitdistr(lens, 'exponential', lower=7) # does not require start, cf. documentation

### gamma
fgam <- fitdistr(lens, 'gamma', lower=7)

### geometric
fgeo <- fitdistr(lens, 'geometric', lower=7)

### lognormal
flog <- fitdistr(lens, 'lognormal', lower=7) # "log-normal" gives same result

### logistic
floi <- fitdistr(lens, 'logistic', lower=7)

### negative binomial
# fbin <- fitdistr(lens, 'negative binomial', lower=7)
# commented out because it throws error

### normal
fnor <- fitdistr(lens, 'normal', lower=7) # does not require start, cf. documentation

### Poisson
fpoi <- fitdistr(lens, 'Poisson', lower=7) # does not require start, cf. documentation

### t
ft <- fitdistr(lens, 't', lower=7)

### weibull
fwei <- fitdistr(lens, 'weibull', lower=7)

# there are 13 distributions
# can proceed by elimination
# after fitting: 12 left


# assess fit with density on histogram ------------------------------------

# change frequencies so we just have

### cauchy
ggplot(inslenf, aes(x=newlyseq_bp, y=freq)) +
  geom_col() +
  stat_function(fun=dcauchy,
                args=list(location=fcau$estimate[1],
                          scale=fcau$estimate[2]),
                colour='red')
# OK but I think no, in contrast to deletions, I think here should start at the maximum

### exponential
ggplot(inslenf, aes(x=newlyseq_bp, y=freq)) +
  geom_col() +
  stat_function(fun=dexp,
                args=list(rate=fexp$estimate[1]),
                colour='red')
# better, but fit is not great, could start higher
# >> 27/01/2026 solved! looks much better
# extendNewlySeq=3 was the issue

### gamma
ggplot(inslenf, aes(x=newlyseq_bp, y=freq)) +
  geom_col() +
  stat_function(
    fun=dgamma,
    args=list(shape=fgam$estimate[1],
              rate=fgam$estimate[2]),
    colour='red'
  )
# >> no!

### geometric
ggplot(inslenf, aes(x=newlyseq_bp, y=freq)) +
  geom_col() +
  stat_function(
    fun=dgeom,
    args=list(prob=fgeo$estimate[1]),
    colour='red'
  )
# >> no!

### log-normal
ggplot(inslenf, aes(x=newlyseq_bp, y=freq)) +
  geom_col() +
  stat_function(
    fun=dlnorm,
    args=list(meanlog=flog$estimate[1],
              sdlog=flog$estimate[2]),
    colour='red'
  )
# >> the fit is good but should start at the maximum

### logistic
ggplot(inslenf, aes(x=newlyseq_bp, y=freq)) +
  geom_col() +
  stat_function(
    fun=dlogis,
    args=list(location=floi$estimate[1],
              scale=floi$estimate[2]),
    colour='red'
  )
# >> not bad but does not start at the maximum


### negative binomial
ggplot(inslenf, aes(x=newlyseq_bp, y=freq)) +
  geom_col() +
  stat_function(
    fun=dnbinom,
    args=list(size=fbin$estimate[1],
              mu=fbin$estimate[2]),
    colour='red'
  )
# >> no!


### normal
ggplot(inslenf, aes(x=newlyseq_bp, y=freq)) +
  geom_col() +
  stat_function(
    fun=dnorm,
    args=list(mean=fnor$estimate[1],
              sd=fnor$estimate[2]),
    colour='red'
  )
# >> same again, fit is good but should start at the maximum


### poisson
ggplot(inslenf, aes(x=newlyseq_bp, y=freq)) +
  geom_col() +
  stat_function(
    fun=dpois,
    args=list(lambda=fpoi$estimate[1]),
    colour='red'
  )
# >> no; discrete


### t
ggplot(inslenf, aes(x=newlyseq_bp, y=freq)) +
  geom_col() +
  stat_function(
    fun=dt,
    args=list(df=ft$estimate[3]),
    colour='red'
  )
# >> no


### weibull
ggplot(inslenf, aes(x=newlyseq_bp, y=freq)) +
  geom_col() +
  stat_function(
    fun=dweibull,
    args=list(shape=fwei$estimate[1],
              scale=fwei$estimate[2]),
    colour='red'
  )
# >> no!

# >>> I think only exponential will start at the maximum
# fit is the same whether I remove lower=7 or not


# QQplot ------------------------------------------------------------------

draws <- rexp(n=10000, rate=fexp$estimate[1])

qqplot(draws, lens)
abline(0,1)


# KS test -----------------------------------------------------------------

# KS test is not happy because of ties in the data though...
# will ignore but may be an issue

ks.test(lens, 'pexp', rate=fexp$estimate[1]) # D = 0.37 / p very low
# 27/01/2026: D = 0.16 / p very low
# lower D = better!


# mode of the draws? ------------------------------------------------------

# mode of many values drawn at random
# (then round to simulate lengths of deletions)
# is 3, so same result
draws <- round(rexp(n=10000, rate=fexp$estimate[1]))

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
getmode(draws)
# 1 bp, but if we just exclude them, we would get 7 bp mode

# in real data:
getmode(lens)
# 7 bp
# 27/01/2026 also 1

# the exponential fit is not amazing
# but should be a decent approximation
# 27/01/2026: now much better!

# summary -----------------------------------------------------------------

# we now have a fitted distribution of newlyseq lengths generated after Cas9 DSBs in zebrafish embryos
# we can use this distribution to simulate 'random' insertions (newly synthesised sequences) at any locus



# polished plot for publication -------------------------------------------

ggdelfit <- ggplot(inslenf, aes(x=newlyseq_bp, y=freq)) +
  geom_col(fill='#8b8b8b') +
  stat_function(fun=dexp,
                args=list(rate=fexp$estimate[1]),
                colour='#cb2a20',
                linewidth=0.5) +
  theme_minimal() +
  theme(
    panel.grid.minor.y=element_blank(),
    axis.title.x=element_text(size=9),
    axis.title.y=element_text(size=9),
    axis.text.y=element_text(size=7),
    axis.text.x=element_text(size=7)
  ) +
  xlab('length of newly synthesised sequence (bp)') +
  ylab('frequency') +
  coord_cartesian(xlim=c(1, 98)) +
  scale_x_continuous(breaks=c(1, 25, 50, 75, 100))
ggdelfit

ggsave('~/Dropbox/cutter/dev_plots/newlyseqfreqCas9db.pdf', ggdelfit, width=75, height=60, units='mm')
