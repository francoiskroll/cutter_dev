#####################################################
# ~ cutter: notes about fitting of deletion lengths ~
#
#
# Francois Kroll 2025
# francois@kroll.be
#####################################################

### 09/12/2024
# currently only using a set of deletions at slc45a2 for test
# may try on other loci later
# *** notes v0
# note v0: was taking *unique* deletions before (using dplyr's distinct())
# either in the entire dataset or within samples
# I think this was a mistake, as we want to model the probabilities as the dataset is
# not (currently) study more general properties of repair
# e.g. if 1-bp deletions are very frequent, when we simulate a sample we should generate many 1-bp deletions
# so I do not think we should lower the frequency of 1-bp deletions because many (within a sample) may be from the same repair event
# this "error" should be random anyways so corrected with sufficient samples
# "error" being for example
# "we overestimate the frequency of 1-bp deletions because some samples happened to have many of those
# because the edit occurred very early, so then many cells ended up with that mutation, or there were just many amplicons with this mutation from chance"
# both sources of error are random, sometimes we will underestimate because mutation occurred late in development or those amplicons were just not sampled often, etc"
# having said this, it could be a good idea to normalise for coverage (depth)
# as in this case we do have access to this (random) error
# (being some samples have higher coverage than others, so give inflated counts),
# so we can actually do something about it

# we have column "freq" in mut which is what we need
# it is number of reads with this mutation out of sample's coverage
# [to explore deeper here: there may be a small error from reads which have multiple deletions, probably rare...]
# as it is now it is basically = mutation's count if total coverage was 1
# a bit counterintuitive
# will * 1000 and round
# so we have read count if total coverage as 1000

### 28/07/2025
# using Cas9MiSeqDB dataset v0

library(ggplot2)
library(MASS)

source('~/Dropbox/cutter/simulateDel.R') # for function sumFreqBins

# function strNthSplit ----------------------------------------------------

strNthSplit <- function(stri,
                        split,
                        nth) {
  
  # confirm we are given string(s)
  stri <- as.character(unlist(stri))
  
  as.character(sapply(strsplit(stri, split=split),
                      function(s) {
                        s[nth]
                      }))
}



# import ------------------------------------------------------------------

mut <- read.csv('~/Dropbox/Cas9MiSeqDB/mutcalls.csv')

## how many reads?
length(unique(mut$urid))

### only keep deletions
del <- mut %>%
  filter(type=='del')
# mut has 233 samples

# but 191 injected samples
mut %>%
  filter(grp!='ni') %>%
  distinct(sample) %>%
  tally()

# only keeping deletion reads removes 10 samples
length(unique(del$sample))

# see *** notes v0

# work here is for function simulateDel
# which simulates only deletions; it does not simulate an entire sample with some reference reads, some insertions, etc.
# we want to give it realistic probabilities for deletions of various lengths
# e.g. if deletions of 55 bp never happen, we should never generate them
# accordingly, what we are interested is frequency of deletion of various lengths out of all reads with deletions
# for example we want to have: "6-bp deletions represent 2% of all deletions"
# so that we know to generate 2% of 6-bp deletions


# delete reads with multiple deletions ------------------------------------
# example
del %>%
  filter(sample=='190805_E11_mitfa-ad')
# many reads have two deletions: del_120_130_11_TGATGCACATA_NA & del_113_114_2_CC_NA

# keeping it as it is creates total frequencies > 1.0
# I would like to keep counts as number of reads with at least one deletion
# as in the end, we will always simulate one deletion per read
# but then we need to make it one read / one deletion here
# only good solution I can think of is to delete them
# alternatives would be to draw at random which deletion we keep but this amounts to dividing the frequencies of each deletion by two
# or to sum their lengths but does not really make sense to generate deletions we never observed

# which reads have multiple deletions?
ndelperread <- del %>%
  group_by(urid) %>%
  tally(name='ndel') %>%
  filter(ndel>1) # directly remove all the 1s (usual case)

nrow(ndelperread) # Cas9DB v0: 3218 reads to delete
# out of
length(unique(del$urid)) # Cas9DB v0: 200,000; so 1.6%

del <- del %>%
  filter(! urid %in% ndelperread$urid)
# now all reads mentioned in del have just one deletion

# but ! we need to update the counts "nreads"
# re-calculate them
new_nreads <- del %>%
  group_by(sample, mutid) %>%
  tally(name='nreads')
# add them
# delete the old one
del$nreads <- NULL
del <- left_join(del, new_nreads, by=c('sample', 'mutid'))


# calculate deletion frequencies ------------------------------------------

# as explained above: within each sample, frequency of each deletion out of all reads with a deletion

# to obtain this, we calculate number of reads with deletions for each sample
# we can simply count unique read ids from del
deltal <- del %>%
  group_by(sample) %>%
  summarise(spl_nreadsdel=n_distinct(urid))

# we add those counts to del
del <- left_join(del, deltal, by='sample')
# check nreads_del is always lower than splcov
del[del$spl_nreadsdel > del$splcov,]
# yes, always the case

# spl_nreadsdel is now, for each sample, total number of reads with at least one deletion

# would be a good idea to delete samples which have a low number of reads with deletions
# so we do not extrapolate probabilities of deletion lengths from low sample size (number of reads with deletions)
# say min. 100
length( unique(del[which(del$spl_nreadsdel < 100), 'sample']) ) # Cas9DB v0: 26 samples to exclude

del <- del %>%
  filter(spl_nreadsdel >= 100)

length(unique(del$sample))
# Cas9DB v0: 155 samples left

# above should have automatically deleted all ni samples
# check that
unique(del$grp) # yes, only inj left

# have a look at current distribution
ggplot(del, aes(x=bp)) + 
  geom_density()

### now to get the frequency of each deletion, out of all deletions in this sample
# we can simply divide nreads / spl_nreadsdel
# nreads is number of reads (in that sample) with that deletion
del <- del %>%
  mutate(delfreq=nreads/spl_nreadsdel)


# simplify dataset --------------------------------------------------------
# we do not need to keep unique reads anymore
# we can simply have one row for each sample, each deletion

# check first that each unique deletion has all the same frequency
check <- del %>%
  group_by(sample, mutid) %>%
  summarise(nfreqs=n_distinct(delfreq), .groups='drop')
unique(check$nfreqs) # all 1, correct

# summarised del: delsum
delsum <- del %>%
  distinct(sample, mutid, .keep_all=TRUE)

# delete useless columns
delsum$edit_type <- NULL
delsum$expedit <- NULL
delsum$pestrand <- NULL
delsum$rhapos <- NULL
delsum$rhadist <- NULL
delsum$ref <- NULL
delsum$ali <- NULL
delsum$rid <- NULL
delsum$urid <- NULL


### check: within each sample, sum of frequencies should be 1.0
check <- delsum %>%
  group_by(sample) %>%
  summarise_at(vars(delfreq),
               list(
                 sum=~sum(.)
               ))
unique(check$sum)


# average deletion frequencies per locus ----------------------------------

# number of samples for each locus vary
# so there would be some biases if we look now at frequency of each deletion length
# e.g. locus X has 20 samples and often has a 6-bp deletion,
# we would conclude that 6-bp deletions are very frequent just because we have many locus X samples
# solution is to use the multiple samples for each locus as replicates
# and calculate, for each unique deletion, its average frequency across loci
# e.g. locus X has 3 samples; and we find deletion #70--#75 at frequencies 0.2, 0.1, 0.5
# then we say: at locus X, frequency of deletion #70--#75 is 0.27
# if deletion is not found in a sample, we should assign it frequency 0.0

# simple group-by / summarise_at works but will not assign frequency 0.0 to samples of that locus which have frequency 0

# first have to expand the dataset so that, within locus, each sample mentions all possible deletions
# including those not observed with frequencies = 0.0

# (with help from GPT)
# ! should update with version from Claude used in forkMMEJ_sumfreqs.R
# first create every option
# i.e. lists, for each locus, all possible deletions
alludel <- delsum %>%
  distinct(locus, mutid) %>%
  left_join(delsum %>% distinct(locus, sample), by = 'locus')

# then add those deletions, with frequencies = 0 when needed
delexp <- alludel %>%
  left_join(delsum, by = c('sample', 'locus', 'mutid')) %>%
  mutate(delfreq = replace_na(delfreq, 0))

# looks correct, e.g.
# slc45a2-aa; del_93_95_3_AGT_NA was observed in all 3 samples in delsum, and that is the same in delexp (frequencies = 0.22, 0.13, 0.14)
# slc45a2-aa; del_88_95_8_CCTCTAGT_NA was observed in only 1 sample (190805_A03_slc24a5-aa, frequency = 0.18), now we have 0.18, 0.0, 0.0
# too many columns are NA (e.g. sid could give the real ones) but may not be an issue

# now ready to calculate average frequency of each deletion at each locus
delav <- delexp %>%
  group_by(locus, mutid) %>%
  summarise_at(vars(delfreq),
               list(
                 avgfreq=~mean(.),
                 nspl= ~length(.)
               ),
               .groups='drop')
# re-create some columns
# I cannot understand how to avoid losing them above
# mainly about splitting mutid
delav <- delav %>%
  mutate(type=strNthSplit(mutid, '_', 1), .before='avgfreq') %>%
  mutate(start=strNthSplit(mutid, '_', 2), .before='avgfreq') %>%
  mutate(stop=strNthSplit(mutid, '_', 3), .before='avgfreq') %>%
  mutate(bp=strNthSplit(mutid, '_', 4), .before='avgfreq') %>%
  mutate(refseq=strNthSplit(mutid, '_', 5), .before='avgfreq') %>%
  mutate(altseq=strNthSplit(mutid, '_', 6), .before='avgfreq')
# some checks

# e.g. cd2ap-ab, del_56_78_23_TGGAGTATGAGTATGAAGCCCTC_NA was observed in only one sample frequency = 0.10
# there are 12 cd2ap-ab samples, so that is (0.10 + 0.0 + 0.0 + ...) / 12 = 0.008
# correct

# csnk1db-aa, del_93_94_2_CA_NA was observed in all 4 samples
# checked average frequency, correct


# average frequency by deletion length ------------------------------------

# we now repeat similar steps as above
# to obtain frequency of deletion of various lengths
# we should add frequencies = 0.0 for loci where a particular deletion length was never observed
# e.g. if a 74-bp deletion is observed at only one locus, we should add frequency of 74-bp deletion at every other locus = 0.0
# so that it impacts the overall frequency of 74-bp deletions (i.e. should be very rare if observed at only one locus)

# there are
length(unique(delav$bp)) # unique deletion lengths

# first create every option of locus & deletion length
# i.e. lists, for each locus, all possible deletion lengths
alldellen <- delav %>%
  distinct(locus, nspl) %>%
  crossing(bp=unique(delav$bp))

# then add those lengths, with frequencies = 0 when not observed for a given locus
dellen <- alldellen %>%
  left_join(delav, by = c('locus', 'bp', 'nspl')) %>%
  mutate(avgfreq=replace_na(avgfreq, 0))
# looks correct, e.g.
# cd2ap-ab, we observed many deletion lengths
# but never 14 or 15 bp, and we added those as 0

# another way to check,
# now each locus should have the same number of unique deletion lengths
tmp <- dellen %>%
  group_by(locus) %>%
  summarise(ndellen=n_distinct(bp))
unique(tmp$ndellen)
# yes, all the number of unique deletion lengths


# sum frequencies of each deletion length within locus --------------------

# e.g. cd2ap-ab, we observe multiple 1-bp deletions
# we should sum to have one overall frequency of 1-bp deletion at cd2ap-ab
dellens <- dellen %>%
  group_by(locus, bp) %>%
  summarise_at(vars(avgfreq),
               list(
                 sumfreq=~sum(.)
               ))
# checked for cd2ap-ab, 1-bp deletions, correct


# average frequencies of deletion lengths across loci ---------------------

# now drop locus labels
# and calculate average frequencies of each length
dellenf <- dellens %>%
  group_by(bp) %>%
  summarise_at(vars(sumfreq),
               list(
                 freq=~mean(.),
                 nloci= ~length(.)
               ),
               .groups='drop')
# nloci is always the same, which makes sense

# bp was character all along, which is fine
# now convert to integer
dellenf$bp <- as.integer(dellenf$bp)
dellenf <- dellenf[order(dellenf$bp),]

# complete with missing lengths
missingLens <- which(! 1:max(dellenf$bp) %in% dellenf$bp ) # gives missing lengths, e.g. 58-bp deletion was never observed

# add those as frequency = 0.0
missingDf <- data.frame(bp=missingLens,
                        freq=rep(0, length(missingLens)),
                        nloci=rep(unique(dellenf$nloci), length(missingLens)))

dellenf <- rbind(dellenf, missingDf)

# order again the dataframe
dellenf <- dellenf[order(dellenf$bp),]


# barplot -----------------------------------------------------------------

# using function sumFreqBins from simulateDel.R
dellen_his <- sumFreqBins(df=dellenf,
                          colVals='bp',
                          colFreqs='freq',
                          every=5,
                          lower=0, # ***
                          last=NA)
# *** 0 bp deletion is not possible
# but putting 1 gives upbounds at 6, 11, ...
# 0 is OK; it gives first bin as 6 values (0, 1, 2, 3, 4, 5) but only 5 can really be observed
# then next bins are 5 values (6, 7, 8, 9, 10), etc.

ggDelLen <- ggplot(dellen_his, aes(x=upbound, y=sumfreq)) +
  geom_col()
ggDelLen
ggsave(here('~/Dropbox/cutter/dev_plots/delfreq_Cas9DB.pdf'), ggDelLen, width=70, height=70, units='mm')


# fit distribution --------------------------------------------------------

# can try all options offered by fitdistr (package MASS)
# uses MLE fitting
# https://www.rdocumentation.org/packages/MASS/versions/7.3-61/topics/fitdistr

# below, we need to fit counts, not frequencies directly
# so we create a sort of sample that represents all loci
# which we will say has 1M reads
# so that lowest frequencies (~ 1e-5) still give a positive read counts
dellenf$counts1M <- round(dellenf$freq * 1000000)
# because of the rounding we get total count 1000005, ignore

# ! some distributions do not want any counts = 0
# put them extremely low
# the smallest observed counts currently is n = 16 for 70-bp deletion (Cas9DB v0)
# will simply put 1
dellenf_nozero <- dellenf
dellenf_nozero[which(dellenf_nozero$counts1M==0), 'counts1M'] <- 1

# from this, we create a vector of 1M values
# e.g. 44,923 observations of 1-bp deletions
# so we write "1" 44,923 times
lens <- rep(dellenf$bp, times=dellenf$counts1M)
lens_nozero <- rep(dellenf_nozero$bp, times=dellenf_nozero$counts1M)

# check
freqBins(lens_nozero, every=1)
# yes correct

# "start" does not mean lower bound but some initial values for the parameters
# so that fitdistr can then iterate from them
# some distributions do accept a "lower" argument to set lower bound
# which we know here is 1 bp

### beta
# some help from GPT to set start arguments
# loc_start <- median(lens)
# scale_start <- IQR(lens) / 2
# fbet <- fitdistr(lens_nozero, 'beta', start = list(shape1=loc_start, shape2=scale_start))
# >> cannot make it fit, will skip

### cauchy
fcau <- fitdistr(lens, 'cauchy', lower=1)

### chi-squared
# similar than for beta, does not want any 0s
# using mean as start was suggested by GPT
start_df <- mean(lens) # I also tried median (as I have integers), very close result
fchi <- fitdistr(lens, 'chi-squared', start=list(df = start_df))

### exponential
fexp <- fitdistr(lens, 'exponential', lower=1) # does not require start, cf. documentation

### gamma
fgam <- fitdistr(lens, 'gamma', lower=1)

### geometric
fgeo <- fitdistr(lens, 'geometric', lower=1)

### lognormal
flog <- fitdistr(lens, 'lognormal', lower=1) # "log-normal" gives same result

### logistic
floi <- fitdistr(lens, 'logistic') # lower=1 makes it fail

### negative binomial
fbin <- fitdistr(lens, 'negative binomial', lower=1)

### normal
fnor <- fitdistr(lens, 'normal', lower=1) # does not require start, cf. documentation

### Poisson
fpoi <- fitdistr(lens, 'Poisson', lower=1) # does not require start, cf. documentation

### t
ft <- fitdistr(lens, 't', lower=1)

### weibull
fwei <- fitdistr(lens, 'weibull', lower=1)

# there are 13 distributions
# can proceed by elimination
# after fitting: 12 left


# assess fit with density on histogram ------------------------------------

# change frequencies so we just have

### cauchy
ggplot(dellenf, aes(x=bp, y=freq)) +
  geom_col() +
  stat_function(fun=dcauchy,
                args=list(location=fcau$estimate[1],
                          scale=fcau$estimate[2]),
                colour='red')
# pretty good

### exponential
ggplot(dellenf, aes(x=bp, y=freq)) +
  geom_col() +
  stat_function(fun=dexp,
                args=list(rate=fexp$estimate[1]),
                colour='red')
# >> I think no; does not catch the lower frequency of very short deletions (1 or 2 bp)
# which is probably real
# distribution should not start at the maximum

### gamma
ggplot(dellenf, aes(x=bp, y=freq)) +
  geom_col() +
  stat_function(
    fun=dgamma,
    args=list(shape=fgam$estimate[1],
              rate=fgam$estimate[2]),
    colour='red'
  )
# >> no; it puts frequency of 1-bp deletion at ~ 0 which is not right

### geometric
ggplot(dellenf, aes(x=bp, y=freq)) +
  geom_col() +
  stat_function(
    fun=dgeom,
    args=list(prob=fgeo$estimate[1]),
    colour='red'
  )
# >> no

### log-normal
ggplot(dellenf, aes(x=bp, y=freq)) +
  geom_col() +
  stat_function(
    fun=dlnorm,
    args=list(meanlog=flog$estimate[1],
              sdlog=flog$estimate[2]),
    colour='red'
  )
# very good!

### logistic
ggplot(dellenf, aes(x=bp, y=freq)) +
  geom_col() +
  stat_function(
    fun=dlogis,
    args=list(location=floi$estimate[1],
              scale=floi$estimate[2]),
    colour='red'
  )
# >> good but should go higher to catch the peak in the distribution


### negative binomial
ggplot(dellenf, aes(x=bp, y=freq)) +
  geom_col() +
  stat_function(
    fun=dnbinom,
    args=list(size=fbin$estimate[1],
              mu=fbin$estimate[2]),
    colour='red'
  )
# >> no; I think it is a discrete distribution


### normal
ggplot(dellenf, aes(x=bp, y=freq)) +
  geom_col() +
  stat_function(
    fun=dnorm,
    args=list(mean=fnor$estimate[1],
              sd=fnor$estimate[2]),
    colour='red'
  )
# >> ok but should go higher at the peak


### poisson
ggplot(dellenf, aes(x=bp, y=freq)) +
  geom_col() +
  stat_function(
    fun=dpois,
    args=list(lambda=fpoi$estimate[1]),
    colour='red'
  )
# >> no; discrete


### t
ggplot(dellenf, aes(x=bp, y=freq)) +
  geom_col() +
  stat_function(
    fun=dt,
    args=list(df=ft$estimate[3]),
    colour='red'
  )
# >> no, starts the maximum


### weibull
ggplot(dellenf, aes(x=bp, y=freq)) +
  geom_col() +
  stat_function(
    fun=dweibull,
    args=list(shape=fwei$estimate[1],
              scale=fwei$estimate[2]),
    colour='red'
  )
# >> not bad, but starts almost at the maximum, should start lower

### left in the race:
# cauchy
# log-normal
# logistic
# normal

# so 4 of them

# but really it is between cauchy and log-normal
# logistic and normal do not go high enough


# QQplot ------------------------------------------------------------------

### cauchy
draws <- rcauchy(n=10000, location=fcau$estimate[1], scale=fcau$estimate[2])

qqplot(draws, lens)
abline(0,1)
# weird...

### log-normal
draws <- rlnorm(n=10000, meanlog=flog$estimate[1], sdlog=flog$estimate[2])

qqplot(draws, lens)
abline(0,1)
# looks OKish

### logistic
draws <- rlogis(n=10000, location=floi$estimate[1], scale=floi$estimate[2])

qqplot(draws, lens)
abline(0,1)
# does not look great

### normal
draws <- rnorm(n=10000, mean=fnor$estimate[1], sd=fnor$estimate[2])

qqplot(draws, lens)
abline(0,1)
# not bad but it draws negative values...
# I think might exclude for this reason
# or should use half a truncated gaussian?
# truncate at 1 could make sense...

# I think will exclude cauchy,
# cannot make sense of qqplot


# KS test -----------------------------------------------------------------

# KS test is not happy because of ties in the data though...
# will ignore but may be an issue

### cauchy
ks.test(lens, 'pgamma', shape=fcau$estimate[1], rate=fcau$estimate[2]) # D = 0.79312 / p very low

### log-normal
ks.test(lens, 'plnorm', meanlog=flog$estimate[1], sdlog=flog$estimate[2]) # D = 0.07 / p very low

### logistic
ks.test(lens, 'plogis', location=floi$estimate[1], scale=floi$estimate[2]) # D = 0.174 / p very low

### normal
ks.test(lens, 'pnorm', mean=fnor$estimate[1], sd=fnor$estimate[2]) # D = 0.23 / p very low

# based on D distance statistic, log-normal is the best by far

# from quick reading wikipedia,
# log-normal looks widely used in biology

# >> log-normal is the one
# histogram was excellent


# where it the peak of log-normal? ----------------------------------------

ggplot(dellenf, aes(x=bp, y=freq)) +
  geom_col() +
  stat_function(
    fun=dlnorm,
    args=list(meanlog=flog$estimate[1],
              sdlog=flog$estimate[2]),
    colour='red'
  ) +
  geom_vline(xintercept=3)
# by trial and error, just after 3

# and mode of many values drawn at random
# (then round to simulate lengths of deletions)
# is 3, so same result
draws <- round(rlnorm(n=10000, meanlog=flog$estimate[1], sdlog=flog$estimate[2]))
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
getmode(draws)

# in real data:
getmode(lens) # also 3 bp is the most common

# I am sure someone better in statistics could do a better job but this seems a good approximation!


# summary -----------------------------------------------------------------

# we now have a fitted distribution of deletion lengths generated after Cas9 DSBs in zebrafish embryos
# we can use this distribution to simulate 'random' deletions at any locus


# polished plot for publication -------------------------------------------

ggdelfit <- ggplot(dellenf, aes(x=bp, y=freq)) +
  geom_col(fill='#8b8b8b') +
  stat_function(
    fun=dlnorm,
    args=list(meanlog=flog$estimate[1],
              sdlog=flog$estimate[2]),
    colour='#cb2a20',
    linewidth=0.5
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor.y=element_blank(),
    axis.title.x=element_text(size=9),
    axis.title.y=element_text(size=9),
    axis.text.y=element_text(size=7),
    axis.text.x=element_text(size=7)
  ) +
  xlab('deletion length (bp)') +
  ylab('frequency') +
  scale_x_continuous(breaks=c(1, 20, 40, 60))
ggdelfit
ggsave(here('~/Dropbox/cutter/dev_plots/delfreqCas9db.pdf'), ggdelfit, width=75, height=60, units='mm')
