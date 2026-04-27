#####################################################
# ~ cutter: simulate random deletions ~
#
# see lots of notes in fitDelLengths_notes.R
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################
library(ggplot2)
library(dplyr)
library(MASS)
# v0

# expects mut, mutation table

# v1
# cf. fitDelLengths_notes on Cas9DB v0: log-normal parameters of deletion lengths are
# flog$estimate[1] = meanlog = 1.955957
# flog$estimate[2] = sdlog = 0.8970051


# simulateDel -------------------------------------------------------------

# 11/02/2026
# can now handle multiple loci

simulateDel <- function(mut,
                        nreads=1000,
                        cutDelbp=3,
                        awayfromCut=4,
                        fit_meanlog=1.955957,
                        fit_sdlog=0.8970051) {
  
  # loop through loci in mut
  loci <- unique(mut$locus)
  
  ### for each locus, simulate reads
  simdelL <- lapply(loci, function(li) {
    
    cat('\t \t \t \t Simulating deletion reads for locus', li, '\n')
    
    # mutation table for one locus
    mutli <- mut %>%
      filter(locus==li)
    
    # check we only have one locus to be safe
    if(length(unique(mutli$locus))>1) stop('\t \t \t \t Error simulateDel: can only simulate a sample for *one* locus at a time.\n')
    # get setting cutpos from the existing mut dataframe
    cutpos <- unique(mutli$cutpos)
    
    ### now simulate `nreads` with each a deletion
    # step below is surprisingly slow...
    simdelL <- lapply(1:nreads, function(i) {
      cat('\t \t \t \t simulating read', i, 'of', nreads, '\n')
      simulateDel_one(mut=mutli,
                      fit_meanlog=fit_meanlog,
                      fit_sdlog=fit_sdlog,
                      cutpos=cutpos,
                      cutDelbp=3,
                      awayfromCut=4)
    })
    simdel <- do.call(rbind, simdelL)
    
    # correct some column types
    simdel$start <- as.integer(simdel$start)
    simdel$stop <- as.integer(simdel$stop)
    simdel$bp <- as.integer(simdel$bp)
    
    ### add rid
    # normally done by alleleToMutation
    # here, will simply be r1.1, r2.1, r3.1, etc.
    # until e.g. r1000.1 if nreads=1000
    simdel$rid <- sprintf('r%i.1', 1:nrow(simdel))
    
    ### add coverage
    # copying from allelesToMutations so we get the same columns
    # calculate total number of reads for this sample, i.e. coverage
    # and add it as column
    # here, we know this is exactly `nreads`
    cat('\t \t \t \t >>> Total number of reads:', nreads, '\n')
    simdel <- simdel %>%
      mutate(splcov=nreads)
    
    ### run mutation table through filterMutations
    # with all filters off
    # just to get the same columns
    simdel <- filterMutations(muttb=simdel,
                              minnreads=NA,
                              cutpos=NA,
                              rhapos=NA,
                              rhadist=NA,
                              controltb=NA,
                              callSubs=TRUE)
    
    ### add rundate, locus, grp, well
    # these we know are present in mut
    # add correct locus name
    # grp: make it locus_simulated
    simdel <- simdel %>%
      mutate(well='sim', .before=1) %>%
      mutate(grp=paste(unique(mutli$locus), 'simulated', sep='_')) %>%
      mutate(locus=unique(mutli$locus), .before=1) %>%
      mutate(rundate=999999, .before=1)
    
    ### also create sample column
    # which is rundate_well_locus
    simdel <- simdel %>%
      mutate(sample=paste(rundate, well, locus, sep='_'), .before=1)
    
    ### add last columns as NA
    # find the ones we are missing from mut
    # these were added from meta
    cols2add <- colnames(mutli) [! colnames(mutli) %in% colnames(simdel)]
    
    for(colnm in cols2add) {
      simdel <- simdel %>%
        mutate(tmp=NA, .before='well')
      # now put the correct column name
      colnames(simdel)[which(colnames(simdel)=='tmp')] <- colnm
    }
    
    ### if cutpos column is present, fill it with actual info given by user
    if('cutpos' %in% colnames(simdel)) {
      simdel$cutpos <- cutpos
    }
    
    return(simdel)
    # simdel gets added to list simdelL
    
  }) # closes lapply
  # we now have one simdel dataframe per locus
  # stored in a list
  
  # rbind the list (we can keep track of locus through locus column)
  simdf <- do.call(rbind, simdelL)
  
  ### add all simulated reads to mut
  # first put the columns of simdf in the same order as in mut
  
  # check columns are identical (before looking at their positions)
  if(!identical( sort(colnames(mut)) , sort(colnames(simdf))))
    stop('\t \t \t \t >>> Error simulateDel: all columns of mut are not in simdf, or vice-versa.\n')
  
  # now put the columns of simdf in the same order
  simdf <- simdf[, colnames(mut) ]
  # check correct
  if(!identical(colnames(mut), colnames(simdf)))
    stop('\t \t \t \t >>> Error simulateDel: columns of simdf are not in the same order as columns of mut.\n')
  
  # can now rbind
  mutsim <- rbind(mut, simdf)
  # and we return result
  return(mutsim)
  
}



# simulateDel_one ---------------------------------------------------------

simulateDel_one <- function(mut,
                            fit_meanlog,
                            fit_sdlog,
                            cutpos,
                            cutDelbp,
                            awayfromCut) {
  
  # if delOK did not get switch to TRUE below
  # then simulate a new deletion
  delOK <- FALSE
  
  while(!delOK) {
    
    ### from the fitted log-normal distribution, we draw one deletion length
    # and round it to closest integer
    # ! if this integer is 0, we draw again
    len <- 0
    while(len==0) {
      len <- round(rlnorm(n=1,
                          meanlog=fit_meanlog,
                          sdlog=fit_sdlog))
    }
    
    ### recover original reference sequence
    # i.e. the reference sequence used for alignment
    # (in alignment, the reference sequence can get --- to indicate insertions)
    # get all the 'true' original reference sequences:
    # which we do by just removing the --- from the aligned reference seauences
    orefs <- gsub('-', '', mut$ref) # "original references"
    # check they are all the same, after we remove NA
    orefs <- orefs[!is.na(orefs)]
    
    if(length(unique(orefs))!=1)
      stop('\t \t \t \t >>> Error simulateDel_one: issue when computing overall reference sequence.\n')
    # we have the original reference sequence
    oref <- unique(orefs)
    
    ### 
    # deletion must remove the 'cutpos' nucleotide
    # which is PAM minus 4
    # cf. https://elifesciences.org/articles/59683#fig2s1, it is indeed the most frequently deleted nucleotide
    # except if 1 or 2 nt (flexible), then we can leave a bit of margin
    # but we still want the deletion to be very close, say ± 4 bp (flexible too)
    
    if(len < cutDelbp) {
      # if deletion is very small, 
      # we do not necessarily have to remove `cutpos` nucleotide
      # deletion most shifted to the left that is acceptable is
      # e.g. if cutpos is #87 and length is 2 bp
      # then most shifted to the left is #83, #84
      # because 4 bp (`awayfromCut`) from #84 includes #87, which is cutpos
      delstart <- cutpos - awayfromCut - len + 1 + 1 # a bit by trial and error to get what I want...
      delposleft <- delstart : (delstart+len-1)
      
      # then we just shift the deletion to the right 1 bp at a time
      # we do it until the left-most deleted position is `awayfromCut` away from cutpos
      # note, window of possible deletions is not symmetrical around `cutpos`
      # this is on purpose, because cutpos labels the nucleotide left of the cut
      # it *is* symmetrical if you consider the exact cut position (i.e. where the DNA chain is broken)
      shiftNtimes <- (cutpos + awayfromCut) - delstart
      
      # starting at 0 allows to get delposleft as first element of the list
      # then we shift one by one to the right
      delposL <- lapply(0:shiftNtimes, function(le) {
        delposleft + le
      })
      
    } else {
      # if deletion is fairly big, it must delete `cutpos` position
      # possible deletion most shifted to the left (so when last deleted nucleotide is cutpos) is
      delposleft <- (cutpos - len + 1) : cutpos # vector of all the deleted positions
      
      # we just shift the deletion to the right 1 bp at a time
      # `len` times
      
      # starting at 0 allows to get delposleft as first element of the list
      # then we shift one by one to the right
      delposL <- lapply(0:len, function(le) {
        delposleft + le
      })
      # so this is our set of possible deletions of length `len`
      # all delete the `cutpos` nucleotide
    }
    
    ### draw one deletion from the set
    drawi <- sample(1:length(delposL), 1)
    delpos <- delposL[[drawi]]
    
    # we now have our simulated deletion
    
    ### ! check the positions !
    # ! the deleted positions may be impossible
    # that is, we may be deleting negative positions (below the start of the read)
    # or positions beyond the end of the read
    # if any of these occur, repeat the procedure
    
    # if smallest deleted position is negative or 0
    # or if largest deleted position is past the end of the read
    # this is not OK and we should repeat the procedure
    if( min(delpos)<1 | max(delpos)>nchar(oref) ) {
      cat('\t \t \t \t impossible deletion, repeating procedure.\n')
      delOK <- FALSE
    } else {
      delOK <- TRUE
    }
  }
  
  # if deletion is OK, finish
  if(delOK) {
    ### create the aligned sequence
    # we split the reference sequence
    refsp <- strsplit(oref, split='')[[1]]
    # then swap the deleted nucleotides to -
    alisp <- refsp
    alisp[delpos] <- '-'
    # the reference sequence stays the same
    
    # check they are the same size just to be safe
    if(length(refsp)!=length(alisp))
      stop('\t \t \t \t >>> Error simulateDel_one: reference sequence and simulated aligned sequence do not have the same length.\n')
    
    # now we can give the reference & aligned sequences directly to recordDel from allelesToMutations.R
    # this will give the row of mut
    # this is what we return
    # overarching function can run simulateDel_one many times to create data from the simulated sample
    return( recordDel(ref=refsp,
                      ali=alisp) )
  }
  
}


# freqBins ----------------------------------------------------------------

# small utilities function
# to turn some data (e.g. deletion lengths) into frequency within bins
# v: a vector of numerics
# every: new bin every x, i.e. controls the breaks

freqBins <- function(v,
                     every,
                     lower=0,
                     last=NA) {
  
  if(is.na(last)) {
    # %% is modulo: e.g. 57 %% 5 gives 2 because 55 can be divided by 5, then remains 2
    # so we can tell how much we need to add to make it divisible based on the remainder
    remainder <- max(v, na.rm=TRUE) %% every
    
    if(remainder != 0) {
      toadd <- every - (max(v, na.rm=TRUE) %% every) + lower
    } else {
      toadd <- 0
    }
    # say every is 5
    # if max(v) was
    # 51, gives 4 > correct
    # 52, gives 3 > correct
    # 53, gives 2 > correct
    # 54, gives 1 > correct
    # 55, gives 5 > incorrect, we should add nothing
    # just catch this case
    lastbreak <- max(v, na.rm=TRUE) + toadd
  } else {
    lastbreak <- last
  }
  
  breaks <- seq(from=lower, to=lastbreak, by=every)
  
  his <- hist(v, breaks=breaks, plot=FALSE)
  hisd <- data.frame(upbound=his$breaks[2:length(his$breaks)],
                     counts=his$counts)
  # add frequencies
  hisd$freq <- hisd$counts / sum(hisd$counts)
  # dataframe
  # breaks: last value of each bin
  # counts: number of observations in this bin
  
  return(hisd)
  
}




# sumFreqBins -------------------------------------------------------------

# similar as freqBins
# but we sum frequencies within each bin
# assumes that frequencies were calculated already (unlike freqBins)
# small utilities function
# to turn some data (e.g. deletion lengths) into frequency within bins
# df: a dataframe with the values to split into bins & the frequencies
# colVals / colFreqs: name of the columns with values to split / with frequencies
# every: new bin every x, i.e. controls the breaks

sumFreqBins <- function(df,
                        colVals,
                        colFreqs,
                        every,
                        lower=0,
                        last=NA) {
  
  # values to split into bins
  v <- as.numeric(unlist(df[colVals]))
  
  # frequencies
  fs <- as.numeric(unlist(df[colFreqs]))
  
  # decide what the last break is
  if(is.na(last)) {
    # %% is modulo: e.g. 57 %% 5 gives 2 because 55 can be divided by 5, then remains 2
    # so we can tell how much we need to add to make it divisible based on the remainder
    remainder <- max(v, na.rm=TRUE) %% every
    
    if(remainder != 0) {
      toadd <- every - (max(v, na.rm=TRUE) %% every) + lower
    } else {
      toadd <- 0
    }
    # say every is 5
    # if max(v) was
    # 51, gives 4 > correct
    # 52, gives 3 > correct
    # 53, gives 2 > correct
    # 54, gives 1 > correct
    # 55, gives 5 > incorrect, we should add nothing
    # just catch this case
    lastbreak <- max(v, na.rm=TRUE) + toadd
  } else {
    lastbreak <- last
  }
  
  breaks <- seq(from=lower, to=lastbreak, by=every)
  
  bins <- cut(v, breaks=breaks, include.lowest=TRUE)
  
  sumfs <- tapply(fs, bins, FUN=sum, default=0)
  
  # prepare small dataframe to plot
  hisd <- data.frame(bin=levels(bins),
                     lowbound=breaks[1 : (length(breaks)-1)],
                     upbound=breaks[2:length(breaks)],
                     sumfreq=sumfs)
  row.names(hisd) <- NULL
  # dataframe
  # upbound: last value of each bin
  # sumfreq: sum of frequencies within this bin
  
  return(hisd)
  
}
