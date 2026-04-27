#####################################################
# ~ ZPRI: mutations to read labels ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

source('~/Dropbox/cutter/isEditPresent.R')
source('~/Dropbox/cutter/removeExpedit.R')

# expect mut, mutation table created by filterMutations

# v1

# v2: scaffold detection

# scaffWin: distance before/after end of RHA to detect scaffold
# for now, detect scaffold just means insertion or substitution that starts with G
# ! definitions of classes are slightly different:
# pure = edit / no mutation / no scaffold
# impure = edit / mutation / no scaffold
# mutated = no edit / mutation / no scaffold
# scaffold = edit or not / mutation / scaffold
# so scaffold class is agnostic to edit or not, i.e. whenever scaffold is detected that is the class the read gets

# v3: correction of scaffold detection, taking into account PE strand

# v4: scaffDetect argument, to turn ON or OFF scaffold detection

# v5: changing name to preciseClassify_one
# will allow multiple modes in classifyReads, mode='precise' to detect classify reads when expected edit (e.g. prime editing)

# v6: added argument unwantedSubs
# changed slightly how detected unwanted mutations as a result

# v7: changes to scaffold detection

preciseClassify_one <- function(mut,
                                scaffDetect,
                                pestrand,
                                whichScaff,
                                scaffWin,
                                unwantedSubs) {
  
  ### check mut is giving data for just one sample
  if(length(unique(mut$sample))>1)
    stop('\t \t \t \t >>> Error: attempting to give mutation data for more than one sample.\n')
  
  # record rundate, sid, locus, well, grp sample, sample's coverage
  rundate <- unique(mut$rundate)
  sid <- unique(mut$sid)
  locnm <- unique(mut$locus)
  wellnm <- unique(mut$well)
  grpnm <- unique(mut$grp)
  splnm <- unique(mut$sample)
  
  ### get the original reference sequence
  # take all reference sequences that are not NA
  # remove any -
  orefs <- gsub('-', '', mut[!is.na(mut$ref), 'ref'])
  # check they are all the same
  if(length(unique(orefs))>1)
    stop('\t \t \t \t >>> Error preciseClassify_one: more than one original reference sequence in this sample\'s mutation table, which does not make sense.\n')
  oref <- unique(orefs)
  
  # also record rhapos & pestrand (if scaffDetect is ON)
  if(scaffDetect) {
    rhapos <- unique(mut$rhapos)
    pestrand <- unique(mut$pestrand)
  }
  
  # for expedit, record it if present or turn to NA if not
  if('expedit' %in% colnames(mut)) {
    expedit <- unique(mut$expedit)
  } else {
    cat('\t \t \t \t >>> classifyReads: column "expedit" not in mut, will assume that there is no expected edit to detect.\n')
    expedit <- NA
  }
  
  # check only one sample's coverage
  splcov <- unique(mut$splcov)
  if(length(splcov)!=1)
    stop('\t \t \t \t >>> Error: 0 or multiple sample coverage.\n')
  
  ###### preparing scaffold detection ######
  
  # note, some of below could be in classifyReads
  # but I prefer keeping the nitty-gritty of scaffold detection here
  # should not make much of a difference in efficiency to run the checks for each sample
  
  ### scaffDetect: we checked in classifyReads that columns rhapos & pestrand look correct
  # we also turned expedit to NA if column "expedit" was not present
  ### if scaffDetect is ON, check whichScaff makes sense
  if(scaffDetect) {
    if(! whichScaff %in% c('std', 'std2', 'altAAAA', 'altUAAA', 'altGAAA'))
      stop('\t \t \t \t This scaffold is not recorded. Options are: "std", "std2", "altAAAA", "altUAAA", "altGAAA".\n')
    
    # then record scaffSeq
    # ! this is not actually sequence of the scaffold
    # this is what the insertion of the entire scaffold would look like
    # if RT had read all of the scaffold
    # ! assuming PE on Forward strand
    if(whichScaff=='std') {
      scaffSeq <- 'GGACCGACTCGGTCCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC'
    } else if(whichScaff=='std2') {
      scaffSeq <- 'GCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC'
    } else if(whichScaff=='altAAAA') {
      scaffSeq <- 'GCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATGCTGTTTCCAGCATAGCTCTAAAAC'
    } else if(whichScaff=='altUAAA') {
      scaffSeq <- 'GCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTAAACTTGCTATGCTGTTTCCAGCATAGCTCTTAAAC'
    } else if(whichScaff=='altGAAA') {
      scaffSeq <- 'GCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTCAACTTGCTATGCTGTTTCCAGCATAGCTCTGAAAC'
    }
    # note, alt scaffolds could also start GG instead of GC
    # we could add alt2AAAA, alt2UAAA, alt2GAAA
    # in practice, unlikely that insertion reach the part that is different in alt scaffold
    # but user should be a little more careful about start being GG or GC
    
    # if PE on Forward strand, sequence above is the one we expect to find
    # if PE on Reverse strand, we actually expect the reverse-complement of sequence above
    if(pestrand=='reverse') {
      scaffSeq <- reversecomplement(scaffSeq)
    }
    # I did reversecomplement to keep it consistent that we read the mutation from left to right on the Forward strand
    # but in practice, scaffold insertion in this case should be read from right to left starting from last nucleotide
  # else = if scaffDetect is FALSE
  } else {
    whichScaff <- NA
    scaffSeq <- NA
    scaffWin <- NA
  }
  
  ### if scaffDetect, set the scaffWin default
  # bit confusing for user to set themselves
  # if PE on Forward strand: theoretical scaffold insertion starts at rhapos
  # if PE on Reverse strand: theoretical scaffold insertion stops at rhapos - 1
  # then it also differs for substitutions, see README figure
  # is.na(scaffWin): if user did not set it themselves
  if( scaffDetect & is.na(scaffWin[1]) & pestrand=='forward' ) {
    scaffWin <- c(0, 0)
  } else if( scaffDetect & is.na(scaffWin[1]) & pestrand=='reverse') {
    scaffWin <- c(-1, -1)
  }
  
  ### scaffWin adjustment ###
  if(scaffDetect) {
    
    # now gets a bit tricky...
    # if nucleotides in reference sequence happen to match the expected scaffold incorporation
    # then we should adjust the detection window
    
    # next steps change slightly whether prime editing is on the Forward or Reverse strand
    # overall logic: we want to search how many nucleotides from the scaffold match the sequence after (if forward) or before (if reverse) rhapos
    # we try one bp at a time
    # nti is nth nucleotide after/before the rhapos
    # will end up being = number of nucleotides from reference sequence that happen to match the expected scaffold incorporation
    # initiate at 0, but we will look first at nucleotide +1 (Forward) or -1 (Reverse)
    # meaning the nucleotide just after/before rhapos (which is last nucleotide read from RTT, assuming RT stops where it should)
    
    # for insertions: it makes sense to extend the search window here because we cannot predict what the alignment algorith will do
    # so no other choice than extending the search window
    # i.e. we cannot be as stringent
    
    # scaffWinAdj is scaffWin adjusted
    # scaffWin stays original: c(0,0) or c(-1,-1) or set by the user
    scaffWinAdj <- scaffWin
    
    ### if FORWARD ###
    if(pestrand=='forward') {
      
      nti <- 0
      # nextmatches will say TRUE if next nucleotide matches / FALSE if next nucleotide does not match
      nextmatches <- TRUE
      while(nextmatches) {
        # we check next nucleotide after rhapos (rhapos + nti + 1)
        # and compare it to next nucleotide in expected scaffold incorporation
        # do they match?
        if( substr(oref, rhapos+nti+1, rhapos+nti+1) == substr(scaffSeq, nti+1, nti+1) ) {
          # if yes, try again with next nucleotide
          nextmatches <- TRUE
          nti <- nti + 1
        } else {
          # if not, we stop and return current nti
          nextmatches <- FALSE
        }
      }
      # works well for forward, checked
      # >> nti is number of nucleotides from reference, after rhapos, that happen to match a scaffold insertion
      # we want to adjust scaffWin based on this; detection window should be bigger on the right if sequence happens to match
      # we simply add this many nucleotides to the end of the window
      if(nti>0) {
        cat('\t \t \t \t \t', nti, 'nucleotide(s) before rhapos happen to match the expected scaffold incorporation, adjusting the search window...\n')
      }
      scaffWinAdj[2] <- scaffWinAdj[2] + nti
      # for example, 2 nt match, so scaffWin goes from c(0,0) to c(0,2)
      
      ### if REVERSE ###
    } else if(pestrand=='reverse') {
      
      # logic is the same, but we just slide to the left
      # above, we changed scaffSeq to its reversecomplement if we were editing the Reverse strand
      # ! we want to read this sequence (longest possible scaffold insertion) from left to right
      nti <- 0
      # nextmatches will say TRUE if next nucleotide matches / FALSE if next nucleotide does not match
      # here, next means next nucleotide to the left
      nextmatches <- TRUE
      while(nextmatches) {
        # first we check if last nucleotide of possible scaffold insertion matches the next nucleotide left of rhapos
        # then we check if before-last nucleotide of possible scaffold insertion matches next nucleotide on the left
        # etc.
        if( substr(oref, rhapos-nti-1, rhapos-nti-1) == substr(scaffSeq, nchar(scaffSeq)-nti, nchar(scaffSeq)-nti) ) {
          nextmatches <- TRUE
          nti <- nti + 1
        } else {
          # if not, we stop and return current nti
          nextmatches <- FALSE
        }
      }
      # works well for reverse, checked
      # >> nti is number of nucleotides from reference, before rhapos, that happen to match a scaffold insertion
      # we want to adjust scaffWin based on this; detection window should be bigger on the left if sequence happens to match
      # we simply substract this many nucleotides to the end of the window
      
      # works well for reverse, checked
      # so nti is number of nucleotides from reference, after rhapos (to the left), that happen to match a scaffold insertion
      # we want to adjust scaffWin based on this; detection window should be bigger on the left if sequence happens to match
      # we simply subtract this many nucleotides to the start of the window
      if(nti>0) {
        cat('\t \t \t \t', nti, 'nucleotide(s) before rhapos happen to match the expected scaffold incorporation, adjusting the search window... \n')
      }
      scaffWinAdj[1] <- scaffWinAdj[1] - nti
      # for example, 2 nt match, so scaffWin goes from c(-1,-1) to c(-3,-1)
      
    }
    
    
    ### interim sanity check/summary
    # Forward strand: scaffold insertion should in theory *start* exactly at rhapos,
    # because start position of insertion is last aligned nucleotide (just before first inserted nucleotide)
    #
    # Reverse strand: scaffold insertion should in theory *end* (and start) exactly at rhapos - 1
    # first extra nucleotide will be nucleotide to the left of rhapos
    
    # >> we already accounted for this minus 1 difference on reverse strand above
    # when setting the defaults for scaffWin
    
    # positions where we want to look are:
    scaffdetectposAdj <- (rhapos + scaffWinAdj[1]) : (rhapos + scaffWinAdj[2])
    # Adj for adjusted
    
    # non-adjusted is:
    scaffdetectpos <- (rhapos + scaffWin[1]) : (rhapos + scaffWin[2])
    
    # e.g. Forward, rhapos is 103, scaffWin is c(0,2)
    # scaffdetectpos becomes 103 + 0 >> 103 + 2
    # i.e. 103, 104, 105
    
    # e.g. Reverse, rhapos is 103, scaffWin is c(-3,-1)
    # scaffdetectpos becomes 103 - 3 >> 103 - 1
    # i.e. 100, 101, 102
    
  } # closes if(scaffDetect)
  # >> we are done with scaffWin adjustment
  # we have scaffdetectpos & scaffdetectposAdj
  
  
  ######################### below is lapply read by read within that sample
  
  # gather unique read IDs
  uids <- unique(mut$rid)
  # loop through unique read IDs
  rlabsL <- lapply(uids, function(ui) {
    
    # for debug
    # cat('ui is', ui, '\n')
    
    ## preallocate vector which is
    # 1) prime-edit present?; 2) other mutation present?; 3) scaffold present?
    # preallocate as FALSE, FALSE, FALSE so = reference read
    rlab <- c(expedit=FALSE, mut=FALSE, scaffold=FALSE)
    
    ## gather all the detected mutations in this read
    muti <- mut %>%
      filter(rid==ui)
    
    ### if mutation "ref" is present ###
    if('ref' %in% unique(muti$mutid)) {
      # there should not be anything else, check that
      if(nrow(muti)!=1)
        stop('\t \t \t \t >>> Error: for read', ui,', there is mutation "ref" and other mutations, which does not make sense.\n')
      # and we can return vector rlab as it is, i.e. FALSE, FALSE, FALSE
      return(rlab)
    } # we are done for this read as we called return, below is `else`
    # so everything below: read has at least one mutation
    # and we know "ref" is not there
    
    # note, read could still be called reference below
    # for example, a read that just has one substitution and no edit will be called reference (assuming unwantedSubs is FALSE)
    
    # collect all the mutids currently in muti
    mutids <- muti$mutid
    
    ### is desired edit present? ###
    # if expedit is NA, just keep expedit as FALSE
    if(isEditPresent(expedit=expedit, mutid=mutids) ) {
      # then switch edit flag to TRUE
      rlab['expedit'] <- TRUE
    }
    
    ### other mutation(s) present? ###
    # should we count substitutions as unwanted mutation or not?
    # if not, then exclude substitutions from mutids
    # copy mutids
    mutids_unw <- mutids
    
    if(!unwantedSubs) {
      mutids_unw <- mutids[which(!str_detect(mutids, 'sub'))]
      # unw for unwanted mutations
    }
    
    ## also need to remove expedit so we do not call it unwanted
    mutids_unw <- removeExpedit(expedit=expedit, mutid=mutids_unw)
    
    ## now count number of mutids left
    # if any, then there is one or more unwanted mutations
    if(length(mutids_unw)>0) {
      rlab['mut'] <- TRUE
    }
    # note, if read only has substitution(s) (no expected edit)
    # and unwantedSubs is FALSE
    # then all three flags are still FALSE and read will be called reference below
    
    #########
    
    ### detect scaffold incorporation ###
    if(scaffDetect) {
      
      # see notes in README, there are many difficult exceptions
      # broad logic is
      
      # PE on Forward:
        # insertion that starts at rhapos and starts with G
        # substitution that starts at rhapos + 1 and starts with G
      
      # PE on Reverse:
        # insertions that stops at rhapos - 1 and ends with C
        # subsitution that stops at rhapos - 1 and ends with C
      if(pestrand=='forward') {
        
        # annoyingly, it is different whether scaffold incorporation is an insertion or a substitution
        # if insertion: we expect scaffold at rhapos
        # if substitution: we expect scaffold at rhapos + 1
        # test each separately, and shift by +1 when we look for substitutions
        
        ### INSERTION
        scaff_ins <- muti %>%
          filter(type == 'ins') %>%
          filter(start %in% scaffdetectposAdj) %>% # we use adjusted window
          # typically just one position: rhapos
          # but may have been extended slightly
          filter(startsWith(altseq,
                            substr(scaffSeq, 1, 1)))
        
        # >>> did we find something? if yes & mutation is longer than 1 bp, we can also check that the second nucleotide matches
        # this will further reduce false positive detections
        if(nrow(scaff_ins)==1) {
          if(scaff_ins$bp>1) {
            scaff_ins <- scaff_ins %>%
              filter( substr(altseq, 2, 2) == substr(scaffSeq, 2, 2) )
          }
        }
        
        ### SUBSTITUTION
        # ! if substitution, it does not make sense to use the adjusted window
        # we know what the alignment algorithm will do so we should stay as stringent (not extend the window)
        # but we do want to *shift* the window if some nucleotides match
        scaff_sub <- muti %>%
          filter(type == 'sub') %>%
          filter( start %in% (scaffdetectpos + 1 + nti) ) %>%
          # we use original window, but shifted if some nucleotides match
          # typically c(0,0); then e.g. nti = 2 would be c(2,2)
          filter(startsWith(altseq,
                            substr(scaffSeq, 1+nti, 1+nti)))
          # e.g. if 2 nt match, the first mismatch will not be the first nucleotide of the scaffold (G)
          # rather, we would expect to find the 3rd nucleotide of the scaffold
          # if nti = 0, we just look for first nucleotide, i.e. G
        
        # >>> did we find something? if yes & mutation is longer than 1 bp, we can also check that the second nucleotide matches
        # this will further reduce false positive detections
        if(nrow(scaff_sub)==1) {
          if(scaff_sub$bp>1) {
            scaff_sub <- scaff_sub %>%
              filter( substr(altseq, 2, 2) == substr(scaffSeq, 2+nti, 2+nti) )
          }
        }
        
        scaff <- rbind(scaff_ins, scaff_sub)
        
      } else if(pestrand=='reverse') {
        # for reverse, it does not make a difference to stop position whether insertion or substitution
        # in case of substitution, last mismatch would also be rhapos - 1
        # but, we still need to separate the two for a few reasons, see below
        
        ### INSERTION
        scaff_ins <- muti %>%
          filter(type == 'ins') %>%
          filter(stop %in% scaffdetectposAdj) %>% # *** here stop, not start; using adjusted window
          filter(endsWith(altseq,
                          substr(scaffSeq, nchar(scaffSeq), nchar(scaffSeq)))) # *** here endsWith (not startsWith) last nucleotide of scaffSeq
        
        # >>> did we find something? if yes & mutation is longer than 1 bp, we can also check that the second nucleotide matches
        # this will further reduce false positive detections
        if(nrow(scaff_ins)==1) {
          if(scaff_ins$bp>1) {
            scaff_ins <- scaff_ins %>%
              filter( substr(altseq, nchar(altseq)-1, nchar(altseq)-1) == substr(scaffSeq, nchar(scaffSeq)-1, nchar(scaffSeq)-1) )
          }
        }
        
        ### SUBSTITUTION
        scaff_sub <- muti %>%
          filter(type == 'sub') %>%
          filter( stop %in% (scaffdetectpos - nti) ) %>% # *** here stop, not start; using original window
          # but window may be shifted if some nucleotides match
          # typically c(-1,-1); then e.g. nti = 2 would be c(-3,-3)
          filter(endsWith(altseq,
                          substr(scaffSeq, nchar(scaffSeq)-nti, nchar(scaffSeq)-nti))) # *** here endsWith, not startsWith
        # e.g. scaffSeq is ...GTCGGTCC and nti is 2 (so CC match reference sequence)
        # then we expect to find T as first mismatch
        # which is correct: nchar would be 8, minus 2 is 6, gives T
        
        # >>> did we find something? if yes & mutation is longer than 1 bp, we can also check that the second nucleotide matches
        # this will further reduce false positive detections
        if(nrow(scaff_sub)==1) {
          if(scaff_sub$bp>1) {
            scaff_sub <- scaff_sub %>%
              filter( substr(altseq, nchar(altseq)-1, nchar(altseq)-1) == substr(scaffSeq, nchar(scaffSeq)-nti-1, nchar(scaffSeq)-nti-1) )
          }
        }
        
        scaff <- rbind(scaff_ins, scaff_sub)
        # ! if e.g. 2 nt match, we expect to find 3rd nucleotide of the scaffold, counting from the last
        # just "nti" not minus 1, as in "nti jumps" from last
        # if nti = 0, we just look for first nucleotide, i.e. C
      }
      # whether forward or reverse, we know have 'scaff' which stores mutations which correspond to scaffold incorporation
      # if 'scaff' has any rows, this means the read has a scaffold incorporation
      
      if(nrow(scaff)>0) {
        
        cat('\t \t \t \t \t >>> scaffold incorporation:', scaff$mutid, '\n')
        
        rlab['scaffold'] <- TRUE
        # note, if unwantedSubs is FALSE (i.e. only call indels)
        # scaffold can be present and mut flag be FALSE
      }
    } # closes scaffold detection
    
    ## return labels for this read
    return(rlab)
  }) # closes lapply to check read by read
  
  # create a simpler dataframe
  rlabs <- data.frame(do.call(rbind, rlabsL))
  
  # add back the read ids
  rlabs <- rlabs %>%
    mutate(rid=uids, .before=1)
  
  # we now have the read-by-read labels
  # assign categories based on those labels
  rcats <- sapply(1:nrow(rlabs), function(r) {
    # from row, only keep the flags, i.e. expedit / mut / scaffold
    # just a vector of three booleans
    rl <- as.logical( rlabs[r, c('expedit', 'mut', 'scaffold')] )
    # add back the names
    names(rl) <- c('expedit', 'mut', 'scaffold')
    
    ### no edit, no unwanted mutations, no scaffold = REFERENCE
    if(!rl['expedit'] & !rl['mut'] & !rl['scaffold']) {
      return('reference')
      
    #### *edit*, no unwanted mutations, no scaffold = PURE EDIT
    } else if(rl['expedit'] & !rl['mut'] & !rl['scaffold']) {
      return('pure')
      
    ### *edit*, but also *unwanted mutations* (that is not scaffold) = IMPURE EDIT
    } else if(rl['expedit'] & rl['mut'] & !rl['scaffold']) {
      return('impure')
      
    ### no edit and *unwanted mutations* (that is not scaffold) = MUTATED
    # if unwantedSubs is FALSE, this category can also be called INDEL
    } else if(!rl['expedit'] & rl['mut'] & !rl['scaffold']) { # if no edit, but mutation (no scaffold)
      return('mutated')
      
    ### category SCAFFOLD trumps everything
    # if we see scaffold, we call that read scaffold, whatever the situation is regarding edit/other unwanted mutations
    } else if(rl['scaffold']) { # if scaffold (agnostic re edit or unwanted mutations)
      return('scaffold')
    }
  })
  
  # add the read categories to rlabs
  rlabs <- rlabs %>%
    mutate(cat=rcats)
  
  ### add meta columns to rlabs
  # add from right to left in final columns
  rlabs <- rlabs %>%
    mutate(grp=grpnm, .before=1) %>%
    mutate(well=wellnm, .before=1) %>%
    mutate(locus=locnm, .before=1) %>%
    mutate(sid=sid, .before=1) %>%
    mutate(rundate=rundate, .before=1) %>%
    mutate(sample=splnm, .before=1) %>%
    mutate(splcov=splcov, .after='grp')
  
  return(rlabs)
  
}