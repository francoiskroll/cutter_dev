#####################################################
# ~ miseqUtils: function callMutations ~
#
#
# Francois Kroll 2024
# francois@kroll.be
#####################################################

library(openxlsx)

# loops through CRISPResso output folders

# copath = folder that contains the CRISPResso results directories

# v0

# v1
# expects columns rundate, well, locus to be in meta file, but will add any other column
# note, callMutations did not change when added scaffold detection

# v2: reads cutpos, cutdist, rhapos, rhadist from meta file.
# This allows them to vary between samples, in case they are not the same amplicon or pegRNA.
# it also writes these settings in the `mut` table.

callMutations <- function(copath,
                          metapath,
                          minnreads=NA,
                          controltb=NA,
                          callSubs=TRUE,
                          exportpath) {
  
  if(!is.na(exportpath)) {
    ### check export path ends with .csv
    frmt <- substr(exportpath, start=nchar(exportpath)-2, stop=nchar(exportpath))
    if(frmt!='csv')
      stop('\t \t \t \t >>> exportpath should end with .csv.\n')
  }

  ### import meta file
  # check file exists
  if(!file.exists(metapath))
    stop('\t \t \t \t >>> Error: cannot find ', metapath, '.\n')
  # check it finishes by .xlsx
  frmt <- substr(metapath, start=nchar(metapath)-3, stop=nchar(metapath))
  if(frmt!='xlsx')
    stop('\t \t \t \t >>> Expecting meta file to be .xlsx.\n')
  # import it
  meta <- read.xlsx(metapath)
  
  # it has to have at least columns rundate & well & locus
  if(!'rundate' %in% colnames(meta))
    stop('\t \t \t \t >>> Error: expecting column "rundate" in meta file.\n')
  if(!'well' %in% colnames(meta))
    stop('\t \t \t \t >>> Error: expecting column "well" in meta file.\n')
  if(!'locus' %in% colnames(meta))
    stop('\t \t \t \t >>> Error: expecting column "locus" in meta file.\n')
  
  ### find CRISPResso output directories
  dirs <- list.dirs(copath)
  # first one is current directory, skip it
  dirs <- dirs[2:length(dirs)]
  # which ones are CRISPResso result directories?
  dirs <- dirs[startsWith(basename(dirs), 'CRISPResso')]
  
  # we have all the CRISPResso output directories
  
  # loop through rows of meta file,
  # for each, find the alleles table
  # then call mutations
  mutL <- lapply(1:nrow(meta), function(metarow) {
    
    # well is
    welli <- meta[metarow, 'well']
    # locus is
    locusi <- meta[metarow, 'locus']
    
    # from the CRISPResso folder name,
    # try to identify which is the one
    # split all the directory names
    dirsplits <- strsplit(basename(dirs), '_')
    
    # find which one has well name & locus name
    foundDir <- which(unlist(lapply(dirsplits, function(dirspli) {
      (welli %in% dirspli) & (locusi %in% dirspli)
    })))
    # there should be one & only one
    if(length(foundDir)>1) stop('\t \t \t \t Error callMutations: multiple CRISPResso directory names have well ', welli, ' and locus ', locusi,
                                '. There should only be one CRISPResso directory for each unique well/locus pair.\n')
    if(length(foundDir)==0) stop('\t \t \t \t Error callMutations: no CRISPResso directory names have well ', welli, ' and locus ', locusi, '.\n')
    
    # so the CRISPResso directory to analyse is
    cridir <- dirs[foundDir]
    
    ### find alleles table
    # path to Alleles_frequency_table.zip should be:
    alzip <- paste(cridir, 'Alleles_frequency_table.zip', sep='/')
    
    # check we found it
    if(!file.exists(alzip)) {
      cat('\t \t \t \t >>> Warning: no Alleles_frequency_table.zip in folder', cridir, '. Skipping this sample.\n')
      return()
    }
    
    # unzip it in same folder
    unzip(alzip, exdir=dirname(alzip))
    # check unzipped file exists
    # should be same as zip, but with .txt
    altxt <- paste0(substr(alzip, start=1, stop=nchar(alzip)-3), 'txt')
    if(!file.exists(altxt))
      stop('\t \t \t \t >>> Error: after unzipping, expecting Alleles_frequency_table.txt in folder', dirs[di], '.\n')
    
    # we now have everything we need
    # row in meta file with settings for filterMutations
    # and corresponding alleles table
    # we can run allelesToMutations & filterMutations
    
    ### convert alleles table to mutation table
    cat('\n \n \t \t \t \t >>> Calling mutations from', altxt,'\n')
    muttb <- allelesToMutations(alpath=altxt)
    
    ### filter the detected mutations
    # get cutpos, cutdist, rhapos, rhadist from meta file
    # (if they are given)
    # below: if not an integer, will give integer(0)
    # can check if length 0, if yes, convert to NA
    cutpos <- as.integer(meta[metarow, 'cutpos'])
    if(length(cutpos)==0) { cutpos <- NA }
    
    cutdist <- as.integer(meta[metarow, 'cutdist'])
    if(length(cutdist)==0) { cutdist <- NA }
    
    rhapos <- as.integer(meta[metarow, 'rhapos'])
    if(length(rhapos)==0) { rhapos <- NA }
    
    rhadist <- as.integer(meta[metarow, 'rhadist'])
    if(length(rhadist)==0) { rhadist <- NA }
    
    cat('\t \t \t \t >>> Filtering mutation calls.\n')
    
    cat('\t \t \t \t', 'cutpos:', cutpos, '\n')
    cat('\t \t \t \t', 'cutdist:', cutdist, '\n')
    cat('\t \t \t \t', 'rhapos:', rhapos, '\n')
    cat('\t \t \t \t', 'rhadist:', rhadist, '\n')
    
    mutf <- filterMutations(muttb=muttb,
                            minnreads=minnreads,
                            cutpos=cutpos,
                            cutdist=cutdist,
                            rhapos=rhapos,
                            rhadist=rhadist,
                            controltb=controltb,
                            callSubs=callSubs)
    
    ## add minimal meta information to mutation table
    # at least well, locus, rundate
    # add from right to left
    mutf <- mutf %>%
      mutate(well=meta[metarow,'well'], .before=1) %>%
      mutate(locus=meta[metarow, 'locus'], .before=1) %>%
      mutate(rundate=meta[metarow,'rundate'], .before=1)
    # also add rundate_well_locus to use as unique sample ID
    # I think this should always be unique, even if pooling multiple runs
    mutf <- mutf %>%
      mutate(sample=paste(rundate, well, locus, sep='_'), .before=1)
    
    # now add other meta information found in meta file
    # this allows it to be any number of other columns so free to add new ones for specific experiments
    metacols <- colnames(meta)
    # the extra meta columns are:
    metacols <- metacols[which(!metacols %in% c('rundate', 'locus', 'well'))]
    
    # add extra columns
    # cannot make mutate work within an sapply loop, despite help from ChatGPT
    # but below works
    for(colnm in metacols) {
      mutf <- mutf %>%
        mutate(tmp=meta[metarow, colnm], .before='well')
      # now put the correct column name
      colnames(mutf)[which(colnames(mutf)=='tmp')] <- colnm
    }
    
    ### return filtered mutations
    return(mutf)
  })
  
  # >>> we get list mutL which is table of filtered mutations for each sample
  # gather everything in one dataframe
  # we have locus & well to keep track of which mutation is from which sample
  mut <- do.call(rbind, mutL)
  
  # save this dataframe
  if(!is.na(exportpath)) {
    write.csv(mut, exportpath, row.names=FALSE)
    cat('\t \t \t \t >>> Wrote', exportpath, '\n')
  } else {
    cat('\t \t \t \t >>> Export is off.')
  }

  # we also return mut
  invisible(mut)
}
