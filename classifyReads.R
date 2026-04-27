#####################################################
# ~ miseqUtils: function classifyReads ~
#
#
# Francois Kroll 2024
# francois@kroll.be
#####################################################

source('~/Dropbox/cutter/reversecomplement.R') # preciseClassify_one needs it

# v1

# v2: with classifyReads_one_v2.R which includes scaffold detection

# v3: added 'pestrand' argument ('forward' or 'reverse') to adjust scaffold detection.
# in v2, would expect starting with G, which is only correct if PE is happening on Forward strand
# see more notes in README.md

# v4: can turn ON or OFF scaffold detection with argument scaffDetect

# v5: takes rhapos from meta file

# v6: added argument unwantedSubs

classifyReads <- function(mut,
                          mode='precise',
                          scaffDetect=FALSE,
                          whichScaff='std',
                          scaffWin=NA,
                          unwantedSubs=FALSE,
                          nreadsSample=NA,
                          exportpath) {
  
  ### check export path ends with .csv
  if(!is.na(exportpath)) {
    frmt <- substr(exportpath, start=nchar(exportpath)-2, stop=nchar(exportpath))
    if(frmt!='csv')
      stop('\t \t \t \t >>> exportpath should end with .csv.\n')
  }
  
  ### check mode given
  if(! mode %in% c('precise', 'precise_simple', 'frameshift', 'insdel'))
    stop('\t \t \t \t >>> Error classifyReads: mode should be one of the following: "precise", "simple_precise", "frameshift".\n')
  
  if(mode=='precise') {
    ### if scaffold detection is ON,
    # check rhapos is given in the mut table so we know where to look
    # check pestrand is giving in the mut table
    if(scaffDetect) {
      ## check rhapos
      if(! 'rhapos' %in% colnames(mut) ) stop('\t \t \t \t >>> Error classifyReads: scaffold detection is ON (scaffDetect=TRUE) but column "rhapos" is not in mut table so do not know where to look.\n')
      
      ## check pestrand
      if(! 'pestrand' %in% colnames(mut) ) stop('\t \t \t \t >>> Error classifyReads: scaffold detection is ON (scaffDetect=TRUE) but column "pestrand" is not in mut table. Make sure it is given in the meta.xlsx file when running callMutations.\n')
      if(! all(mut$pestrand %in% c('forward', 'reverse'))) stop('\t \t \t \t >>> Error classifyReads: in column "pestrand", only acceptable values are "forward" and "reverse".\n')
    }
  }
  
  ### import mut
  # if is a character, assume we are given a path,
  # if not, assume we are given R object directly
  if(is.character(mut)) {
    # check it is .csv
    frmt <- substr(mut, start=nchar(mut)-2, stop=nchar(mut))
    if(frmt!='csv')
      stop('\t \t \t \t >>> If giving a path for mut, it should end with .csv.\n')
    mut <- read.csv(mut)
  }
  
  ### split into mutations by sample
  # to get a list where each slot is one sample
  mutL <- split(mut, mut$sample)
  # then we pass each mutation table to classifyReads_one
  rlabL <- lapply(1:length(mutL), function(muti) {
    cat('\n \t \t \t \t >>> Sample', muti, 'out of', length(mutL), '\n')
    cat('\t \t \t \t', names(mutL)[muti], '\n')
    mut <- mutL[[muti]]
    
    ### sample fewer reads to analyse?
    if(!is.na(nreadsSample)) {
      
      # check integer or double
      if(!is.integer(nreadsSample) & !is.double(nreadsSample))
        stop('\t \t \t \t Error classifyReads: nreadsSample should be NA or an integer, e.g. 1000.\n')
      
      # ! if asking to sample more reads than there are in this sample, just keep all reads
      # how many reads do we have?
      nrids <- length(unique(mut$rid))
      
      # if trying to sample more, tell user and keep every read
      if(nreadsSample > nrids) {
        cat('\t \t \t \t >>> Trying to sample', nreadsSample, ', but there are only', nrids, 'reads in this sample. We will analyse all.\n')
        splcov_new <- unique(mut$splcov)
      } else {
        # sample nreads from the unique reads in mut
        rid2keep <- sample(unique(mut$rid), nreadsSample, replace=FALSE)
        
        mut <- mut %>%
          filter(rid %in% rid2keep)
        splcov_new <- nreadsSample
      }
    }
    
    ### mode PRECISE
    if(mode=='precise') {
      
      ## about 'pestrand' argument
      pestrand <- unique(mut$pestrand)
      
      # we should not have more than one 'pestrand'
      if(length(pestrand)>1)
        stop('\t \t \t \t >>> Error classifyReads: more than one "pestrand" value within one sample, which does not make sense.\n')
      
      # ! 'pestrand' may not be not given, make sure it is NA in that case
      if(length(pestrand)==0) {
        pestrand <- NA
      }
      
      rlab <- preciseClassify_one(mut,
                                  scaffDetect=scaffDetect,
                                  pestrand=pestrand,
                                  whichScaff=whichScaff,
                                  scaffWin=scaffWin,
                                  unwantedSubs=unwantedSubs)
      
    ### mode PRECISE_SIMPLE
    } else if(mode=='precise_simple') {
      rlab <- simpleClassify_one(mut,
                                 unwantedSubs=unwantedSubs)
    
    ### mode FRAMESHIFT
    } else if(mode=='frameshift') {
      rlab <- frameshiftClassify_one(mut=mut)
    
      
    ### mode INSDEL
    } else if(mode=='insdel') {
      rlab <- insdelClassify_one(mut=mut)
    }
    
    ### if we used nreadsSample, we need to update splcov
    if(exists('splcov_new')) {
      rlab$splcov <- splcov_new
    }
    
    ### return into rlabL
    return(rlab)

  })
  
  ### gather in one dataframe
  rlab <- do.call(rbind, rlabL)

  # make sure it does not add row names
  row.names(rlab) <- NULL
  
  # add mode used as a column so we can use that info later
  rlab <- rlab %>%
    mutate(mode=mode, .after='rid')
  
  # export
  if(!is.na(exportpath)) {
    write.csv(rlab, exportpath, row.names=FALSE)
    cat('\t \t \t \t >>> Wrote', exportpath, '\n')
  } else {
    cat('\t \t \t \t >>> Export is OFF.')
  }

  # return
  invisible(rlab)
  
}
