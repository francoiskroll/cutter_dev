#####################################################
# ~ cutter: function isEditPresent ~
#
#
# Francois Kroll 2025
# francois@kroll.be
#####################################################

source('~/Dropbox/cutter/utilities.R')


# subToSingleNt -----------------------------------------------------------
# small helper function
# split a substitution mutid (or expedit) as multiple single-nucleotide mutid
subToSingleNt <- function(mutid) {
  
  # we can start by writing mutid as a df
  startpos <- strNthSplit(mutid, '_', 2)
  stoppos <- strNthSplit(mutid, '_', 3)
  refseq_s <- strsplit(strNthSplit(mutid, '_', 5), '')[[1]]
  altseq_s <- strsplit(strNthSplit(mutid, '_', 6), '')[[1]]
  mutid_df <- as.data.frame(rbind(refseq_s, altseq_s))
  colnames(mutid_df) <- startpos:stoppos
  rownames(mutid_df) <- NULL
  
  # then loop through columns
  singleNts <- sapply(1:ncol(mutid_df), function(mucol) {
    paste('sub',
          colnames(mutid_df)[mucol], # start
          colnames(mutid_df)[mucol], # stop (same as start as single-nt)
          1, # length
          mutid_df[1, mucol], # refseq
          mutid_df[2, mucol], # altseq
          sep='_')
  })
  
  return(singleNts)
}


# isEditPresent -----------------------------------------------------------

# should have expedit / mutations (mutid) in one read
# both can be multiple

isEditPresent <- function(expedit,
                          mutid) {
  
  ### exception: expedit may be NA
  if(is.na(expedit)) {
    return(FALSE)
  } # below: expedit is not NA
  
  ### when expedit is an insertion or a deletion:
  # simple case, just look if expedit is in the mutids
  if(str_detect(expedit, 'ins') | str_detect(expedit, 'del')) {
    return( expedit %in% mutid )
    
  ### when expedit is a substitution,
  # we need to be more careful
  # e.g. expedit is sub_185_185_1_G_A and mutid is sub_184_185_2_GG_AA
  # we should call edit (possibly call it impure edit) as second G>A is the expected edit
  # same if expedit was sub_184_184_1_G_A
  } else {
    
    ## split expedit as multiple 1-nt mutids
    expedit_si <- subToSingleNt(mutid=expedit)
    
    ## now loop through mutids
    # for each, check if it may be the expedit
    isexpedit <- sapply(mutid, function(mutone) {
      
      ### if insertion or deletion, just check
      # (will be FALSE as we know expedit is a substitution)
      if(str_detect(mutone, 'ins') | str_detect(mutone, 'del')) {
        return( expedit==mutone )
        
        ### substitution, need to check by single positions
      } else if (str_detect(mutone, 'sub')) {
        
        ## split mutone into 1-nt mutids, as for expedit
        mutone_si <- subToSingleNt(mutid=mutone)
        
        ## now we simply need to check that all 1-nt mutids of expedit
        # are in the 1-nt mutids of mutone_si
        return( all(expedit_si %in% mutone_si) )
      }
    })
    # vector isexpedit is, for each mutid, TRUE (it is expedit) or FALSE (is not expedit)
    # if any TRUE, return TRUE
    if(TRUE %in% isexpedit) {
      return(TRUE)
    } else {
      return(FALSE)
    }

  }
}
