#####################################################
# ~ cutter: function removeExpedit ~
#
#
# Francois Kroll 2025
# francois@kroll.be
#####################################################

source('~/Dropbox/cutter/isEditPresent.R')

# vector of mutids, remove expedit
# usually simple but can be more complicated for substitutions
# imagine
# expedit is sub_185_185_2_G_A
# mutid is sub_184_185_2_GG_AA

removeExpedit <- function(expedit,
                          mutid) {
  
  ### loop through mutids
  mutid_noexp <- sapply(mutid, function(mutone) {
    
    ## if it is exactly the expedit, simply remove it
    if(mutone==expedit) {
      return()
      
      ## if insertion or deletion, return it directly
      # (we checked it is not the expedit above)
    } else if(str_detect(mutone, 'ins') | str_detect(mutone, 'del')) {
      return(mutone)
      
      ## if substitution but expedit is not a substitution
      # just return it
    } else if(str_detect(mutone, 'sub') & !str_detect(expedit, 'sub')) {
      return(mutone)
      
      ## if substitution & expedit is also a substitution
      # need to be careful
    } else if(str_detect(mutone, 'sub') & str_detect(expedit, 'sub')) {
      
      ## is the expedit included in the substitution?
      # if not, leave mutone unchanged
      if(!isEditPresent(expedit=expedit, mutid=mutone)) {
        return(mutone)
        
        ## but if expedit is present,
        # need to proceed by position
      } else {
        
        # split the expedit into 1-nt mutids
        expedit_si <- subToSingleNt(mutid=expedit)
        
        # split the mutone into 1-nt mutids
        mutone_si <- subToSingleNt(mutid=mutone)
        
        # note, we tested above (by isEditPresent) that expedit is *entirely* present
        # so, if some positions correspond to the edit but expedit is not entirely present,
        # these will be left and called unwanted
        # (which I think is how it should be!)
        # remove the expedit
        # i.e. keep from mutone_si the positions that are not included in expedit_si
        # there are the surplus (unwanted) substitutions
        mutone_sinoexp <- mutone_si[!mutone_si %in% expedit_si]
        cat('\t \t \t \t \t >>> removeExpedit: case where mutid is expedit & other substitution(s).\n')
        cat('\t \t \t \t \t >>> expedit is', expedit, '; mutid is', mutone, '\n')
        cat('\t \t \t \t \t >>> removing expedit excises the substitution to', mutone_sinoexp, '\n')
        # we keep mutone_sinoexp as separated substitutions, I think is easier to deal with at this stage
        return( mutone_sinoexp )
      }
    } else {
      stop('\t \t \t \t >>> Error removeExpedit: unexpected exception.\n')
    }
  })
  
  # now we have mutid_noexp
  # these are mutids without the expedit
  # substitution may have been split
  mutid_noexp <- as.character(unlist(mutid_noexp)) # just making sure it is a simple vector of characters
  return(mutid_noexp)
  
}
