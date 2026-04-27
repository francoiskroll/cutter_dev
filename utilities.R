#####################################################
# ~ cutter: common utilities ~
#
#
# Francois Kroll 2025
# francois@kroll.be
#####################################################

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