#####################################################
# ~ ZPRI: function to plot as deletionsbyStrand.R ~
#
#
# Francois Kroll 2025
# francois@kroll.be
#####################################################

# v0, first version of function

# v1
# option to export png, pdf is heavy and difficult to use


# other functions ---------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)

source('~/Dropbox/cutter/detectMHdel.R')

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


# -------------------------------------------------------------------------



ggIGVdel <- function(mut,
                     rangeStart=NA,
                     rangeStop=NA,
                     bpSides=NA,
                     titleOrNo=TRUE,
                     orderdels='start',
                     ligPos=FALSE,
                     forceligPos=NA,
                     avgLigPos=FALSE,
                     addBottomMargin=0,
                     exportPath,
                     width,
                     height) {
  
  ### delete any edit
  # i.e. keep only when mutid is different than expedit
  # we directly keep only deletions
  mut <- mut %>%
    filter(is.na(expedit) | expedit!=mutid) %>% # ! can sometimes be NA, e.g. standard Cas9
    filter(type=='del')
  
  ### keep only reads with a single deletion
  # create new column urid for sample_rid
  mut <- mut %>%
    mutate(urid=paste(sample, rid, sep='_'), .after='rid')
  # count number of deletions by urid
  deltal <- mut %>%
    group_by(urid) %>%
    tally(name='ndel')
  # delete any where more than one deletions, keeps things simple
  urid2del <- unique(as.character(unlist( deltal[deltal$ndel>1, 'urid'] )))
  
  # remove those reads
  mut <- mut %>%
    filter(! urid %in% urid2del)
  
  # v0 was running detectMHdel here but actually not using the output!
  
  ### keep only unique deletion alleles
  # ! within groups
  # e.g. if we see the same deletion in D10A and H840A, we plot once for each
  del <- mut %>%
    filter(type=='del') %>%
    mutate(grpmut=paste(grp, mutid, sep='.'), .before='mutid') %>% # e.g. Cas9.del_68_72_5_CGTAC_NA
    distinct(grpmut, .keep_all=TRUE)
  
  ### prepare data to plot
  # each row is one allele
  # each column is one position on reference
  
  # preallocate it as all TRUE, i.e. nucleotide is present
  # nrow(del) gives number of unique alleles
  
  # original reference sequence is
  oref <- unique(gsub('-', '', del$ref))
  
  # alls for alleles
  alls <- as.data.frame(matrix(nrow=nrow(del),
                               ncol=nchar(oref)))
  
  colnames(alls) <- 1:ncol(alls)
  
  # fill it all with TRUE
  # i.e. nucleotide is present
  alls[,] <- TRUE
  
  # to fill it, we check each allele in del
  # find deleted positions, switch those to FALSE in alls matrix, i.e. not present (deleted)
  for(al_i in 1:nrow(del)) {
    
    # deleted positions are:
    delpos <- del[al_i, 'start'] : del[al_i, 'stop']
    
    # switch those positions in alls to TRUE
    alls[al_i, delpos] <- FALSE
    
    # done!
    
  }
  
  ## add some useful columns
  # grpmut
  alls <- alls %>%
    mutate(grpmut=del$grpmut, .before=1)
  
  
  ## order deletion alleles
  # can do by start position of deletion
  # split by '_' but best is to count from the end because group name might include _
  # it is 5th from the last
  # how many '_' splits?
  if(orderdels=='start') {
    nsplits <- unique(unlist(lapply(strsplit(alls$grpmut, '_'), length)))
    alls <- alls[order(as.numeric(strNthSplit(alls$grpmut, '_', nsplits-4))),]
    grpmutOrder <- alls$grpmut
    
  # can do by end position
  } else if(orderdels=='stop') {
    nsplits <- unique(unlist(lapply(strsplit(alls$grpmut, '_'), length)))
    alls <- alls[order(as.numeric(strNthSplit(alls$grpmut, '_', nsplits-3))),]
    grpmutOrder <- alls$grpmut
    
  # can do by total length of deletion
  } else if(orderdels=='length') {
    nsplits <- unique(unlist(lapply(strsplit(alls$grpmut, '_'), length)))
    alls <- alls[order(as.numeric(strNthSplit(alls$grpmut, '_', nsplits-2)), decreasing=TRUE),]
    grpmutOrder <- alls$grpmut
  
  # length by longest at the top
  } else if(orderdels=='lengthUp') {
    nsplits <- unique(unlist(lapply(strsplit(alls$grpmut, '_'), length)))
    alls <- alls[order(as.numeric(strNthSplit(alls$grpmut, '_', nsplits-2)), decreasing=FALSE),]
    grpmutOrder <- alls$grpmut
    
  # can do random
  } else if(orderdels=='random') {
    alls <- alls[sample(1:nrow(alls), replace=FALSE) ,]
    grpmutOrder <- alls$grpmut
  }
  
  # pivot
  alls_l <- alls %>%
    pivot_longer(-grpmut,
                 names_to='pos',
                 values_to='present')
  
  # recover grp
  alls_l <- alls_l %>%
    mutate(grp=strNthSplit(grpmut, '\\.', 1), .before=1)
  
  # positions
  alls_l$pos <- as.integer(alls_l$pos)
  alls_l <- alls_l[order(alls_l$pos),]
  
  # remove ni grp
  alls_l <- alls_l %>%
    filter(grp != 'ni')
  
  # crop a little, if given rangeStart/Stop
  # if not, leave everything
  
  # ! may have given bpSides
  # i.e. number of bp to keep on either side
  if(is.na(rangeStart) & is.na(rangeStop) & !is.na(bpSides)) {
    ## get cutpos from mut
    cutpos <- unique(mut$cutpos)
    # at present, can only work for one locus
    if(length(cutpos)>1) stop('\t \t \t \t >>> Error ggIGVdel: sorry, bpSides only works for one locus at a time for now.\n')
    rangeStart <- cutpos - bpSides
    rangeStop <- cutpos + bpSides
    cat('\t \t \t \t >>> cutpos is', cutpos, '; plotted', rangeStart, 'to', rangeStop, 'i.e.', bpSides, 'on either side.\n')
  }
  
  # now crop
  if(!is.na(rangeStart) & !is.na(rangeStop)) {
    alls_l <- alls_l %>%
      filter(pos %in% rangeStart : rangeStop)
  }
  
  ## order again using order we got above
  alls_l$grpmut <- factor(alls_l$grpmut, levels=grpmutOrder)
  
  cat('\t \t \t \t >>> Last position is', max(as.numeric(alls_l$pos)), '\n')
  cat('\t \t \t \t >>>', length(unique(alls_l$pos)), 'positions plotted \n')
  
  ## shortest & longest deletion
  nsplits <- unique(unlist(lapply(strsplit(as.character(alls_l$grpmut), '_'), length)))
  lens <- as.numeric(strNthSplit(alls_l$grpmut, '_', nsplits-2))
  cat('\t \t \t \t >>> Shortest deletion is', min(lens), 'bp. \n')
  cat('\t \t \t \t >>> Longest deletion is', max(lens), 'bp. \n')
  
  ### ligPos to plot?
  if(ligPos & is.na(forceligPos)) {
    ligxintercept <- unique(mut$rhapos)
    
    if(avgLigPos) {
      # in this case, maybe avgLigPos is ON
      # in which case we should use average of lig positions
      # (if there is only, will not change anything to do avergae)
      ligxintercept <- round(mean(ligxintercept))
    }
    
  } else if(ligPos & is.numeric(forceligPos)) {
    ligxintercept <- forceligPos
  }
  cat('\t \t \t \t >>> Ligation position at', ligxintercept, '\n')
   
  # now we can plot
  ggigv <- ggplot(alls_l, aes(x=pos, y=grpmut)) +
    geom_tile(aes(fill=present), colour='white', linewidth=0.01) +
    facet_grid(grp ~ . , scales='free', switch='y') + # switch puts title on the left
    # scale_fill_manual(values=c('#d53e4f', '#abdda4')) +
    scale_fill_manual(values=c('white', '#4d4d4d')) +
    theme_void() +
    theme(
      #axis.text.y=element_blank(),
      axis.text.x=element_blank(),
      axis.title.y=element_blank(),
      axis.title.x=element_blank(),
      legend.position='none',
      strip.text.y.left=element_text(size=9, angle=90, margin=margin(r=5)),
      plot.margin=margin(t=0, r=0, b=addBottomMargin, l=0)
    ) +
    {if(!titleOrNo) theme(strip.text.y.left=element_blank())} +
    {if(ligPos) geom_vline(xintercept=ligxintercept, linetype=5, linewidth=0.5, colour='#fcb505') } +
    geom_vline(xintercept=unique(mut$cutpos), linetype=5, linewidth=0.5, colour='#cb2a20')


  print(ggigv)
  
  if(endsWith(exportPath, '.pdf')) {
    ggsave(filename=exportPath, plot=ggigv, width=width, height=height, units='mm')
  } else if(endsWith(exportPath, '.png')) {
    ggsave(filename=exportPath, plot=ggigv, width=width, height=height, units='mm', dpi=1000)
    # output is ~ 0.5 Mb
  }
  
  ### return plot
  return(ggigv)
  
}