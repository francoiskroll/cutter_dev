#####################################################
# ~ cutter: plot stacked barplot of deletions with microhomology ~
#
#
# Francois Kroll 2025
# francois@kroll.be
#####################################################

library(colorspace)
library(ggplot2)
library(dplyr)
library(tidyr)


# ggMHdel -----------------------------------------------------------------
# colours correspond to length of MH

# ! colourManual: need to give colours as 0 bp, then reverse order
# e.g. 0 bp, 4 bp, 3 bp, 2 bp, 1 bp

# colourManual:
# paper used colourManual=c('#E9EBEC', '#C26682', '#C77F93', '#CD96A4', '#FFE4EA'),
# here are 6 darker if need more:
# #C04871" "#AB3861" "#962952" "#821843" "#6E0335" "#560027"

ggMHdel <- function(mut,
                    min_del_nreads=50, # arbitrary
                    colourLight='#f1b9c7',
                    colourManual=NA,
                    splitby='grp',
                    onlygrp=NA,
                    grporder=NA,
                    legendOrNo=TRUE,
                    titleOrNo=TRUE,
                    xtextOrNo=TRUE,
                    ytextOrNo=TRUE,
                    ynameOrNo=TRUE,
                    panelSpacing=NA,
                    exportpath,
                    width=110,
                    height=65) {
  
  ### if onlygrp is given, keep only those groups
  if(!is.na(onlygrp[1])) {
    cat('\t \t \t \t >>> ggMHdel: not plotting groups', unique(mut$grp)[!unique(mut$grp) %in% onlygrp], '\n')
    mut <- mut %>%
      filter(grp %in% onlygrp)
  }
  
  ### keep only deletions
  del <- mut %>%
    filter(type=='del')
  
  ### ! one read can have more than one deletions
  # there is no great way to deal with those, best is probably to remove them
  # otherwise we count one read multiple times below, so axis as % of reads is not perfectly correct
  
  # add unique read id
  del <- del %>%
    mutate(urid=paste(sample, rid, sep='_'))
  
  readstodel <- del %>%
    group_by(urid) %>%
    tally(name='ndel') %>%
    filter(ndel>1)
  cat('\t \t \t \t >>>', nrow(readstodel), 'reads, out of', length(unique(del$urid)), 'reads, have more than one deletion, we will delete them.\n')
  
  del <- del %>%
    filter(! urid %in% readstodel$urid)
  
  # now, each read mentioned in del has one deletion, easy
  # also the case that each read is mentioned only once, check:
  if(sum(duplicated(del$urid))>0) stop('\t \t \t \t Error ggMHdel: some reads are still mentioned more than once after removing reads with more than one deletion, something is not right!\n')
  
  # remove urid column, was not there originally
  del$urid <- NULL
  
  ### tally, by sample & length of MH (MHbp), number of reads
  mhtal <- del %>%
    group_by(sample, rundate, locus, sid, grp, well, splcov, MHbp,
             .drop=FALSE) %>%
    tally(name='nreads')
  
  # replace NA in MHbp by 0, easier to deal with later
  mhtal[which(is.na(mhtal$MHbp)), 'MHbp'] <- 0
  
  # same tally but total number of reads with a deletion
  deltal <- del %>%
    group_by(sample,
             .drop=FALSE) %>%
    tally(name='del_nreads')
  
  # merge both
  mhtal <- left_join(mhtal, deltal, by='sample')
  
  # calculate frequencies & plot
  mhtal <- mhtal %>%
    mutate(catpro=nreads/del_nreads)
  
  # exclude samples with less than min_del_nreads reads with deletion (arbitrary)
  mhtal <- mhtal %>%
    filter(del_nreads >= min_del_nreads)
  
  # set order of MHbp
  # we want 0 first, then decreasing
  # e.g. 0, 6, 4, 3, 2, 1
  MHbp_ord <- sort(unique(mhtal$MHbp), decreasing=TRUE)
  # remove 0 if present
  MHbp_ord <- MHbp_ord[-which(MHbp_ord==0)]
  # and add it as first
  MHbp_ord <- c(0, MHbp_ord)
  mhtal$MHbp <- factor(mhtal$MHbp, levels=MHbp_ord)
  
  # prepare the colours
  # either user user gave us the lightest colour >> then we make it darker for each additional bp of MH
  # or user gave us all colours with colourManual
  if(!is.na(colourLight) & is.na(colourManual[1])) {
    MHbp_ord_no0 <- MHbp_ord[-which(MHbp_ord==0)]
    catcols <- sapply(MHbp_ord_no0, function(bp) {
      darken(col=colourLight, amount=0.08*bp)
    })
    # add colour for 0
    catcols <- c('#d7d9da', catcols)
  } else if(is.na(colourLight) & !is.na(colourManual[1])) {
    catcols <- colourManual
  } else {
    stop('\t \t \t \t >>> Error: set colourLight as NA and give colourManual, or vice-versa.\n')
  }
  cat('\t \t \t \t >>> Colours starting from 0 bp MH are: ', catcols, '\n')
  
  # set the order of the groups (plot's facet)
  if(!is.na(grporder[1])) {
    mhtal$grp <- factor(mhtal$grp, levels=grporder)
  }
  
  #catcols <- c('#5a6974', '#c54867', '#db5072', '#e2738e', '#e996aa', '#f1b9c7')

  ggMhbp <- ggplot(mhtal, aes(x=sample, y=catpro, fill=MHbp)) +
    
    {if(splitby=='rundate') facet_grid(~rundate, scales='free_x', space='free')} +
    {if(splitby=='grp') facet_grid(~grp, scales='free_x', space='free')} +
    {if(splitby=='locus') facet_grid(~locus, scales='free_x', space='free')} +
    
    geom_col(width=0.8) +
    scale_fill_manual(drop=FALSE, values=catcols) +
    theme_minimal() +
    theme(
      panel.grid.minor=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_text(size=9),
      # axis.text.x=element_text(size=7, angle=90, hjust=1, vjust=0.5, margin=margin(t=-10, r=0, b=0, l=0)),
      axis.text.y=element_text(size=7, margin=margin(t=0, r=-1, b=0, l=0)),
      axis.text.x=element_text(size=7, angle=90, hjust=1, vjust=0.5, margin=margin(t=-10, r=0, b=0, l=0))) +
    scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1.0),
                       labels=c('0', '25', '50', '75', '100')) +
    
    {if(!legendOrNo) theme(legend.position='none')} +
    {if(ynameOrNo) ylab('% of reads with deletion')} +
    {if(!ynameOrNo) ylab('')} +
    
    {if(!xtextOrNo) theme(axis.text.x=element_blank())} +
    {if(!ytextOrNo) theme(axis.text.y=element_blank())} +
    {if(!titleOrNo) theme(strip.text.x=element_blank())} +
    {if(!is.na(panelSpacing)) theme(panel.spacing = unit(panelSpacing, 'pt'))}
  
  print(ggMhbp)
  
  ### report to user which samples got included
  cat('\t \t \t \t', length(unique(mhtal$sample)), 'samples plotted, including potential simulated sample.\n')
  cat('\t \t \t \t groups plotted:', unique(mhtal$grp), '\n')
  
  if(!is.na(exportpath)) {
    ggsave(exportpath, ggMhbp, width=width, height=height, units='mm')
  }

  return(mhtal)
  
}
