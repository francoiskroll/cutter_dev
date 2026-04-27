#####################################################
# ~ cutter: ggTemplIns ~
#
# ... ... ...
# Francois Kroll 2025
# francois@kroll.be
#####################################################

library(colorspace)
library(ggplot2)
library(dplyr)
library(tidyr)

library(forcats)
library(scales)


# ggTemplIns --------------------------------------------------------------
# colours correspond to length of LCS

ggTemplIns <- function(mut,
                       min_ins_nreads=50, # arbitrary
                       colourLow='#e4eee0',
                       colourHigh='#304528',
                       grporder=NA,
                       legendOrNo=TRUE,
                       titleOrNo=TRUE,
                       xtextOrNo=TRUE,
                       ytextOrNo=TRUE,
                       ynameOrNo=TRUE,
                       exportpath,
                       width=110,
                       height=65) {
  
  ### keep only insertions
  ins <- mut %>%
    filter(type=='ins')
  
  ### ! one read can have more than one insertion
  # but usually it is the same newly synthesised sequence, and therefore the same output in terms of LCS length
  # so we do not necessairly have to delete them
  
  # add unique read id
  ins <- ins %>%
    mutate(urid=paste(sample, rid, sep='_'))
  
  # these reads have more than one insertion
  nins_perread <- ins %>%
    group_by(urid) %>%
    tally(name='nins') %>%
    filter(nins>1)
  
  # look if they have different newly synthesised sequences or not
  check <- ins %>%
    filter(urid %in% nins_perread) %>%
    group_by(urid) %>%
    summarise(nnewlyseq=n_distinct(newlyseq))
  # for now, will throw an error if not the case
  # then will deal with it if it ever happens
  if(nrow(check)>0)
    stop('\t \t \t \t Error ggTemplIns: some reads have multiple insertions which do not have the same newly synthesised sequence.\n')
  
  # we still need to deal with reads which have multiple insertions
  # we can just remove duplicated read IDs
  if(sum(duplicated(ins$urid))>0) {
    ins <- ins[-which(duplicated(ins$urid)),]
  }
  
  # check to be safe
  nins_perread <- ins %>%
    group_by(urid) %>%
    tally(name='nins') %>%
    filter(nins>1)
  if(nrow(nins_perread>0))
    stop('\t \t \t \t Error ggTemplIns: some reads still have multiple insertions even after removing duplicated urids.\n')
  # another way to check, each read is mentioned only once:
  if(sum(duplicated(ins$urid))>0) stop('\t \t \t \t Error ggTemplIns: some reads are still mentioned more than once after removing duplicated urids, something is not right!\n')
  
  # remove urid column, was not there originally
  ins$urid <- NULL
  
  ### tally, by sample & length of LCS (lcbp), number of reads
  lctal <- ins %>%
    group_by(sample, rundate, locus, sid, grp, well, splcov, lcbp,
             .drop=FALSE) %>%
    tally(name='nreads')
  
  # same tally but total number of reads with insertion
  instal <- ins %>%
    group_by(sample,
             .drop=FALSE) %>%
    tally(name='ins_nreads')
  
  # merge both
  lctal <- left_join(lctal, instal, by='sample')
  
  # calculate frequencies & plot
  lctal <- lctal %>%
    mutate(catpro=nreads/ins_nreads)
  
  # exclude samples with less than min_ins_nreads reads with insertion
  lctal <- lctal %>%
    filter(ins_nreads >= min_ins_nreads)
  
  # replace lcbp by NA
  # to be clearer it is not part of the gradient
  lctal[which(lctal$lcbp==0), 'lcbp'] <- NA
  
  # turn lcbp into integer
  lctal$lcbp <- as.integer(lctal$lcbp)
  
  # have a new column lcbp_gg
  # to control levels just for ggplot2
  lctal <- lctal %>%
    mutate(
      lcbp_gg = fct_na_value_to_level(factor(lcbp), level='NA'),
    )
  
  # prepare the levels
  # slightly complicated because we want NA on top, then higher to lower values (going top to bottom)
  lvls <- levels(lctal$lcbp_gg)
  num_lvls <- sort(as.numeric(lvls[lvls!='NA']), decreasing=TRUE)
  lctal <- lctal %>%
    mutate(
      lcbp_gg=factor(lcbp_gg, levels=c('NA', as.character(num_lvls)))
    )
  
  # prepare colour gradient
  # again, slightly complicated as we need one colour for NA, then gradient for the other colours
  # create a gradient function (e.g. using blue to green)
  gradient_fn <- scales::gradient_n_pal(c(colourLow, colourHigh))
  
  # map levels to colors
  fillCols <- c(
    'NA' = '#E9EBEC',
    setNames(gradient_fn(rescale(num_lvls)), as.character(num_lvls))
  )
  
  # set the order of the groups (plot's facet)
  if(!is.na(grporder[1])) {
    lctal$grp <- factor(lctal$grp, levels=grporder)
  }
  
  ### ready to plot
  ggInsbp <- ggplot(lctal, aes(x=sample, y=catpro, fill=lcbp_gg)) +
    facet_grid(~grp, scales='free_x', space='free') +
    geom_col(width=0.8) +
    scale_fill_manual(values=fillCols) +
    theme_minimal() +
    theme(
      panel.grid.minor=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_text(size=9),
      axis.text.x=element_text(size=7, angle=90, hjust=1, vjust=0.5, margin=margin(t=-10, r=0, b=0, l=0)),
      axis.text.y=element_text(size=7, margin=margin(t=0, r=-1, b=0, l=0))) +
    scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1.0),
                       labels=c('0', '25', '50', '75', '100')) +
    {if(!legendOrNo) theme(legend.position='none')} +
    {if(ynameOrNo) ylab('% of reads with insertion')} +
    {if(!ynameOrNo) ylab('')} +
    
    {if(!xtextOrNo) theme(axis.text.x=element_blank())} +
    {if(!ytextOrNo) theme(axis.text.y=element_blank())} +
    {if(!titleOrNo) theme(strip.text.x=element_blank())}
  
  print(ggInsbp)
  
  if(!is.na(exportpath)) {
    ggsave(exportpath, ggInsbp, width=width, height=height, units='mm')
  }

  return(lctal)
  
}
