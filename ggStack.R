#####################################################
# ~ ZPRI: reads labels to stack plot ~
#
#
# Francois Kroll 2023
# francois@kroll.be
#####################################################

# v1
# v2: with classifyReads_one_v2.R which includes scaffold detection

# heavily based on function ggFrameshift created for previous projects

# ggFrameshift(...)

# given ampliCan results (config_summary.csv), plot frameshift plots per sample
# in the style of Figure 2A https://elifesciences.org/articles/59683/figures#content

# ! expects sample ID to be exactly: gene.locus_scr/ko_num
# e.g. psen2.1_ko_8
# e.g. psen2.3_scr_1

# arguments

# amplican = full path to main results from ampliCan (typically config_summary.csv)
# or R dataframe of it
# should guess if given the path or the dataframe directly

# onlygene = plot all samples of that gene
# default 'all' = plot all samples
# ! currently only supports all genes or one

# onlysource = only plot samples of a certain source
# 'source' column added to amplican results file
# default 'all' = plot all samples
# if no source column, just leave default

# covFilter = whether or not to filter low-coverage samples
# default = TRUE, i.e. yes filter out low coverage samples (see below)

# mincov = minimum coverage (above or equal to), in number of pairs of reads
# default = 30 x (i.e. 30x and above are okay)
# pairs of reads, so 1 pair = 1x, not 2x (i.e. what ampliCan reports, but not IGV)
# exclude the sample from the plot if coverage is lower

# mutated_col = colour to use for mutated proportion
# default = #5a6974, grey

# frameshift_col = colour to use for frameshifted proportion
# default = #f1876b, orange

# exportOrNo = whether to export the plot to pdf or no
# default = TRUE = yes, export to pdf

# width = width of the pdf plot (mm)
# by default = 159.1 mm = A4 minus margins

# height = height of the pdf plot (mm)
# by default = 82.03 mm = 1/3rd of A4 minus margins

# exportfull = full output path for plot export


# functions & packages ----------------------------------------------------

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(stringr)

sem <- function(x) sd(x)/sqrt(length(x))


ggStack <- function(rlab,
                    splitby='locus',
                    onlyrundate='all',
                    onlylocus='all',
                    onlygrp='all',
                    onlycat='all',
                    grporder=NA,
                    locusorder=NA,
                    catorder=NA,
                    catcols=NA,
                    mincov=NA,
                    ymax=1.0,
                    titleOrNo=TRUE,
                    xtextOrNo=TRUE,
                    ytextOrNo=TRUE,
                    ynameOrNo=TRUE,
                    legendOrNo=TRUE,
                    splcovOrNo=FALSE,
                    labelSamples=NA,
                    xgrid=FALSE,
                    titleSize=9,
                    ytitleSize=9,
                    ytextSize=7,
                    panelSpacing=NA,
                    exportOrNo=TRUE,
                    logOrNo=FALSE,
                    width=159.1,
                    height=82.03,
                    exportpath) {
  
  ### check exportpath is a pdf
  if(!is.na(exportpath)) {
    frmt <- substr(exportpath, start=nchar(exportpath)-2, stop=nchar(exportpath))
    if(frmt!='pdf')
      stop('\t \t \t \t >>> exportpath should end with .pdf.\n')
  }
  
  ### import reads labels, if necessary
  if(is.character(rlab)) { # if string then assume we are being given the path
    rlab <- read.csv(rlab)
  }
  # if not a string, we are given the dataframe so we do not have to import it
  
  
  # if user did not give catorder,
  # set defaults based on the mode
  if(is.na(catorder[1])) {
    
    if(unique(rlab$mode)=='precise') {
      catorder <- c('reference', 'scaffold', 'mutated', 'impure', 'pure')
      
    } else if(unique(rlab$mode)=='precise_simple') {
      catorder <- c('reference', 'mutated', 'edit')
      
    } else if(unique(rlab$mode)=='frameshift') {
      catorder <- c('reference', 'indel_inframe', 'indel_frameshift')
      
    } else if(unique(rlab$mode)=='insdel') {
      catorder <- c('noindel', 'both', 'insertion', 'deletion')
    }
  }

  
  # if user did not give catcols,
  # set defaults based on the mode
  if(is.na(catcols[1])) {
    
    if(unique(rlab$mode)=='precise') {
      catcols <- c('#aeb3b4', '#982150', '#fcb505', '#a3bbdb', '#417dcd')
      
    } else if(unique(rlab$mode)=='precise_simple') {
      catcols <- c('#aeb3b4', '#fcb505', '#417dcd')
      
    } else if(unique(rlab$mode)=='frameshift') {
      catcols <- c('#aeb3b4', '#f7b8b1', '#EE7163')
      
    } else if(unique(rlab$mode)=='insdel') {
      catcols <- c('#e7e8e9', '#982150', '#b0beab', '#c77f93')
    }
  }
  
  ### tally by sample and read category
  # e.g. for sample A01, we count all reference reads
  # force categories here, this (in combination with .drop=FALSE) will create e.g. pure edit for uninjected sample = 0
  rlab$cat <- factor(rlab$cat, levels=catorder)
  rtal <- rlab %>%
    group_by(sample, cat, .drop=FALSE) %>%
    tally(name='nreads')

  # then add back meta information
  rmeta <- rlab %>%
    distinct(sample, .keep_all=TRUE) %>%
    dplyr::select(sample, rundate, sid, locus, well, grp, splcov)
  
  rtal <- left_join(rmeta, rtal, by='sample')
  
  ### as check, separately, count total number of reads for each sample
  rtot <- rlab %>%
    group_by(sample, .drop=FALSE) %>%
    tally(name='ntot')
  # add the total to tally from above
  rtal <- right_join(rtal, rtot, by='sample')
  # check that splcov is always = ntot
  if(!identical(rtal$splcov, rtal$ntot))
    stop('\t \t \t \t >>> Sample\'s coverage calculated at the start is not equal to adding reads from the different categories.\n')
  # can delete ntot now
  rtal$ntot <- NULL
  
  ### calculate proportions
  rtal <- rtal %>%
    mutate(catpro=nreads/splcov)
  
  ### keep only rundate we want to plot
  if (onlyrundate[1] != 'all') {
    rtal <- rtal %>%
      filter(rundate %in% onlyrundate)
  }
  
  ### keep only locus we want to plot
  if (onlylocus[1] != 'all') {
    rtal <- rtal %>%
      filter(locus %in% onlylocus)
  }
  
  ### keep only grp we want to plot
  if (onlygrp[1] != 'all') {
    rtal <- rtal %>%
      filter(grp %in% onlygrp)
  }
  
  ### keep only categories we want to plot
  if (onlycat[1] != 'all') {
    rtal <- rtal %>%
      filter(cat %in% onlycat)
  }
  
  ### filter low coverage samples
  if (!is.na(mincov)) {
    
    # keep only one row for each sample
    rtalu <- distinct(rtal, sample, .keep_all=TRUE)
    
    # find positions of low coverage samples
    lcis <- which(rtalu$splcov < mincov) # low coverage indices
    
    # tell user we are throwing them out
    if (length(lcis) != 0) { # are there any samples to throw?
      for (i in lcis) {
        cat('\t \t \t \t >>> Excluding sample', as.character(rtalu[i, 'sample']),
            'because its coverage is', as.numeric(rtalu[i, 'splcov']), 'x.\n')
      }
    }
    
    # exclude them
    rtal <- rtal %>%
      filter(splcov >= mincov)
    
  }
  
  # pause: it is possible we excluded every sample!
  # in this case, warn the user and returns nothing
  if(nrow(rtal)==0) {
    cat('\t \t \t \t >>> STOP: no more reads left to count after excluding low-coverage samples. Stopping execution. \n')
    return()
  }
  
  # give minimum coverage to user
  cat('\t \t \t \t >>> Minimum coverage is', rtal[which.min(rtal$splcov), 'splcov'], 'x',
      'for sample', rtal[which.min(rtal$splcov), 'sample'], '\n')
  
  ### control order of loci, if given by the user
  if(!is.na(locusorder[1])) {
    rtal$locus <- factor(rtal$locus, levels=locusorder)
  }
  
  ### control order of groups, if given by the user
  if(!is.na(grporder[1])) {
    # keep only groups given by user, as if onlygrp was also given
    rtal <- rtal %>%
      filter(grp %in% grporder)
    rtal$grp <- factor(rtal$grp, levels=grporder)
  }
  
  ### control order of type in plot
  # categories will in stacked barplot from top to bottom
  rtal$cat <- factor(rtal$cat, levels=catorder)
  
  ### prepare coverage to be written like sample size on top of bar
  # can do it whether or not user asked for it
  covwrite <- rtal %>%
    distinct(sample, .keep_all=TRUE)
  
  ### if user asked to label some samples, prepare data now
  if(!is.na(labelSamples[1])) {
    labeldf <- data.frame(sample=labelSamples, label='*')
    sampledf <- rtal %>%
      distinct(sample, .keep_all=TRUE)
    labeldf <- left_join(sampledf, labeldf, by='sample')
  }

  
  ### plot
  ggstack <- ggplot(data=rtal, aes (x=sample, y=catpro, fill=cat)) +
    
    # split plot by...
    {if(splitby=='rundate') facet_grid(~rundate, scales='free_x', space='free')} +
    {if(splitby=='grp') facet_grid(~grp, scales='free_x', space='free')} +
    {if(splitby=='locus') facet_grid(~locus, scales='free_x', space='free')} +
    # space='free' lets subplots being different width,
    # which results in bar width staying equal regardless of number of samples for each grouping
    
    geom_col(width=0.8) +
    scale_fill_manual(drop=FALSE, values=catcols) +
    theme_minimal() +
    # write coverage on top of each bar?
    {if(splcovOrNo) geom_text(data=covwrite, aes(x=sample, label=splcov, y=1.0), size=2)} +
    {if(!is.na(labelSamples[1])) geom_text(data=labeldf, aes(x=sample, label=label, y=0.98), size=2.5)} +
    
    theme(

      panel.grid.minor=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_text(size=ytitleSize, margin=margin(t=0, r=-1, b=0, l=0)),
      axis.text.x=element_text(size=7, angle=90, hjust=1, vjust=0.5, margin=margin(t=-10, r=0, b=0, l=0)),
      axis.text.y=element_text(size=ytextSize, margin=margin(t=0, r=-1, b=0, l=0)),
      strip.text.x=element_text(size=titleSize)
      ) +
    coord_cartesian(ylim=c(0,ymax)) +
    scale_y_continuous(breaks=seq(0.0, 1.0, 0.1),
                       labels=seq(0, 100, 10)) +
    
    {if(ynameOrNo) ylab('% reads')} +
    {if(!ynameOrNo) ylab('')} +
    
    {if(!xtextOrNo) theme(axis.text.x=element_blank())} +
    {if(!ytextOrNo) theme(axis.text.y=element_blank())} +
    {if(!titleOrNo) theme(strip.text.x=element_blank())} +
    {if(!xgrid) theme(panel.grid.major.x=element_blank())} +
    {if(!legendOrNo) theme(legend.position='none')} +
    {if(!is.na(panelSpacing)) theme(panel.spacing = unit(panelSpacing, 'pt'))}
  
  print(ggstack)
  
  
  ### export plot
  if (!is.na(exportpath)) {
    cat('\t \t \t \t >>> Exporting to', exportpath, '.\n')
    ggsave(exportpath, plot=ggstack, width=width, height=height, unit='mm')
  } else {
    cat('\t \t \t \t >>> Export is off.\n')
  }
  
  
  ### prepare & export log about samples plotted
  if(logOrNo) {
    logspl <- rtal %>%
      distinct(sample, .keep_all=TRUE)
    logspl$cat <- NULL
    logspl$catpro <- NULL
    # sort by whatever factor user gave as 'splitby'
    logspl[[splitby]] <- as.factor(logspl[[splitby]])
    logspl <- logspl[order(logspl[[splitby]]) , ]
    
    # prepare export
    logtxtPath <- paste0(dirname(exportpath), '/',
                         tools::file_path_sans_ext(basename(exportpath)), '_LOG.csv')
    write.csv(logspl, logtxtPath, row.names=FALSE)
  }
  
  
  ### return processed tallies
  invisible(rtal)
  
} # closes the function