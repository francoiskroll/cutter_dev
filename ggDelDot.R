#####################################################
# ~ ZPRI/utilities: plot deletion asymmetry vs. length ~
# function ggDelDot
#
# Francois Kroll 2025
# francois@kroll.be
#####################################################

source('~/Dropbox/cutter/detectMHdel.R')

# take mut mutation table
# to quantify bias, we calculate
# leftDel = | start - cutpos | > extent on left side
# rightDel = | stop - cutpos | > extent on right side

# then can divide left side / right side = biasDiv
# will give
# 1.0 = symmetrical
# < 1.0 = more left
# > 1.0 = more right

# can do right side - left side = biasDif
# so gives bp "protruding" on one side
# 0 = symmetrical, both exactly same length
# neg = biased left (left is bigger)
# pos = biased right (right is bigger)
# then divide by length of deletion, gives extra on one side in % = biasDifPro
# can keep sign
# or not = biasDifProAbs

ggDelDot <- function(mut,
                     colMHdel='#c77f93',
                     minMHlen=1,
                     xtextOrNo=TRUE,
                     xnameOrNo=TRUE,
                     ynameOrNo=TRUE,
                     legendOrNo=TRUE,
                     marginTB=NA,
                     exportPath=NA,
                     width=60,
                     height=50) {
  
  ### keep only unique deletions
  mubia <- mut %>%
    filter(type=='del') %>%
    distinct(mutid, .keep_all=TRUE)
  
  ### run detectMHdel
  if(!is.na(minMHlen)) {
    mubia <- detectMHdel(mut=mubia,
                         minMHlen=minMHlen)
    # add a boolean column for MH or not
    mubia <- mubia %>%
      mutate(MHthere=!is.na(MHbp), .before='MHside')
    
  }
  
  ### calculate bias as explained above
  # leftDel #
  mubia <- mubia %>%
    mutate(leftDel=cutpos-start)
  # 0 means first deleted nucleotide is cutpos
  # if pos: cutpos is after start (e.g. cutpos #50, start #45), as expected from del that removes cutpos
  # if neg: deletion starts *after* cutpos (e.g. cutpos #50, start #55), so deletion is entirely right of cutpos
  # > need to correct and set to 0, i.e. 0 bp on left of cutpos
  mubia[which(mubia$leftDel<0), 'leftDel'] <- 0
  # e.g. tyr base editing 166/228 73%, which makes sense from plot
  
  # rightDel #
  mubia <- mubia %>%
    mutate(rightDel=stop-cutpos)
  # 0 means last deleted nucleotide is cutpos
  # if pos: cutpos is before stop (e.g. cutpos #50, stop #55), as expected from del that removes cutpos
  # if neg: deletion stops *before* cutpos (e.g. cutpos #50, stop #45), i.e. deletion is entirely left of cutpos
  # > need to correct and set to 0, i.e. 0 bp on right of cutpos
  mubia[which(mubia$rightDel<0), 'rightDel'] <- 0
  # e.g. tyr base editing 4/228, which makes sense from plot
  
  # ! also cannot write leftDel or rightDel that is greater than length of deletion
  # does not make sense
  # maximum bias on one side is size of deletion
  # below: whenever leftDel or rightDel is greater than bp, just edit to bp
  mubia[which(mubia$leftDel > mubia$bp), 'leftDel'] <- mubia[which(mubia$leftDel > mubia$bp), 'bp']
  mubia[which(mubia$rightDel > mubia$bp), 'rightDel'] <- mubia[which(mubia$rightDel > mubia$bp), 'bp']
  
  # in summary, leftDel / rightDel never negative, as expected
  
  # now normalised metrics
  mubia <- mubia %>%
    mutate(biasDiv=leftDel/rightDel) %>%
    # ! this might generate 0 or Inf, will need to think about it if use it
    
    mutate(biasDif=rightDel-leftDel) %>%
    mutate(biasDifPro=biasDif/bp) %>%
    # +1.0 is deletion entirely on right, e.g. 50 - 0 = 50 >> 50 / 50 bp = +1.0
    # positive is biased right
    
    # -1.0 is deletion entirely on left, e.g. 0 - 50 = -50 >> -50 / 50 bp = -1.0
    # negative is biased left
    
    # 0.0 is perfectly symmetrical, same length on either side
    # (to be precise it is +1 bp on one side because cutpos is nucleotide just before cut, but will not matter much)
    
    # can also make asymmetry absolute measure in case useful
    mutate(biasDifProAbs=abs(biasDifPro))
  
  ### prepare text to write
  # number of deletions on left or right / total number of deletions
  # then line below: %
  # lefttxt <- paste0(sum(mubia$biasDifPro<0), '/', nrow(mubia), '\n',
  #                   round((sum(mubia$biasDifPro<0) / nrow(mubia))*100), '%')
  # 
  # righttxt <- paste0(sum(mubia$biasDifPro>0), '/', nrow(mubia), '\n',
  #                    round((sum(mubia$biasDifPro>0) / nrow(mubia))*100), '%')
  
  # alternative version in one line instead
  lefttxt <- paste0(sum(mubia$biasDifPro<0), '/', nrow(mubia),
                    ' ', round((sum(mubia$biasDifPro<0) / nrow(mubia))*100), '%')
  
  righttxt <- paste0(sum(mubia$biasDifPro>0), '/', nrow(mubia),
                     ' ', round((sum(mubia$biasDifPro>0) / nrow(mubia))*100), '%')
  
  txtdf <- data.frame(biasDifPro=c(-0.5, 0.5),
                      txt=c(lefttxt, righttxt),
                      MHthere=NA) # complains without
  
  ### plot
  if(is.na(colMHdel)) {
    ggdeldot <- ggplot(mubia, aes(x=biasDifPro, y=bp)) +
      geom_vline(xintercept=0, linetype=2) +
      geom_point(shape=21, size=1, colour='white', fill='#4d4d4d', stroke=0.2) +
      geom_text(data=txtdf, aes(x=biasDifPro, label=txt, y=txtH), size=1.9)
  } else {
    ggdeldot <- ggplot(mubia, aes(x=biasDifPro, y=bp, fill=MHthere)) +
      geom_vline(xintercept=0, linetype=2) +
      geom_point(shape=21, size=1, colour='white', stroke=0.2) +
      geom_text(data=txtdf, aes(x=biasDifPro, label=txt, y=txtH), size=1.9) +
      scale_fill_manual(values=c('#5a6974', colMHdel))
  }
  
  
  ## height at which we should write
  txtH <- max(mubia$bp) + 3
  
  ggdeldot <- ggdeldot +
    theme_minimal() +
    coord_cartesian(xlim=c(-1, 1),
                    ylim=c(0, txtH+1)) +
    
    theme(
      axis.text.x=element_text(size=7),
      axis.text.y=element_text(size=7, margin=margin(r=-1)),
      panel.grid.minor.x=element_blank(),
      panel.grid.minor.y=element_blank()
    ) +
    
    {if(!xtextOrNo) theme(axis.text.x=element_blank())} +
    
    {if(!xnameOrNo) theme(axis.title.x=element_blank())} +
    {if(xnameOrNo) theme(axis.title.x=element_text(size=8))} +
    
    {if(!ynameOrNo) theme(axis.title.y=element_blank())} +
    {if(ynameOrNo) theme(axis.title.y=element_text(size=8))} +
    
    {if(!legendOrNo) theme(legend.position='none')} +
    
    {if(!is.na(marginTB)) theme(plot.margin=margin(t=marginTB, r=5.5, b=marginTB, l=5.5))} +
    
    scale_x_continuous(breaks=c(-1.0, -0.5, 0, 0.5, 1.0),
                       labels=c('-1.0', '', '0', '', '+1.0')) +
    
    xlab('asymmetry') +
    ylab('deletion length (bp)')
  
  print(ggdeldot)
  
  if(!is.na(exportPath)) {
    ggsave(exportPath, ggdeldot, width=width, height=height, units='mm')
  }
  
  return(ggdeldot)
  
}