## this script includes the functions to plot mirrored Manhattan plots


## pos_to_cumpos
## this function convert chromosome-position pairs to continuous positions
## chr: a vector with the length of
## pos: a vector with the same length of chr
## gap_bp: numeric, determine when to decrease distance between two SNPs that are far away
pos_to_cumpos <- function(chr, pos, gap_bp = 5e6) {
  cumpos = pos
  current = pos[1]
  for (i in 2:length(pos)) {
    if (chr[i] != chr[i-1] | pos[i] - pos[i-1] > gap_bp) {
      cumpos[i:length(cumpos)] = cumpos[i:length(cumpos)] - cumpos[i] + cumpos[i-1] + gap_bp
    }
  }
  return(cumpos)
}


## point_cols
points_cols <- function(chr, POS, up = c('#7876B1FF', '#EE4C97FF'),
                        bottom = c('grey', 'skyblue')) {
  chr_cut <- c()
  colors <- c()
  colors_val <- c()
  for (ch in unique(chr)) {
    i=1
    tmp <- POS[chr == ch]
    cut <- (min(tmp) + max(tmp))/2
    chr_cut <- c(chr_cut, cut)
    if (as.numeric(ch)%%2 == 1) {
      co = up[1]
      co_val = bottom[1]
    }
    if (as.numeric(ch)%%2 == 0) {
      co = up[2]
      co_val = bottom[2]
    }
    col <- rep(co, times = length(tmp))
    col_val <- rep(co_val, times = length(tmp))
    colors <- c(colors, col)
    colors_val <- c(colors_val, col_val)
    i = i+1
  }
  return(list(chr_cut = chr_cut, colors = colors, colors_val = colors_val))
}




colors <- function(chr, pos, gap_bp = 5e6, up.col = c('#7876B1FF', '#EE4C97FF'),
                   btm.col = c('grey', 'skyblue')) {
  cumpos <- pos_to_cumpos(chr, pos, gap_bp)
  res <- points_cols(chr, cumpos, up = up.col, bottom = btm.col)
  #res <- c(res, POS = cumpos)
  return(res)
}







## scale the val columns
mirror_val <- function(x, sca=20000, dis=0.025) {
  xx <- -(x - max(x, na.rm = TRUE))
  xx <- xx - max(xx, na.rm = TRUE)
  xx <- xx/sca - dis
  return(xx)
}




## scale values into [0,1]

To1 <- function(x) {
  x[x>300] <- 300
  max.x <- max(x, na.rm = TRUE)
  #max.x <- ifelse(max.x > 300, 300, max.x)
  min.x <- min(x, na.rm = TRUE)
  xx <- (x-min.x)/(max.x-min.x)
  return(xx)
}





runrun <- function(dd, chr, pos, up, bottom, chr.lab = NA, up.col = c('#7876B1FF', '#EE4C97FF'), btm.col = c('grey', 'skyblue'), file) {
  ylims <- c(-1.1, 1)
  dd <- dd[order(dd[, chr], dd[, pos], decreasing = F),]

  if(is.na(chr.lab) == TRUE) {
    chr_lab <- unique(dd[,chr])
  } else {
    chr_lab <- chr.lab
  }

  #chr_lab <- ifelse(is.na(chr.lab), c(unique(dd[,chr])), chr.lab)

  POS <- pos_to_cumpos(dd[, chr], dd[, pos])
  mycolors <- colors(dd[, chr], dd[, pos], up.col = up.col, btm.col = btm.col)


  upval <- To1(dd[, up])
  upmin <- round(min(dd[, up], na.rm = TRUE), digits = 4)
  upmax <- round(max(dd[, up], na.rm = TRUE), digits = 4)
  upmax <- ifelse(upmax > 300, 300, upmax)
  upbrk <- seq(0,1,.2)
  uplab <- round(seq(upmin, upmax, length = 6), digits = 2)


  btmval <- -To1(dd[, bottom])-0.1
  btmin <- round(min(dd[, bottom], na.rm = TRUE), digits = 4)
  btmax <- round(max(dd[, bottom], na.rm = TRUE), digits = 4)
  btmax <- ifelse(btmax > 300, 300, btmax)
  btbrk <- seq(0,-1,-.2)-0.1
  btlab <- round(seq(btmin, btmax, length = 6), digits = 2)

  #sigline <- -To1(-log10(5e-8)) -0.1

  png(file, width = 240, height = 120, units = 'mm', res = 300, pointsize = 8)
  par(mar = c(1,3,1,1), oma = c(3,4,1,1))
  plot(POS, upval, pch = 16, cex = .4, axes = F, col = mycolors$colors, ylim = ylims, xlab = '', ylab = '')
  par(new = T)
  plot(POS, btmval, pch = 16, cex = .4, axes = F, col = mycolors$colors_val, ylim = ylims, xlab = '', ylab = '')
  axis(1, at = mycolors$chr_cut, las = 1, cex.axis = 1.1, labels = chr_lab)
  axis(2, at = upbrk, labels = uplab, las = 1, cex.axis = 1.2)
  axis(2, at = btbrk, labels = btlab, las = 1, cex.axis = 1.2)
  mtext('Chromosome', side=1, line=2, cex=1.1, col="black", outer=TRUE, padj = -.8)
  mtext(expression(tau), side=2, line=2, cex=2, col="black", outer=TRUE, adj = 0.78, padj = .7)
  mtext(expression(-log[10]~italic(P)), side=2, line=2, cex=1, col="black", outer=TRUE, adj = 0.33, padj = .45)
  abline(h = sigline, col=2)
  dev.off()
}











