####################################################################
# Author: Jafta Stomp
# Date: 28-02-2018
# Description: 
#   script for building circos up plot crammed in 1 reusable function
####################################################################

circler <- function(g, d3Gos, df, go_groups,tests){
  # set colors for sectors
  col1 <- colorRampPalette(c("magenta", "dodgerblue3", "gray87", "tan1", "firebrick3")) (n=length(d3Gos))
  names(col1) = d3Gos
  col2 = brewer.pal(3, "Set1")
  names(col2) = tests
  
  df[[1]] = as.character(df[[1]])
  df[[2]] = as.character(df[[2]])
  
  # sector set-up
  sector = NULL
  sector_xlim = NULL
  for(t in unique(df[[1]])) {
    sector = c(sector, t)
    sector_xlim = rbind(sector_xlim, c(0, sum(df[df[[1]] == t, 3])))
  }
  for(t in unique(df[[2]])) {
    sector = c(sector, t)
    sector_xlim = rbind(sector_xlim, c(0, sum(df[df[[2]] == t, 4])))
  }
  
  # Initiate the circos frame
  circos.par(cell.padding = c(0, 0, 0, 0), start.degree = 270, gap.degree = c(rep(1,length(tests)-1),10,rep(1,length(d3Gos)-1),10))
  
  #  Initialize circos using factors, 3 comparisons and x depth 3 GOs also set sector width with a total sum of 2
  circos.initialize(factors = factor(sector, levels = sector), xlim = sector_xlim,
                    sector.width = c(rep(1/3,3), c(rep(1/length(d3Gos),length(d3Gos)))))
  
  # First track on the right side of the plot containing the GO group names
  circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    if(sector.index %in% sector[4:length(sector)]) {
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      l = unique(g[[3]]) == sector.index
      x = seq(0, by = 3, length = sum(l))
      x = x + mean(xlim) - mean(x)
      circos.rect(xlim[1], ylim[1], xlim[2], ylim[2], col = col1[sector.index], border = NA)
      circos.text(x, rep(0, sum(l)), sector.index, col = "black", facing = "bending.inside", adj = c(0.5, -1), cex = 0.7)
      
    }
  }, bg.border = NA, track.height = 0.08)
  
# Third track full circle, right side GO group names again and on the left side the 3 tests
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.6, facing = "bending.inside")
}, track.height = 0.05, bg.border = NA)
  
  # Fourth track, this will be filles with the proportions of either genes or gos from left to right or right to left
  circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  }, track.height = 0.02, bg.col = c(col2, col1), track.margin = c(0, 0.01))
  circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  }, track.height = 0.02, track.margin = c(0, 0.01))
  
  # get accumulative data for the final part of the plot (linking)
  accum_gos = sapply(d3Gos, function(x) get.cell.meta.data("xrange", sector.index = x)); names(accum_gos) = d3Gos
  accum_genes = sapply(tests, function(x) get.cell.meta.data("xrange", sector.index = x)); names(accum_genes) = tests
  
  for(i in seq_len(nrow(df))) {
    # Make a link between tests and GO groups
    circos.link(df[i,1], c(accum_genes[df[i,1]], accum_genes[df[i,1]] - df[i, 3]),
                df[i,2], c(accum_gos[df[i,2]], accum_gos[df[i,2]] - df[i, 4]),
                col = paste0(col1[df[i,2]], "80"), border = NA)
    
    # Add colors for the fourth track proportions.
    circos.rect(accum_genes[df[i,1]], 0, accum_genes[df[i,1]] - df[i, 3], 1, sector.index = df[i,1], col = col1[df[i,2]])
    circos.rect(accum_gos[df[i,2]], 0, accum_gos[df[i,2]] - df[i, 4], 1, sector.index = df[i,2], col = col2[df[i,1]])
  
    accum_genes[df[i,1]] = accum_genes[df[i,1]] - df[i, 3]
    accum_gos[df[i,2]] = accum_gos[df[i,2]] - df[i, 4]
  }
  
  #clear the plot.
  circos.clear()
}
