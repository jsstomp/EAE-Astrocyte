####################################################################
# Author: Jafta Stomp
# Date: 26-03-2018
# Description: 
#   This script makes a circos plot using DEA/GO results
####################################################################
library(circlize)

set.seed(999)
n = 1000
df = data.frame(factors = sample(letters[1:8], n, replace = TRUE),
                x = rnorm(n), y = runif(n))
df <- data.frame(factors = c(sample(letters[1],500, replace=T), sample(letters[2:8], 500, replace=T))
                 ,x = rnorm(n), y=runif(n))


circos.par("track.height" = 0.1)
circos.initialize(factors = df$factors, x = df$x, sector.width = c(7,1,1,1,1,1,1,1))

circos.track(factors = df$factors, y = df$y,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(6, "mm"), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             })
col = rep(c("#FF0000", "#00FF00"), 4)
circos.trackPoints(df$factors, df$x, df$y, col = col, pch = 16, cex = 0.5)
circos.text(-1, -0.5, "text", sector.index = "a", track.index = 1)

bgcol = rep(c("#EFEFEF", "#CCCCCC"), 4)
circos.trackHist(df$factors, df$x, bin.size = 0.2, bg.col = bgcol, col = NA)

circos.track(factors = df$factors, x = df$x, y = df$y,
             panel.fun = function(x, y) {
               ind = sample(length(x), 10)
               x2 = x[ind]
               y2 = y[ind]
               od = order(x2)
               circos.lines(x2[od], y2[od])
             })

circos.update(sector.index = "a", track.index = 2, 
              bg.col = "#FF8080", bg.border = "black")
circos.points(x = -2:2, y = rep(0.5, 5), col = "white")
circos.text(CELL_META$xcenter, CELL_META$ycenter, "updated", col = "white")

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  breaks = seq(xlim[1], xlim[2], by = 0.1)
  n_breaks = length(breaks)
  circos.rect(breaks[-n_breaks], rep(ylim[1], n_breaks - 1),
              breaks[-1], rep(ylim[2], n_breaks - 1),
              col = rand_color(n_breaks), border = NA)
})

circos.link("a", 0, "b", 0, h = 0.4)
circos.link("c", c(-0.5, 0.5), "d", c(-0.5,0.5), col = "red",
            border = "blue", h = 0.2)
circos.link("e", 0, "g", c(-1,1), col = "green", border = "black", lwd = 2, lty = 2)
circos.link("c", c(0,1), "h", c(-1,1), col = "yellow", lwd = 2, lty = 2)

circos.clear()

##########################################################################################################

#op = par(no.readonly = TRUE)

#par(mar = c(2, 2, 2, 2), mfrow = c(1, 3))

factors = letters[1:4]
factors = c("HBA","HBAG","SCA",letters[1:26])
factors1 = c(factors[1:3],paste(factors[1:3],".1"))
factors2 = sample(rownames(countData),200)
# factors2 = c(genes, paste(genes,".1"))
circos.par("canvas.xlim" = c(-1, 1.5), "canvas.ylim" = c(-1, 1.5), start.degree = -45)
circos.initialize(factors = factors1, xlim = c(0, 1))
circos.trackPlotRegion(ylim = c(0, 1), bg.col = NA, bg.border = NA)
circos.updatePlotRegion(sector.index = factors[1])
circos.text(0.5, 0.5, factors[1], niceFacing = T)
circos.updatePlotRegion(sector.index = factors[2])
circos.text(0.5, 0.5, factors[2], niceFacing = T)
circos.updatePlotRegion(sector.index = factors[3])
circos.text(0.5, 0.5, factors[3], niceFacing = T)

circos.clear()
box()
axis(side = 1)
axis(side = 2)


circos.clear()

par(new = TRUE)
circos.par("canvas.xlim" = c(-1.5, 1), "canvas.ylim" = c(-1.5, 1), start.degree = -45)
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(factors = factors2, xlim = c(0, 1))
circos.trackPlotRegion(ylim = c(0, 1), bg.col = NA, bg.border = NA)
for(i in factors2[101:length(factors2)]){
  circos.updatePlotRegion(sector.index = i)
  # circos.text(0.5, 0.5, i, facing = "clockwise")
}

circos.link("SCA", 0, "ENSMUSG0000040034", 0, h = 0.4)

circos.clear()

#par(op)

