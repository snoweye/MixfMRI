library(MixfMRI, quietly = TRUE)
library(RColorBrewer, quietly = TRUE)

##
## 3 phantoms with  proportion of activated voxels:
## shepp0fMRI:  0.917%
## shepp1fMRI:  2.249%
## shepp2fMRI:  3.968%
##
## plotting the phantom
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")
myPalette <- cbPalette[c(6, 7, 2, 8, 3, 5, 4, 1)]
myCol <- c(myPalette[2:1],
           rgb(red = 230/255, green = 159/255, blue = 0, alpha = 0.1))

xlim <- ylim <- c(0.19, 0.82)
# pdf(file = "maitra_phantom.pdf", width = 6, height = 2)
# par(mfrow = c(1, 3), mar = c(0, 0, 1, 0))
par(mfrow = c(2, 2), mar = c(1, 1, 3, 1))
image(-shepp2fMRI, axes = FALSE, xlim = xlim, ylim = ylim, col = myCol,
      main = "shepp2fMRI: 3.968%")
contour(matrix(sheppAnat, ncol = 256), drawlabels = FALSE, lwd = 0.15,
        nlevels = length(unique(as.vector(sheppAnat))), add = TRUE)

image(-shepp1fMRI, axes = FALSE, xlim = xlim, ylim = ylim, col = myCol,
      main = "shepp1fMRI: 2.249%")
contour(matrix(sheppAnat, ncol = 256), drawlabels = FALSE,lwd = 0.15,
        nlevels = length(unique(as.vector(sheppAnat))), add = TRUE)

image(-shepp0fMRI, axes = FALSE, xlim = xlim, ylim = ylim, col = myCol,
      main = "shepp0fMRI: 0.917%")
contour(matrix(sheppAnat, ncol = 256), drawlabels = FALSE, lwd = 0.15,
        nlevels = length(unique(as.vector(sheppAnat))), add = TRUE)
# dev.off()
