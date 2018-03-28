library(MixfMRI, quietly = TRUE)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")
myPalette <- cbPalette[c(6, 7, 2, 8, 3, 5, 4, 1)]

par(mar = rep(0.1,4), mfrow = c(2, 2))

image(-shepp0fMRI, axes = F, col = c(myPalette[2:1],
      rgb(red = 230 / 255, green = 159 / 255, blue = 0, alpha = 0.1)))
contour(sheppAnat, drawlabels = F, lwd = 0.1,
        nlevels = length(unique(as.vector(sheppAnat))), add = T)

image(-shepp1fMRI, axes = F, col = c(myPalette[2:1],
      rgb(red = 230 / 255, green = 159 / 255, blue = 0, alpha = 0.1)))
contour(sheppAnat, drawlabels = F, lwd = 0.1,
        nlevels = length(unique(as.vector(sheppAnat))), add = T)

image(-shepp2fMRI, axes = F, col = c(myPalette[2:1],
      rgb(red = 230 / 255, green = 159 / 255, blue = 0, alpha = 0.1)))
contour(sheppAnat, drawlabels = F, lwd = 0.1,
        nlevels = length(unique(as.vector(sheppAnat))), add = T)

image(-shepp2fMRI[66:191, 49:208], axes = F, col = c(myPalette[2:1],
      rgb(red = 230 / 255, green = 159 / 255, blue = 0, alpha = 0.1)))
contour(sheppAnat[66:191, 49:208], drawlabels = F, lwd = 0.15,
        nlevels = length(unique(as.vector(sheppAnat))), add = T)
