library(MixfMRI, quietly = TRUE)
set.seed(1234)
da <- gendataset(phantom = shepp1fMRI, overlap = 0.01)$pval
set.seed(1234)
da2 <- gendataset(phantom = shepp2fMRI, overlap = 0.01)$pval

xlim <- ylim <- c(0.19, 0.82)
# pdf(file = "maitra_2d.pdf", width = 4, height = 4)
# par(mfrow = c(2, 2), mar = c(0, 0, 1, 0))
par(mfrow = c(2, 2), mar = c(1, 1, 3, 1))
image(shepp1fMRI, xlim = xlim, ylim = ylim, axes = FALSE,
      col = gray(0:15 / 15), main = "shepp1fMRI")
image(da, xlim = xlim, ylim = ylim, axes = FALSE,
      col = rev(my.YlOrRd(10)), main = "overlap 0.01")
image(shepp2fMRI, xlim = xlim, ylim = ylim, axes = FALSE,
      col = gray(0:15 / 15), main = "shepp2fMRI")
image(da2, xlim = xlim, ylim = ylim, axes = FALSE,
      col = rev(my.YlOrRd(10)), main = "overlap 0.01")
# dev.off()
