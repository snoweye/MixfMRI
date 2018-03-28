case <- "simu.2d"

.FC.CT$algorithm <- "apecma"
.FC.CT$model.X <- "V"
.FC.CT$ignore.X <- FALSE
.FC.CT$CONTROL$debug <- 1
.FC.CT$CONTROL$RndEM.iter <- 50
.FC.CT$INIT$min.1st.prop <- 0.95
.FC.CT$INIT$class.method <- "prob.extend"
K.min <- 2
K.max <- 8
max.try.fclust <- 3
seed <- 1234
smooth <- FALSE

if(!file.exists("output")){
  dir.create("output")
}
if(!file.exists("plot")){
  dir.create("plot")
}

set.seed(seed)
simu.2d <- gendataset(phantom = shepp1fMRI, overlap = 0.01,
                      smooth = smooth)$pval
