case <- "ds.pval"

.FC.CT$algorithm <- "apecma"
# .FC.CT$model.X <- "V"
.FC.CT$ignore.X <- TRUE
.FC.CT$CONTROL$debug <- 1
.FC.CT$CONTROL$RndEM.iter <- 50
.FC.CT$INIT$min.1st.prop <- 0.99
.FC.CT$INIT$max.PV <- 0.05
.FC.CT$INIT$class.method <- "simple"
K.min <- 2
K.max <- 10
max.try.fclust <- 1
seed <- 1234

if(!file.exists("output")){
  dir.create("output")
}
if(!file.exists("plot")){
  dir.create("plot")
}
