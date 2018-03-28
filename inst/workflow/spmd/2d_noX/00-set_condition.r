case <- "pval.olddata.complex"
# case <- "pval.olddata.mag"

.FC.CT$algorithm <- "apecma"
# .FC.CT$model.X <- "I"
.FC.CT$ignore.X <- TRUE
.FC.CT$CONTROL$debug <- 1
.FC.CT$CONTROL$RndEM.iter <- 100
.FC.CT$INIT$min.1st.prop <- 0.8
.FC.CT$INIT$class.method <- "simple"
K.min <- 2
K.max <- 10
max.try.fclust <- 3
seed <- 1234

if(!file.exists("output")){
  dir.create("output")
}
if(!file.exists("plot")){
  dir.create("plot")
}
