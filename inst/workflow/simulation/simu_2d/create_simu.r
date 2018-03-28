### SHELL> Rscript create_simu.r

phantom <- c("shepp0fMRI", "shepp1fMRI", "shepp2fMRI")
overlap <- c("0.01", "0.1", "0.25", "0.50", "0.75", "0.95")
min.1st.prop <- c("0.975")
max.PV <- c("0.05")
model <- c("I")
seed <- 1234

# system <- function(x) cat(x, "\n", sep = "")
mkdir.prefix <- TRUE
cp.prefix <- TRUE
mkdir.output <- TRUE

prefix <- "simu"
system(paste("mkdir ", prefix, sep = ""))
for(i.phantom in 1:length(phantom)){
  for(i.overlap in 1:length(overlap)){
    for(i.min.1st.prop in 1:length(min.1st.prop)){
      for(i.max.PV in 1:length(max.PV)){
        for(i.model in 1:length(model)){
          if(mkdir.prefix){
            # mkdir prefix
            cmd <- paste("cp -R script simu/",
                         phantom[i.phantom], "_",
                         overlap[i.overlap], "_",
                         min.1st.prop[i.min.1st.prop], "_",
                         max.PV[i.max.PV], "_",
                         model[i.model],
                         sep = "")
            system(cmd)
          }

          if(cp.prefix){
            # mkdir prefix
            cmd <- paste("cp -f script/* simu/",
                         phantom[i.phantom], "_",
                         overlap[i.overlap], "_",
                         min.1st.prop[i.min.1st.prop], "_",
                         max.PV[i.max.PV], "_",
                         model[i.model], "/",
                         sep = "")
            system(cmd)
          }

          if(mkdir.output){
            # mkdir output
            cmd <- paste("mkdir simu/",
                         phantom[i.phantom], "_",
                         overlap[i.overlap], "_",
                         min.1st.prop[i.min.1st.prop], "_",
                         max.PV[i.max.PV], "_",
                         model[i.model], "/output",
                         sep = "")
            system(cmd)

            # mkdir plot
            cmd <- paste("mkdir simu/",
                         phantom[i.phantom], "_",
                         overlap[i.overlap], "_",
                         min.1st.prop[i.min.1st.prop], "_",
                         max.PV[i.max.PV], "_",
                         model[i.model], "/plot",
                         sep = "")
            # system(cmd)
          }

          # phantom
          cmd <- paste("sed -i ",
                       "'s#shepp1fMRI#", phantom[i.phantom], "#'",
                       " simu/",
                       phantom[i.phantom], "_",
                       overlap[i.overlap], "_",
                       min.1st.prop[i.min.1st.prop], "_",
                       max.PV[i.max.PV], "_",
                       model[i.model], "/",
                       "00-set_condition.r",
                       sep = "")
          system(cmd)

          # overlap
          cmd <- paste("sed -i ",
                       "'s#overlap = 0.01#overlap = ", overlap[i.overlap], "#'",
                       " simu/",
                       phantom[i.phantom], "_",
                       overlap[i.overlap], "_",
                       min.1st.prop[i.min.1st.prop], "_",
                       max.PV[i.max.PV], "_",
                       model[i.model], "/",
                       "00-set_condition.r",
                       sep = "")
          system(cmd)

          # min.1st.prop
          cmd <- paste("sed -i ",
                       "'s#.FC.CT$INIT$min.1st.prop <- 0.95#.FC.CT$INIT$min.1st.prop <- ", min.1st.prop[i.min.1st.prop], "#'",
                       " simu/",
                       phantom[i.phantom], "_",
                       overlap[i.overlap], "_",
                       min.1st.prop[i.min.1st.prop], "_",
                       max.PV[i.max.PV], "_",
                       model[i.model], "/",
                       "00-set_condition.r",
                       sep = "")
          system(cmd)

          # set.seed
          cmd <- paste("sed -i ",
                       "'s#seed <- 1234#seed <- ", seed, "#'",
                       " simu/",
                       phantom[i.phantom], "_",
                       overlap[i.overlap], "_",
                       min.1st.prop[i.min.1st.prop], "_",
                       max.PV[i.max.PV], "_",
                       model[i.model], "/",
                       "00-set_condition.r",
                       sep = "")
          system(cmd)

          # model
          if(model[i.model] %in% c("I", "noX")){
            cmd <- paste("sed -i ",
                         "'s#.FC.CT$model.X <- \"V\"#.FC.CT$model.X <- \"I\"#'",
                         " simu/",
                         phantom[i.phantom], "_",
                         overlap[i.overlap], "_",
                         min.1st.prop[i.min.1st.prop], "_",
                         max.PV[i.max.PV], "_",
                         model[i.model], "/",
                         "00-set_condition.r",
                         sep = "")
            system(cmd)
          }

          # noX
          if(model[i.model] == "noX"){
            cmd <- paste("sed -i ",
                         "'s#.FC.CT$ignore.X <- FALSE#.FC.CT$ignore.X <- TRUE#'",
                         " simu/",
                         phantom[i.phantom], "_",
                         overlap[i.overlap], "_",
                         min.1st.prop[i.min.1st.prop], "_",
                         max.PV[i.max.PV], "_",
                         model[i.model], "/",
                         "00-set_condition.r",
                         sep = "")
            system(cmd)

            cmd <- paste("sed -i ",
                         "'s#.FC.CT$INIT$class.method <- \"prob.extend\"#.FC.CT$INIT$class.method <- \"extend\"#'",
                         " simu/",
                         phantom[i.phantom], "_",
                         overlap[i.overlap], "_",
                         min.1st.prop[i.min.1st.prop], "_",
                         max.PV[i.max.PV], "_",
                         model[i.model], "/",
                         "00-set_condition.r",
                         sep = "")
            system(cmd)
          }

          # run_1.sh
          cmd <- paste("sed -i ",
                       "'s#simu_2d#",
                       "s", i.phantom, i.overlap, i.min.1st.prop, i.model,
                       "#'",
                       " simu/",
                       phantom[i.phantom], "_",
                       overlap[i.overlap], "_",
                       min.1st.prop[i.min.1st.prop], "_",
                       max.PV[i.max.PV], "_",
                       model[i.model], "/",
                       "run_1.sh",
                       sep = "")
          system(cmd)

          # run_2.sh
          cmd <- paste("sed -i ",
                       "'s#simu_2d#",
                       "s", i.phantom, i.overlap, i.min.1st.prop, i.model,
                       "#'",
                       " simu/",
                       phantom[i.phantom], "_",
                       overlap[i.overlap], "_",
                       min.1st.prop[i.min.1st.prop], "_",
                       max.PV[i.max.PV], "_",
                       model[i.model], "/",
                       "run_2.sh",
                       sep = "")
          system(cmd)

          # create script run_1.sh
          cmd <- paste("cd simu/",
                       phantom[i.phantom], "_",
                       overlap[i.overlap], "_",
                       min.1st.prop[i.min.1st.prop], "_",
                       max.PV[i.max.PV], "_",
                       model[i.model], "/;",
                       " source run_1.sh;", 
                       " cd ../../",
                       sep = "")
          cat(cmd, "\n", sep = "", file = "submit_run_1.sh", append = TRUE)

          # create script run_2.sh
          cmd <- paste("cd simu/",
                       phantom[i.phantom], "_",
                       overlap[i.overlap], "_",
                       min.1st.prop[i.min.1st.prop], "_",
                       max.PV[i.max.PV], "_",
                       model[i.model], "/;",
                       " source run_2.sh;", 
                       " cd ../../",
                       sep = "")
          cat(cmd, "\n", sep = "", file = "submit_run_2.sh", append = TRUE)
        }
      }
    }
  }
}

# script_summary
cmd <- paste("sed -i ",
             "'s#set.seed(1234)#set.seed(", seed, ")#'",
             " script_summary/*",
             sep = "")
system(cmd)
