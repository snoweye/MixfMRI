### SHELL> cd ../; Rscript simu_2d/create_rep.r

seed <- 1234

for(i.rep in 1:25){
  # mkdir rep
  cmd <- paste("cp -R simu_2d/ ",
               "simu_2d_", sprintf("%02d", i.rep),
               sep = "")
  system(cmd)

  # change seed
  cmd <- paste("cd simu_2d_", sprintf("%02d", i.rep), "/; ",
               "sed -i ",
               "'s#^seed <- ", seed, "#seed <- ", (seed + i.rep), "#' ",
               "create_simu.r",
               sep = "")
  system(cmd)

  # create simu
  cmd <- paste("cd simu_2d_", sprintf("%02d", i.rep), "/; ",
               "Rscript create_simu.r &",
               sep = "")
  system(cmd)

  # create script
  cmd <- paste("cd simu_2d_", sprintf("%02d", i.rep), "/; ",
               "source submit_run_1.sh; cd ../",
               sep = "")
  cat(cmd, "\n", sep = "", file = "q_run_1.sh", append = TRUE)

  cmd <- paste("cd simu_2d_", sprintf("%02d", i.rep), "/; ",
               "source submit_run_2.sh; cd ../",
               sep = "")
  cat(cmd, "\n", sep = "", file = "q_run_2.sh", append = TRUE)

  cmd <- paste("cd simu_2d_", sprintf("%02d", i.rep), "/; ",
               "source run_summary.sh; cd ../",
               sep = "")
  cat(cmd, "\n", sep = "", file = "q_run_3.sh", append = TRUE)
}
