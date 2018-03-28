### Controls

.FC.CT <- list(
  algorithm = c("apecma", "em", "ecm"),
  optim.method = c("BFGS", "Nelder-Mead"),
  model.X = c("I", "V"),    # independent or unstructure.
  ignore.X = FALSE,         # if using location information.
  check.X.unit = TRUE,      # if checking X in [0, 1].
  CONTROL = list(
              max.iter = 1000,
              abs.err = 1e-4,
              rel.err = 1e-6,
              debug = 1,
              RndEM.iter = 10,
              exp.min = log(.Machine$double.xmin),
              exp.max = log(.Machine$double.xmax),
              sigma.ill = 1e-6,
              DS.max = 1e+4,
              DS.min = 1e-6
            ),
  INIT = list(
              min.1st.prop = 0.8,
              max.PV = 0.1,
              BETA.alpha.min = 0 + 1e-6,
              BETA.alpha.max = 1 - 1e-6,
              BETA.beta.min = 1 + 1e-6,
              BETA.beta.max = 1e+6,
              max.try.iter = 10,
              class.method = c("prob.extend", "prob.simple",
                               "qnorm.extend", "qnorm.simple",
                               "extend", "simple")
             ),
  LRT = list(
             H0.alpha = 1,
             H0.beta = 1,
             H0.mean = 0.05
            ),
  ### For parallel
  MPI.gbd = FALSE,
  common.gbd = TRUE
)
