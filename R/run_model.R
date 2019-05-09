run_model <- function() {
  inputs <- readRDS(system.file("direct_incidence_model_input.rds", package = "cpptest"))
  tst_vector <- c(1,2,3)
  modify_vector(tst_vector)
  tst <- list(one = c(1,2,3),
              two = c(4,5,6))
  modify_list(tst)
}


test_models <- function() {
  array1 <- array(1:24, c(3,4,2))
  array2 <- array(31:54, c(3,4,2))

  arrays <- list(first = array1,
                 second = array2)

  dims <- lapply(arrays, dim)

  push_list_arrays(arrays)

}

test_run_model <- function(timesteps) {
  inputs <- readRDS(system.file("direct_incidence_model_input.rds", package = "cpptest"))
  pop_dims <- list("ages" = 15:80,
                   "sexes" = c("Male", "Female"),
                   "hiv_status" = c("HIVN", "HIVP"),
                   "proj_years" = 1970:2025)
  inputs$eppmodInt <- match(inputs$eppmod, c("rtrend", "directincid"), nomatch = 0)
  inputs$incidmodInt <- match(inputs$incidmod, c("eppspectrum")) - 1L
  model_run <- runModel(inputs$basepop, inputs$ss$h.ag.span, inputs$tARTstart,
                        inputs$entrantprev, inputs$verttrans_lag,
                        inputs$paedsurv_lag, inputs$popadjust, inputs$entrantpop,
                        inputs$birthslag, inputs$cumsurv, inputs$cumnetmigr,
                        inputs$netmig_hivprob, inputs$paedsurv_cd4dist,
                        inputs$entrantartcov, inputs$paedsurv_artcd4dist,
                        inputs$Sx, inputs$netmigr, inputs$asfr, inputs$srb,
                        inputs$ss$hiv_steps_per_year, inputs$cd4_prog,
                        inputs$cd4_initdist, inputs$cd4_mort, inputs$incrr_age,
                        inputs$relinfectART, inputs$iota, inputs$incrr_sex,
                        inputs$incidmodInt, inputs$eppmodInt, inputs$scale_cd4_mort,
                        inputs$proj.steps, inputs$tsEpidemicStart,
                        inputs$artcd4elig_idx, inputs$specpop_percelig,
                        inputs$pw_artelig, inputs$who34percelig, inputs$rvec,
                        inputs$rtrend$beta, inputs$rtrend$tStabilize,
                        inputs$rtrend$r0, timesteps)
  array(
    model_run, lengths(pop_dims), pop_dims
  )
}
