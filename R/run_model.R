run_model <- function() {
  inputs <- readRDS(system.file("direct_incidence_model_input.rds", package = "cpptest"))
  tst_vector <- c(1,2,3)
  modify_vector(tst_vector)
  tst <- list(one = c(1,2,3),
              two = c(4,5,6))
  modify_list(tst)
}
