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

  test_3d_arrays()


  array1 <- 1:10
  array2 <- 11:20

  arrays <- list(first = array1,
                 second = array2)
  test_1d_arrays(arrays)
}
