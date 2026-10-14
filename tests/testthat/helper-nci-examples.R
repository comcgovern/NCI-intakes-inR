# Helpers for the NCI SAS example parity tests (test-sas-examples.R).
#
# NCI distributes the analytic dataset used by its MIXTRAN/DISTRIB v2.1
# example programs at
#   https://epi.grants.cancer.gov/diet/usualintakes/dataset.html
# The SAS listing output of those examples is frozen in
# fixtures/sas_examples_reference.csv.

nci_examples_enabled <- function() {
  identical(tolower(Sys.getenv("NCIUSUAL_SAS_PARITY")), "true")
}

# Path to indiv_ims.sas7bdat: $NCIUSUAL_NCI_DATA if set, otherwise download
# (about 14 MB) into a per-session cache.
nci_indiv_ims_path <- function() {
  p <- Sys.getenv("NCIUSUAL_NCI_DATA")
  if (nzchar(p) && file.exists(p)) return(p)
  dir <- file.path(tempdir(), "nciusual_nci_data")
  dir.create(dir, showWarnings = FALSE)
  f <- file.path(dir, "indiv_ims.sas7bdat")
  if (!file.exists(f)) {
    z <- file.path(dir, "indiv_ims.zip")
    utils::download.file(
      "https://epi.grants.cancer.gov/diet/usualintakes/indiv_ims.zip",
      z, mode = "wb", quiet = TRUE)
    utils::unzip(z, exdir = dir)
  }
  f
}

# Example 1: males 9+ (stratum 2), added sugars, amount model.
# Example 2: children 1-8 (stratum 1), total fruit, correlated model.
nci_example_data <- function(example) {
  d <- haven::read_sas(nci_indiv_ims_path())
  names(d) <- tolower(names(d))
  keep <- if (example == 1) {
    c("seqn", "day", "add_sug", "age4", "age5", "age6", "age7", "age8",
      "weekend", "rndw1", "agegrp", "race")
  } else {
    c("seqn", "day", "f_total", "sex1", "age2", "weekend", "rndw1",
      "agegrp", "race")
  }
  x <- as.data.frame(d[d$stratum == c(2, 1)[example], keep])
  x$race1 <- as.numeric(x$race == 1)
  x$race3 <- as.numeric(x$race == 3)
  as.data.frame(lapply(x, as.numeric))
}

nci_example_fit <- function(example) {
  x <- nci_example_data(example)
  if (example == 1) {
    mixtran(x, "add_sug", "seqn", "day", model_type = "amount",
            covariates = c("age4", "age5", "age6", "age7", "age8", "race1", "race3"),
            weekend_var = "weekend", weight_var = "rndw1",
            subgroup_var = "agegrp", verbose = FALSE)
  } else {
    mixtran(x, "f_total", "seqn", "day", model_type = "corr",
            covariates = c("sex1", "age2", "race1", "race3"),
            weekend_var = "weekend", weight_var = "rndw1",
            subgroup_var = "agegrp", verbose = FALSE)
  }
}

nci_example_cutpoints <- list(c(10, 12, 14, 17, 22, 25, 33, 38, 42, 47),
                              c(0.7, 1.1, 1.6))

sas_reference <- function() {
  read.csv(test_path("fixtures", "sas_examples_reference.csv"),
           stringsAsFactors = FALSE)
}
