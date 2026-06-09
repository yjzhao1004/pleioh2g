## File: R/globals.R
## Purpose: Declare global variables to silence R CMD check
##          “no visible binding for global variable …” warnings

utils::globalVariables(c(
  ".data","N", "SNP", "Z", "A1", "CHR", "h2_Z",
  "phenotype_path_update", "save_path", ".x"
))

