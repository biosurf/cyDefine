#' Inverted version of in
#'
#' @noRd
#' @examples
#' 1 %!in% 1:10
`%!in%` <- Negate(`%in%`)

#' weighted Mahalanobis distance
#'
#' @noRd
wmahalanobis <- function(x, center, cov, weight) {
  if (is.vector(x)) {
    x <- matrix(x, ncol = length(x))
  } else {
    x <- as.matrix(x)
  }
  x <- sweep(x, 2, center)
  cov <- weight %*% solve(cov)
  retval <- diag(x %*% cov %*% t(x))
  retval
}

#' Check for ID column
#'
#' @noRd
check_id <- function(df) {
  if (!"id" %in% colnames(df)) df$id <- seq_len(nrow(df))
  return(df)
}

#' Check package
#'
#' @param package Name of package
#' @param repo Repository of package, choose from "CRAN", "Bioc", and "github"
#' @param git_repo If repo is github, name of the GitHub repository
#' @param required Whether the package is required - will break the program if
#' required but not installed.
#'
#' @return Boolean - TRUE if package is present, FALSE if not.
#'
check_package <- function(package,
                          repo = "CRAN",
                          git_repo = "",
                          required = TRUE) {
  # if package not installed
  if (!requireNamespace(package, quietly = TRUE)) {
    if (required) {
      if (repo == "CRAN") {
        install_function <- "install.packages('"
      } else if (repo == "github") {
        install_function <- "remotes::install_github('"
        package <- git_repo
      } else if (repo == "Bioc") {
        install_function <- "BiocManager::install('"
      }

      stop(
        "Package ", package, " is not installed.\n",
        "Please run: ", install_function, package, "')"
      )
    }
    return(FALSE)
  }
  return(TRUE)
}



#' Check column names
#'
#' @param avail_colnames Available column names to check if include
#' 'req_colnames'
#' @param req_colnames Required column names to check if are among
#' 'avail_colnames'
#'
#' @return NULL
#'
check_colnames <- function(avail_colnames, req_colnames) {
  if (length(setdiff(req_colnames, avail_colnames)) > 0) {
    stop(
      "Required column(s) ",
      paste0(setdiff(req_colnames, avail_colnames), collapse = ", "),
      " were not found."
    )
  }
}
