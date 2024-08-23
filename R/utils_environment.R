#' Check if there is a loaded \code{conda} environment.
#'
#' @description Function to check if a \code{conda} environment is already loaded.
#'
#' @return Boolean, \code{TRUE} if an environment is loaded, \code{FALSE} otherwise.
#'
#' @importFrom reticulate py_discover_config
#' @export have_loaded_env

have_loaded_env = function() {
  return(grepl("envs", reticulate::py_discover_config()$python))
}


#' Check if \code{Python} packages are installed in the environment.
#'
#' @description Function to check if one or more \code{Python} packages are
#' present in a \code{conda} environment.
#'
#' @param envname the name of the environment to check. If empty, the function will check the
#' currently loaded environment.
#' @param py_pkgs a list or vector of \code{Python} packages.
#'
#' @return a list of Boolean. For each input package, \code{TRUE} if the package is
#' installed, \code{FALSE} otherwise.
#'
#' @importFrom reticulate py_list_packages conda_list py_discover_config
#' @export have_python_deps

have_python_deps = function(envname="", py_pkgs=c("pybascule")) {
  if (envname == "")
    envname = which_conda_env()

  tryCatch(
    expr = {
      ll = py_pkgs %in% reticulate::py_list_packages(envname)$package
      names(ll) = py_pkgs
      return(ll)
    },
    error = function(e) FALSE)
}


#' Retrieve the name of the currently loaded environment.
#'
#' @description Function to retrieve the name of the currently loaded \code{conda} environment.
#'
#' @return a \code{character} corresponding to the name of the loaded environment.
#'
#' @importFrom reticulate conda_list py_discover_config
#' @export which_conda_env

which_conda_env = function() {
  if (have_loaded_env())
    return(
      sapply(reticulate::conda_list()$name, grepl, reticulate::py_discover_config()$python) %>%
        which() %>%
        names()
    )

  return(cat("No loaded environments!"))
}


#' Load the input \code{conda} environment.
#'
#' @description Function to load the input \code{conda} environment. The function
#' will raise an error if a \code{Python} version has already been attached to the
#' \code{reticulate} package.
#' In that case, it will be necessary to restart the \code{R} session and to load the
#' desired environment **before** calling \code{lineaGT} function interfacing with
#' \code{Python} - \code{filter_dataset() and fit()}.
#'
#' @importFrom reticulate use_condaenv
#' @export load_conda_env

load_conda_env = function(envname="bascule-env") {
  Sys.unsetenv("RETICULATE_PYTHON")
  tryCatch(
    expr = {
      cli::cli_process_start(paste0("Loading the `", envname, "` environment."))
      reticulate::use_condaenv(envname, required=TRUE)
      cli::cli_process_done()
      },
    error = function(e) {
      cli::cli_alert_warning(paste0("The `", which_conda_env(), "` environment is already loaded!)"))
      cli::cli_alert_warning("To change the loaded environment, you need to restart the R session!")
    }
  )
}


have_conda_env = function(envname="bascule-env"){
  tryCatch(expr = envname %in% reticulate::conda_list()$name,
           error = function(e) FALSE )
}


# find out whether the user has conda installed and visible
have_conda = function() {
  conda_bin = tryCatch(reticulate::conda_binary("auto"),
                       error = function(e) NULL)
  !is.null(conda_bin)
}


create_conda_env = function(envname="bascule-env") {
  cli::cli_alert_warning(paste0("A new conda environment named `", envname, "` is being installed."))
  reticulate::conda_create(envname=envname, python_version="3.9")
}


using_conda_env = function(envname="bascule-env") {
  config = reticulate::py_discover_config()
  grepl(envname, config$python)
}


install_miniconda_bascule = function() {
  reticulate::install_miniconda()
}


install_python_deps = function(envname="bascule-env", pip=TRUE) {
  if (pip) {
    cli::cli_process_start("Installing the `pybascule` package from PyPI.")
    pkg = "pybascule"
  } else {
    cli::cli_process_start("Installing the `pybascule` package from GitHub. Insert the branch to use: ")
    branch = readline()
    pkg = paste0("git+https://github.com/caravagnalab/pybascule@", branch)
  }

  system(paste0(reticulate::py_config()$python, " -m pip install ", pkg, " --upgrade"))

  cli::cli_process_done()
}
