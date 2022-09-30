#' Configure the reticulate environment
#'
#' @description Function to configure the Python dependencies on R.
#' If a Python environment is not available, the function will check if there is a version of
#' \code{conda} or \code{miniconda}, otherwise it will install \code{miniconda}, on which
#' install the Python package \code{pybasilica}.
#'
#' @param env_name name of the \code{conda} environment to use, if available.
#'
#' @importFrom reticulate import conda_create conda_install install_miniconda miniconda_path conda_binary
#' @export configure_environment

configure_environment = function(envname="basilica-env") {
  # install miniconda if needed
  check_conda()

  # check a conda environment is present or load one
  envname = check_conda_env(envname=envname)

  # check if the required python packages are installed in the environment
  check_python_deps(envname=envname)
}


check_conda = function() {
  if (!have_conda()) {
    cat("It is not possibe to find a Anaconda or Miniconda version installed.\n A Miniconda installation will be prompted.")
    install_miniconda_lineagt()
  }
}


check_conda_env = function(envname="basilica-env") {
  if (have_loaded_env()) {
    envname = sapply(reticulate::conda_list()$name, grepl, reticulate::py_discover_config()$python) %>%
      which() %>% names()
    return(envname)
  }

  if (!have_conda_env("basilica-env")) {
    cat("The environment 'basilica-env' is not present.\n")
    envname = "basilica-env"
    create_conda_env(envname=envname)
  }

  load_conda_env(envname)

  return(envname)
}


check_python_deps = function(envname="basilica-env") {
  try(install_python_deps(envname), silent=T)
}

