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

  load_conda_env(envname)

  # check if the required python packages are installed in the environment
  check_python_deps(envname=envname, pip=TRUE)
}


check_conda = function(use_default=F) {
  if (!have_conda()) {
    cat("It is not prossibe to find a Anaconda or Miniconda version installed.\n")

    if (use_default) answ = "yes"
    else {
      cli::cli_alert_warning("The Miniconda installer will be downloaded and used to install Miniconda. Proceed? (yes/no)\n")
      answ = readline()
    }

    if (answ == "yes") {
      install_miniconda_lineagt()
      create_conda_env()
      load_conda_env()
    }
    else { cli::cli_alert_info("Miniconda will not be installed."); return() }
  }

}



check_conda_env = function(envname="basilica-env", use_default=F) {
  if (have_loaded_env()) {
    envname = sapply(reticulate::conda_list()$name, grepl, reticulate::py_discover_config()$python) %>%
      which() %>% names()
    cli::cli_alert_warning(paste0("The '", envname, "' environment is already loaded!"))
    return(envname)
  }

  if (!have_conda_env("basilica-env")) {
    cli::cli_alert_info("The environment 'basilica-env' is not present.\n")

    if (use_default) answ = "create" else {
      cli::cli_alert_warning("Do you want to load an existing environment, to create a new one named 'basilica-env' or to cancel? (load/create/cancel)\n")
      answ = readline()
    }

    if (answ == "create") {
      envname = "basilica-env"
      create_conda_env()
    } else if (answ == "load") {
      cli::cli_alert_info("Insert the environment name: ")
      envname = readline()
    } else {
      cli::cli_alert_info("No environment will be loaded nor created.")
      return()
    }

  } else {
    cli::cli_alert_info("The environment 'basilica-env' is already present and will be loaded!\n")
    envname = "basilica-env"
  }

  if (!have_loaded_env())
    load_conda_env(envname)

  return(envname)
}



check_python_deps = function(envname="basilica-env", pip=FALSE) {
  tryCatch(
    expr = install_python_deps(envname, pip=pip),
    error = function(e) cli::cli_alert_warning("Not able to install the Python package.")
  )
}


plot_catalogue = function(catalogue) {
  catalogue %>%
    tibble::rownames_to_column() %>% reshape2::melt() %>% dplyr::mutate(variable=stringr::str_replace_all(variable, "\\[|\\]", "_")) %>% tidyr::separate(variable, into=c("left","subs","right"), sep="_") %>% dplyr::mutate(context=paste0(left,"_",right), left=NULL, right=NULL) %>% ggplot() + geom_bar(aes(x=context, y=value, fill=rowname), stat="identity") + facet_grid(rowname~subs)
}

