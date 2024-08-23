.onLoad <- function(libname, pkgname) {

  pk = 'bascule'
  www = "https://caravagnalab.github.io/bascule/"

  cli::cli_alert_success( 'Loading {.field {pk}}. Support : {.url { www}}' )

  # configure_environment()
}
