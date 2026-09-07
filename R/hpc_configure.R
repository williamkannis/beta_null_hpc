hpc_config <- function(
    host,
    hpc_dir,
    username = NULL
) {
  
  options(
    beta_hpc_host = host,
    beta_hpc_dir = hpc_dir,
    beta_hpc_username = username
  )
  
  invisible(NULL)
}