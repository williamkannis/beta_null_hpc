# Null model iteration function  -----------------------------------------------
null_iterations <- function(
    type,  # obs or null
    change = F,
    label = "",
    lcbd = NULL, # T or F
    beta_comps = NULL,  # all or vector of components matching BAT naming
    fun,  
    fun_args,
    algorithm = NULL,
    null.iter = NULL,
    null.cores = NULL
){
  
  if(!change){
    if(type == "obs"){
      obs <- .div_fun(
        label = label,
        lcbd = lcbd,
        beta_comps = beta_comps,
        fun = fun,
        fun_args = fun_args)
      return(obs)
    }
    
    if (type == "null"){
      
      # Set up null algorithm options
      args <- list(...)
      alg_input <- args[names(args) %in% c("comm","trait","tree")]
      alg_input <- c(alg_input,list(algorithm=algorithm))
      
      # estimate null iterations across processing cores
      null_iters <- parallel::mclapply(
        X = 1:null.iter,
        mc.cores = null.cores,
        FUN = function(i) {
          
          # Create null iteration
          null_input <- do.call(.null_algorithm,alg_input)
          null_args <- c(
            null_input,
            args[!names(args) %in% c("comm","trait","tree")]
          )
          
          # Estimate diversity
          .div_fun(
            label = label,
            lcbd = lcbd,
            beta_comps = beta_comps,
            fun = fun,
            fun_args = null_args
          )
        }
      )
    }
    return(null_iters)
  }
  # if(change){
  #   
  # }
}


# Helper functions  ------------------------------------------------------------

# estimate diversity using function of choice
.div_fun <- function(label = NULL,lcbd,beta_comps,fun,fun_args) {
  
  # Prepare custom labels if provided
  if(!is.null(label)) label <- paste0(label,"_")
  
  # Estimate diversity with specified method
  helper_fun <- get(fun)
  div <- do.call(helper_fun,fun_args)
  
  # Determine if alpha or beta
  metric = stringr::str_extract(fun, "(?<=_).*")
  
  if(metric == "alpha"){
    colnames(div) <- lapply(colnames(div),function(i) paste0(label,i))
    return(div)
  }
  
  # Select components and lcbd for beta diversity
  if (metric == "beta") {
    beta <- div[beta_comps]
    names(beta) <- lapply(names(beta),function(i) paste0(label,i))
    if (!lcbd) return(beta)
    
    # Create data frame with LCBD of all components
    lcbd <- purrr::map2(beta,names(beta),.lcbd_batch) %>% 
      purrr::reduce(left_join,by = join_by(COMID)) %>% 
      tibble::column_to_rownames("COMID")
    
    return(list(beta = beta,lcbd = lcbd))
  }
}

# Batch convert LCBD output to dataframe
.lcbd_batch <- function(D,name) {
  
  # check if euclidean
  sqrt <- ifelse(ade4::is.euclid(D) == T, F, T)
  
  # estimate LCBD
  lcbd <- adespatial::LCBD.comp(D,sqrt.D = sqrt)$LCBD
  
  # Create output dataframe
  out <- data.frame(COMID = labels(D))
  out[,name] <- lcbd
  row.names(out) <- out$COMID
  out
}



