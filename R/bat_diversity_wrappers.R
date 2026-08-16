# Wrappers to easily call BAT functions

# FOR TESTING - TEST THAT VALUES EQUAL NON WRAPPER OUTPUTS
# AND THAT FUNCTIONS TAKE THE ARGUMENTS CORRECTLY

.kernel_alpha <- function(...){
  args <- list(...)
  kernel <- BAT::kernel.build(...)
  BAT::kernel.alpha(kernel)
}
.kernel_beta<- function(...){
  args <- list(...)
  kern_args <- args[names(args) != "func"]
  kernel <- do.call(BAT::kernel.build,kern_args)
  BAT::kernel.beta(kernel,func = args$func)
}

.hull_alpha <- function(...){
  hull <- BAT::hull.build(...)
  BAT::hull.alpha(hull)
}
.hull_beta <- function(...){
  args <- list(...)
  hull_args <- args[names(args) != "func"]
  hull <- do.call(BAT::hull.build,hull_args)
  BAT::hull.beta(hull,func = args$func)
}

.dendrogram_alpha <- function(...) BAT::alpha(...)
.dendrogram_beta <- function(...) BAT::beta(...)

.tax_alpha <-function(...) BAT::alpha(...)
.tax_beta <-function(...) BAT::beta(...)

# Reorder community and trait matrices to match
## REQUIRE MOST UPTODATE BAT PACKAGE, CHECK THAT BUG
# IS FIXED, THEN REOVE THIS HELPER
.trait_match <- function(comm,trait){
  sp <- colnames(comm)[order(colnames(comm))]
  com_match <- comm[,sp]
  trait_match <- trait[sp,]
  
  stopifnot(
    "Species in com do not match species in trait" = 
    all(colnames(com_match) == row.names(trait_match)))
  
  list(com_match,trait_match)
}