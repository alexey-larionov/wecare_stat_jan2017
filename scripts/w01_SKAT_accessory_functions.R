# SKAT_accessory_functions.R
# Author: AL
# Started: 19Jul2016
# Last update: 24Jul2016

# Convinience functions to run multiple SKAT tests

library(tools) # for assertWarning()

#-------------------------------------------------------------------#
#                                                                   #
#           function_to_run_tests_for_multiple_variants             #
#                                                                   #
#-------------------------------------------------------------------#
#set.name <- gene

run_tests_for_multiple_variants <- function(
  set.name, phenotypes, covariates, genotypes){
  
  # Calculate logistic weights
  if( sum(genotypes, na.rm=TRUE) > 0 ){
    lw <- Get_Logistic_Weights(genotypes, 0.07, 150)
  }else{
    # Get_Logistic_Weights generates error in case of uniform genotypes
    lw <- rep(1,ncol(genotypes))
  }
  
  # Calculate null-model
  # Keep adjustment = TRUE because it will be used bu SKAT_RareCommon
  skat.null <- tryCatch(SKAT_Null_Model(phenotypes ~ covariates, out_type="D"), error=function(e) e)
  
  # Run Burden tests
  burden.test.df_wt <- tryCatch(SKATBinary(genotypes, skat.null, method="Burden", weights.beta=c(1,25)), error=function(e) e)
  burden.test.mb_wt <- tryCatch(SKATBinary(genotypes, skat.null, method="Burden", weights.beta=c(0.5,0.5)), error=function(e) e)
  burden.test.lw_wt <- tryCatch(SKATBinary(genotypes, skat.null, method="Burden", weights=lw), error=function(e) e)
  
  # Run SKAT tests
  skat.test.df_wt <- tryCatch(SKATBinary(genotypes, skat.null, method="SKAT", weights.beta=c(1,25)), error=function(e) e)
  skat.test.mb_wt <- tryCatch(SKATBinary(genotypes, skat.null, method="SKAT", weights.beta=c(0.5,0.5)), error=function(e) e)
  skat.test.lw_wt <- tryCatch(SKATBinary(genotypes, skat.null, method="SKAT", weights=lw), error=function(e) e)
  
  # Run SKATO tests with beta weights
  skato.test.df_wt <- tryCatch(SKATBinary(genotypes, skat.null, method="SKATO", weights.beta=c(1,25)), error=function(e) e)
  skato.test.mb_wt <- tryCatch(SKATBinary(genotypes, skat.null, method="SKATO", weights.beta=c(0.5,0.5)), error=function(e) e)
  
  # For ~0.5% of genes SKATO crashed with memory misallocation when weights were explicitly set to lw.
  # Because this crached R completely it was triky to handle with tryCatch or other R error-handling procedures. 
  
  # I noted that some genes causing SKATO crashes had similar features: two variants, with one het for each variant,
  # lw=c(1,1) etc. In particular, it looked like crashes happened when the previous SKATO calls with beta-weighting
  # produced warning: "Rank of the genotype matrix is one! SKAT is used instead of SKAT-O!"
  
  # So I desided NOT to run this skato(lw) in such cases, using warnings from SKATO calls with beta-weights. 
  # To get warnings I used tools::assertWarning. Unfortunately, it also generated errors in ~0.1% genes.
  # However, these were not complete crashes, so I managed it with tryCatch. 
  
  a <- tryCatch(assertWarning(SKATBinary(genotypes, skat.null, method="SKATO", weights.beta=c(0.5,0.5))),
                error=function(e) list("Rank of the genotype matrix is one! SKAT is used instead of SKAT-O!",e))
  
  b <- grepl("Rank of the genotype matrix is one! SKAT is used instead of SKAT-O!", unlist(a))
  if(any(b)){
    skato.test.lw_wt <- list(NA)
  }else{
    skato.test.lw_wt <- tryCatch(SKATBinary(genotypes, skat.null, method="SKATO", weights=lw), error=function(e) e)
  }

  # Run SKAT_CommonRare tests
  cr.test.skat <- tryCatch(SKAT_CommonRare(genotypes, skat.null, r.corr.rare=0, r.corr.common=0), error=function(e) e)
  cr.test.burden <- tryCatch(SKAT_CommonRare(genotypes, skat.null, r.corr.rare=1, r.corr.common=1), error=function(e) e)
  
  # Compile results list
  
  result <- list(
    set.name, 
    phenotypes, 
    covariates, 
    genotypes, 
    skat.null, 
    burden.test.df_wt,
    burden.test.mb_wt,
    burden.test.lw_wt,
    skat.test.df_wt,
    skat.test.mb_wt,
    skat.test.lw_wt,
    skato.test.df_wt,
    skato.test.mb_wt,
    skato.test.lw_wt,
    cr.test.skat,
    cr.test.burden)
  
  names(result) <- c(
    "set.name", 
    "phenotypes", 
    "covariates", 
    "genotypes", 
    "skat.null",
    "burden.skatw",
    "burden.mbw",
    "burden.lw",
    "skat.skatw",
    "skat.mbw",
    "skat.lw",
    "skato.skatw",
    "skato.mbw",
    "skato.lw",
    "cr.skat",
    "cr.burden")
  
  # Return result
  return(result)
  
}

#-------------------------------------------------------------------#
#                                                                   #
#               function_to_prepare_tables_for_results              #
#                                                                   #
#-------------------------------------------------------------------#

# Makes empty tables with headers only. 
# Deletes previously existing files, if any. 

prepare_tables_for_results <- function(
  prefix, report_type = "summary_only", path=""){
  
  # Prepare headers
  summary_table_header <- c("set", "n_cases_init", "n_vars_init", "AC_UBC_init", "AC_CBC_init", 
                            "missed_genotypes_fraction_init", "n_cases_retained", "n_vars_retained", "n_vars_retained_cr", 
                            "burden_skatw", "burden_mbw", "burden_lw", "skat_skatw", "skat_mbw", "skat_lw", 
                            "skato_skatw", "skato_mbw", "skato_lw", "cr_skat",  "cr_burden", "min_p")
  
  details_table_header <- c("set", "n_cases_init", "n_vars_init", 
                            "missed_genotypes_fraction_init", "n_cases_retained", 
                            "p_val", "n_vars_total", "n_vars_tested", 
                            "ma_count", "n_cases_with_ma", "method_for_p_assessment",
                            "rare_vars_tested", "common_vars_tested", "cutoff")
  
  # Write header to summary table
  write(paste(summary_table_header, sep="", collapse="\t"),
        file=paste(path, prefix, "_summary.txt", sep="", collapse=""))
  
  # Write headers to other tables, if required
  if(report_type == "full"){
    
    write(paste(details_table_header, sep="", collapse="\t"),
          file=paste(path, prefix, "_burden_skatw.txt", sep="", collapse=""))
    
    write(paste(details_table_header, sep="", collapse="\t"),
          file=paste(path, prefix, "_burden_mbw.txt", sep="", collapse=""))
    
    write(paste(details_table_header, sep="", collapse="\t"),
          file=paste(path, prefix, "_burden_lw.txt", sep="", collapse=""))
    
    write(paste(details_table_header, sep="", collapse="\t"),
          file=paste(path, prefix, "_skat_skatw.txt", sep="", collapse=""))
    
    write(paste(details_table_header, sep="", collapse="\t"),
          file=paste(path, prefix, "_skat_mbw.txt", sep="", collapse=""))
    
    write(paste(details_table_header, sep="", collapse="\t"),
          file=paste(path, prefix, "_skat_lw.txt", sep="", collapse=""))
    
    write(paste(details_table_header, sep="", collapse="\t"),
          file=paste(path, prefix, "_skato_skatw.txt", sep="", collapse=""))
    
    write(paste(details_table_header, sep="", collapse="\t"),
          file=paste(path, prefix, "_skato_mbw.txt", sep="", collapse=""))
    
    write(paste(details_table_header, sep="", collapse="\t"),
          file=paste(path, prefix, "_skato_lw.txt", sep="", collapse=""))
    
    write(paste(details_table_header, sep="", collapse="\t"),
          file=paste(path, prefix, "_cr_skat.txt", sep="", collapse=""))
    
    write(paste(details_table_header, sep="", collapse="\t"),
          file=paste(path, prefix, "_cr_burden.txt", sep="", collapse=""))
    
  }
}

#-------------------------------------------------------------------#
#                                                                   #
#       functions_to_check_existance_of_fields_in_skat_output       #
#                                                                   #
#-------------------------------------------------------------------#

# SKAT functions may choose different algorithms and output fields depending on data. 
# This function is used to avoid errors if some fields do not exist. 
# If the field does not exist, the function outputs NA instead of error.
# It works well with skat output objects in required contect. 
# It is not intended and it may not work with other objects or context. 
# In particular, it does not work with atomic input objects. 

nm <- function(x){if(is.numeric(x)){x}else{NA}}
ch <- function(x){if(is.character(x)){x}else{NA}}


#-------------------------------------------------------------------#
#                                                                   #
#        function_to_extract_p_values_from_skat_test_object         #
#                                                                   #
#-------------------------------------------------------------------#

get_p_values <- function(x){
  
  #summary_table_header <- c("set", "n_cases_init", "n_vars_init", 
  #  "missed_genotypes_fraction_init", "n_cases_retained", "n_vars_retained", "n_vars_retained_cr", 
  #  "burden_skatw", "burden_mbw", "burden_lw", "skat_skatw", "skat_mbw", "skat_lw", 
  #  "skato_skatw", "skato_mbw", "skato_lw", "cr_skat",  "cr_burden", "min_p")
  
  c(nm(x$burden.skatw$p.value), # burden_skatw
    nm(x$burden.mbw$p.value), # burden_mbw
    nm(x$burden.lw$p.value), # burden_lw
    nm(x$skat.skatw$p.value), # skat_skatw
    nm(x$skat.mbw$p.value), # skat_mbw
    nm(x$skat.lw$p.value), # skat_lw
    nm(x$skato.skatw$p.value), # skato_skatw
    nm(x$skato.mbw$p.value), # skato_mbw
    nm(x$skato.lw$p.value), # skato_lw
    nm(x$cr.skat$p.value), # cr_skat
    nm(x$cr.burden$p.value)) # cr_burden
  
}

#-------------------------------------------------------------------#
#                                                                   #
#        function_to_extract_detailes_from_skat_test_object         #
#                                                                   #
#-------------------------------------------------------------------#

get_details <- function(x){
  
  #details_table_header <- c("set", "n_cases_init", "n_vars_init", 
  #  "missed_genotypes_fraction_init", "n_cases_retained", 
  #  "p_val", "n_vars_total", "n_vars_tested", 
  #  "ma_count", "n_cases_with_ma", "method_for_p_assessment",
  #  "rare_vars_tested", "common_vars_tested", "cutoff")
  
  c(nm(x$p.value), # p_val
    nm(x$param$n.marker), # n_vars_total
    nm(x$param$n.marker.test), # n_vars_tested
    nm(x$MAC), # ma_count
    nm(x$m), # n_cases_with_ma
    ch(x$method.bin), # method_for_p_assessment
    nm(x$n.rare), # rare_vars_tested
    nm(x$n.common), # common_vars_tested
    nm(x$Cutoff)) # common_vars_tested
  
}

#-------------------------------------------------------------------#
#                                                                   #
#         function_to_print_results_for_multiple_variants           #
#                                                                   #
#-------------------------------------------------------------------#
#data <- skat.test

print_results_for_multiple_variants <- function(
  data, prefix, AC_UBC_init, AC_CBC_init, report_type = "summary_only", path=""){
  
  # Set name
  set <- data$set.name
  
  # Total num of variants and fraction of missed genotypes
  n_cases_init <- nrow(data$genotypes)
  n_vars_init <- ncol(data$genotypes)
  missed_genotypes_fraction_init <- sum(is.na(data$genotypes)) / (n_cases_init * n_vars_init)
  
  # Num of retained cases
  n_cases_retained <- length(data$skat.null$re1$id_include)
  
  # Num of retained variants
  n_vars_retained <- nm(data$burden.skatw$param$n.marker.test)

  if (is.na(nm(data$burden.mbw$param$n.marker.test)) | 
      is.na(nm(data$burden.lw$param$n.marker.test)) | 
      is.na(nm(data$skat.skatw$param$n.marker.test)) | 
      is.na(nm(data$skat.mbw$param$n.marker.test)) | 
      is.na(nm(data$skat.lw$param$n.marker.test)) | 
      is.na((data$skato.skatw$param$n.marker.test)) | 
      is.na(nm(data$skato.mbw$param$n.marker.test)) | 
      is.na(nm(data$skato.lw$param$n.marker.test))){
    NA -> n_vars_retained
    }else{
      if (n_vars_retained != nm(data$burden.mbw$param$n.marker.test) | 
        n_vars_retained != nm(data$burden.lw$param$n.marker.test) | 
        n_vars_retained != nm(data$skat.skatw$param$n.marker.test) | 
        n_vars_retained != nm(data$skat.mbw$param$n.marker.test) | 
        n_vars_retained != nm(data$skat.lw$param$n.marker.test) | 
        n_vars_retained != nm(data$skato.skatw$param$n.marker.test) | 
        n_vars_retained != nm(data$skato.mbw$param$n.marker.test) | 
        n_vars_retained != nm(data$skato.lw$param$n.marker.test)){
      NA -> n_vars_retained
      }
    }
  
  # Num of retained variants in cr tests
  n_vars_retained_cr <- nm(data$cr.skat$param$n.marker.test)
  
  if (is.na(nm(data$cr.skat$param$n.marker.test)) | 
      is.na(nm(data$cr.burden$param$n.marker.test))){
    NA -> n_vars_retained_cr
  }else{
    if (n_vars_retained_cr != nm(data$cr.burden$param$n.marker.test)){
      NA -> n_vars_retained_cr
    }
  }

  # p-values
  p_values <- get_p_values(data)
  min_p <- min(p_values, na.rm=TRUE)
  
  #summary_table_header <- c("set", "n_cases_init", "n_vars_init", "AC_UBC_init", "AC_CBC_init", 
  #  "missed_genotypes_fraction_init", "n_cases_retained", "n_vars_retained", "n_vars_retained_cr", 
  #  "burden_skatw", "burden_mbw", "burden_lw", "skat_skatw", "skat_mbw", "skat_lw", 
  #  "skato_skatw", "skato_mbw", "skato_lw", "cr_skat",  "cr_burden", "min_p")
  
  # Compile data for summary table 
  # NB: should be consistent with the header!!!
  summary <- c(set, n_cases_init, n_vars_init, AC_UBC_init, AC_CBC_init, 
               missed_genotypes_fraction_init, n_cases_retained, 
               n_vars_retained, n_vars_retained_cr, p_values, min_p)
  
  # Write data to summary table
  write(paste(summary, sep="", collapse="\t"),
        file=paste(path, prefix, "_summary.txt", sep="", collapse=""),
        append = TRUE)
  
  # Write data to other tables, if required
  if (report_type == "full") {
    
    # Compile data
    
    #details_table_header <- c("set", "n_cases_init", "n_vars_init", 
    #  "missed_genotypes_fraction_init", "n_cases_retained", 
    #  "p_val", "n_vars_total", "n_vars_tested", 
    #  "ma_count", "n_cases_with_ma", "method_for_p_assessment",
    #  "rare_vars_tested", "common_vars_tested", "cutoff")
    
    initial_fields <- c(set, n_cases_init, n_vars_init, 
                        missed_genotypes_fraction_init, 
                        n_cases_retained)
    
    burden_skatw <- get_details(data$burden.skatw)
    burden_skatw <- c(initial_fields, burden_skatw)
    
    burden_mbw <- get_details(data$burden.mbw)
    burden_mbw <- c(initial_fields, burden_mbw)
    
    burden_lw <- get_details(data$burden.lw)
    burden_lw <- c(initial_fields, burden_lw)
    
    skat_skatw <- get_details(data$skat.skatw)
    skat_skatw <- c(initial_fields, skat_skatw)
    
    skat_mbw <- get_details(data$skat.mbw)
    skat_mbw <- c(initial_fields, skat_mbw)
    
    skat_lw <- get_details(data$skat.lw)
    skat_lw <- c(initial_fields, skat_lw)
    
    skato_skatw <- get_details(data$skato.skatw)
    skato_skatw <- c(initial_fields, skato_skatw)
    
    skato_mbw <- get_details(data$skato.mbw)
    skato_mbw <- c(initial_fields, skato_mbw)
    
    skato_lw <- get_details(data$skato.lw)
    skato_lw <- c(initial_fields, skato_lw)
    
    cr_skat <- get_details(data$cr.skat)
    cr_skat <- c(initial_fields, cr_skat) 
    
    cr_burden <- get_details(data$cr.burden)
    cr_burden <- c(initial_fields, cr_burden)
    
    # Write data to files
    
    write(paste(burden_skatw, sep="", collapse="\t"),
          file=paste(path, prefix, "_burden_skatw.txt", sep="", collapse=""),
          append = TRUE)
    
    write(paste(burden_mbw, sep="", collapse="\t"),
          file=paste(path, prefix, "_burden_mbw.txt", sep="", collapse=""),
          append = TRUE)
    
    write(paste(burden_lw, sep="", collapse="\t"),
          file=paste(path, prefix, "_burden_lw.txt", sep="", collapse=""),
          append = TRUE)
    
    write(paste(skat_skatw, sep="", collapse="\t"),
          file=paste(path, prefix, "_skat_skatw.txt", sep="", collapse=""),
          append = TRUE)
    
    write(paste(skat_mbw, sep="", collapse="\t"),
          file=paste(path, prefix, "_skat_mbw.txt", sep="", collapse=""),
          append = TRUE)
    
    write(paste(skat_lw, sep="", collapse="\t"),
          file=paste(path, prefix, "_skat_lw.txt", sep="", collapse=""),
          append = TRUE)
    
    write(paste(skato_skatw, sep="", collapse="\t"),
          file=paste(path, prefix, "_skato_skatw.txt", sep="", collapse=""),
          append = TRUE)
    
    write(paste(skato_mbw, sep="", collapse="\t"),
          file=paste(path, prefix, "_skato_mbw.txt", sep="", collapse=""),
          append = TRUE)
    
    write(paste(skato_lw, sep="", collapse="\t"),
          file=paste(path, prefix, "_skato_lw.txt", sep="", collapse=""),
          append = TRUE)
    
    write(paste(cr_skat, sep="", collapse="\t"),
          file=paste(path, prefix, "_cr_skat.txt", sep="", collapse=""),
          append = TRUE)
    
    write(paste(cr_burden, sep="", collapse="\t"),
          file=paste(path, prefix, "_cr_burden.txt", sep="", collapse=""),
          append = TRUE)
    
  }
}

#-------------------------------------------------------------------#
#                                                                   #
#              function_to_run_SKATBinary_Single                    #
#                                                                   #
#-------------------------------------------------------------------#

# SKATBinary_Single behaves slightly different from SKATBinary.  
# In particular:
# - it does not support missing_cutoff option 
# (which is defaulting to 0.15 in SKATBinary) and 
# - it requires explicit handling of cases with missed 
# covariates before calculating the null model

run_single_variant_test <- function(
  set.name, phenotypes, covariates, variant){
  
  # Initial num of cases and fraction of missed genotypes
  n_cases_init <- length(variant)
  missed_genotypes_fraction_init <- sum(is.na(variant)) / n_cases_init
  
  # In some tests SKATBinary_Single generated errors 
  # when there were cases with incomplete data on covariates
  incomplete.covariates <- is.na(apply(covariates,1,sum))
  
  phenotypes <- phenotypes[!incomplete.covariates]
  covariates <- covariates[!incomplete.covariates,]
  variant <- variant[!incomplete.covariates]
  
  # Num of retained cases
  n_cases_retained <- sum(!incomplete.covariates)
  
  # Calculate null model
  skat.null <- tryCatch(SKAT_Null_Model(phenotypes ~ covariates, Adjustment = FALSE, out_type="D"), error=function(e) e)
  
  # Calculate test
  skat.test <- tryCatch(SKATBinary_Single(variant, skat.null), error=function(e) e)
  
  # Compile result
  result <- list(
    set.name, 
    phenotypes, 
    covariates, 
    variant, 
    n_cases_init, 
    missed_genotypes_fraction_init, 
    n_cases_retained, 
    skat.null, 
    skat.test)
  
  names(result) <- c(
    "set.name", 
    "phenotypes", 
    "covariates", 
    "variant", 
    "n_cases_init",
    "missed_genotypes_fraction_init",
    "n_cases_retained",
    "skat.null",
    "single.variant.test")
  
  # Return result
  return(result)
  
}

#-------------------------------------------------------------------#
#                                                                   #
#          function_to_print_results_for_single_variant             #
#                                                                   #
#-------------------------------------------------------------------#

print_results_for_single_variant <- function(
  data, prefix, AC_UBC_init, AC_CBC_init, report_type = "summary_only", path=""){
  
  # Set name
  set <- data$set.name
  
  # Nums of initial variants, cases and fraction of missed genotypes
  n_cases_init <- data$n_cases_init
  n_vars_init <- 1
  missed_genotypes_fraction_init <- data$missed_genotypes_fraction_init
  
  # Num of retained cases
  n_cases_retained <- data$n_cases_retained
  
  # p-values
  p_value <- nm(data$single.variant.test$p.value)
  n_vars_tested <- nm(data$single.variant.test$param$n.marker)
  
  #summary_table_header <- c("set", "n_cases_init", "n_vars_init", "AC_UBC_init", "AC_CBC_init", 
  #  "missed_genotypes_fraction_init", "n_cases_retained", "n_vars_retained", "n_vars_retained_cr", 
  #  "burden_skatw", "burden_mbw", "burden_lw", "skat_skatw", "skat_mbw", "skat_lw", 
  #  "skato_skatw", "skato_mbw", "skato_lw", "cr_skat",  "cr_burden", "min_p")
  
  # Compile data for summary table 
  # NB: same p-value written to each test, as there is no need in aggregation,
  # so there should be no differences between the tests
  summary <- c(set, n_cases_init, n_vars_init, AC_UBC_init, AC_CBC_init, 
               missed_genotypes_fraction_init, 
               n_cases_retained, n_vars_tested, n_vars_tested, 
               p_value, p_value, p_value, 
               p_value, p_value, p_value, 
               p_value, p_value, p_value,
               p_value, p_value, p_value)
  
  # Write data to summary table
  write(paste(summary, sep="", collapse="\t"),
        file=paste(path, prefix, "_summary.txt", sep="", collapse=""),
        append = TRUE)
  
  # Write data to other tables, if required
  if (report_type == "full") {
    
    #details_table_header <- c("set", "n_cases_init", "n_vars_init", 
    #  "missed_genotypes_fraction_init", "n_cases_retained", 
    #  "p_val", "n_vars_total", "n_vars_tested", 
    #  "ma_count", "n_cases_with_ma", "method_for_p_assessment",
    #  "rare_vars_tested", "common_vars_tested", "cutoff")
    
    # Compile data
    details <- get_details(data$single.variant.test)
    single_variant_test <- c(set, n_cases_init, n_vars_init, 
                             missed_genotypes_fraction_init, 
                             n_cases_retained, details)
    
    # Write data to files
    
    write(paste(single_variant_test, sep="", collapse="\t"),
          file=paste(path, prefix, "_burden_skatw.txt", sep="", collapse=""),
          append = TRUE)
    
    write(paste(single_variant_test, sep="", collapse="\t"),
          file=paste(path, prefix, "_burden_mbw.txt", sep="", collapse=""),
          append = TRUE)
    
    write(paste(single_variant_test, sep="", collapse="\t"),
          file=paste(path, prefix, "_burden_lw.txt", sep="", collapse=""),
          append = TRUE)
    
    write(paste(single_variant_test, sep="", collapse="\t"),
          file=paste(path, prefix, "_skat_skatw.txt", sep="", collapse=""),
          append = TRUE)
    
    write(paste(single_variant_test, sep="", collapse="\t"),
          file=paste(path, prefix, "_skat_mbw.txt", sep="", collapse=""),
          append = TRUE)
    
    write(paste(single_variant_test, sep="", collapse="\t"),
          file=paste(path, prefix, "_skat_lw.txt", sep="", collapse=""),
          append = TRUE)
    
    write(paste(single_variant_test, sep="", collapse="\t"),
          file=paste(path, prefix, "_skato_skatw.txt", sep="", collapse=""),
          append = TRUE)
    
    write(paste(single_variant_test, sep="", collapse="\t"),
          file=paste(path, prefix, "_skato_mbw.txt", sep="", collapse=""),
          append = TRUE)
    
    write(paste(single_variant_test, sep="", collapse="\t"),
          file=paste(path, prefix, "_skato_lw.txt", sep="", collapse=""),
          append = TRUE)
    
    write(paste(single_variant_test, sep="", collapse="\t"),
          file=paste(path, prefix, "_cr_skat.txt", sep="", collapse=""),
          append = TRUE)
    
    write(paste(single_variant_test, sep="", collapse="\t"),
          file=paste(path, prefix, "_cr_burden.txt", sep="", collapse=""),
          append = TRUE)
    
  }
}
