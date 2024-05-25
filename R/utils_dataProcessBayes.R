#' @importFrom protDP dpc
#' @importFrom data.table ":=" uniqueN setnames
#' @importFrom MASS fitdistr
#' @import tidyverse
#' @import cmdstanr
#' @export
calculateMNARCurve = function(data){

  if ("H" %in% unique(data$LABEL)){
    data = data[data$LABEL != "H", ]
  }

  wide_data = data.table::dcast(data, FEATURE~RUN, value.var = "ABUNDANCE")
  wide_data[, FEATURE := NULL]

  dpc_params = dpc(wide_data)

  return(dpc_params$beta)
}



#' @export
MSstatsBayesSummarize = function(data, dpc_betas,
                                 bayes_method="MCMC",
                                 n_iterations=1000,
                                 chains=4, cores=4,
                                 elbo_samples=500,
                                 tol_rel_obj=.00001){

  keep = c("PROTEIN", "FEATURE", "RUN", "ABUNDANCE")
  subset_data = data[, ..keep]

  model_data = prepare_for_bayes(subset_data, dpc_betas)
  feature_data = model_data[[1]]
  stan_input = model_data[[2]]
  missing_runs = model_data[[3]]

  stan_file = system.file("stan", "missing_model.stan", package="MSstatsBayes")

  # Model
  # suppressWarnings({

  # path_to_opencl_lib <- "C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v10.2\\lib\\x64\\"
  # cpp_options = list(
  #   paste0("LDFLAGS+= -L\"",path_to_opencl_lib,"\" -lOpenCL")
  # )
  #
  # cmdstanr::cmdstan_make_local(cpp_options = cpp_options)
  # cmdstanr::rebuild_cmdstan()

  model = cmdstan_model(stan_file)#, force_recompile=TRUE,
                          # cpp_options = list(stan_opencl = TRUE,
                          #                    paste0("LDFLAGS+= -L\"",path_to_opencl_lib,"\" -lOpenCL")))
  # })

  fit = model$sample(
    data = stan_input,
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    refresh = 500 # print update every 500 iters
    # opencl_ids = c(0, 0), refresh = 0
    )

  # if (bayes_method == "MCMC"){
  #   fit = stan(file = stan_file, data = stan_input,
  #              chains = chains, iter = n_iterations, cores = cores,
  #              seed=100, model_name="MSstats_model")
  # } else if (bayes_method == "VB"){
  #   model = stan_model(stan_file)
  #   fit = vb(model, data = stan_input, iter = n_iterations,
  #            elbo_samples = elbo_samples, tol_rel_obj = tol_rel_obj)
  # }
  model_results = recover_data(fit, feature_data, missing_runs)

  return(model_results)
}



MSstatsBayesSummarizationOutput = function(summarized, input){

  # Define feature and protein data
  feature_data = input
  protein_data = summarized

  # Add missing value info to feature level data
  feature_data = merge(feature_data,
                       unique(protein_data[,c("PROTEIN_original",
                                              "RUN_original",
                                              "FEATURE_original",
                                              "imp_mean", "imp_sd"
                       )]),
                       by.x = c("PROTEIN", "RUN", "FEATURE"),
                       by.y = c("PROTEIN_original", "RUN_original",
                                "FEATURE_original"),
                       all.x=TRUE, all.y=FALSE)
  feature_data[imp_mean == 0, imp_mean:=NA]
  feature_data[imp_sd == 0, imp_sd:=NA]

  # Calculate summary stats for protein level data
  feature_data[, total_features := uniqueN(FEATURE), by = PROTEIN]
  feature_data[,"censored"] = !is.na(feature_data[,"imp_mean"])
  summary_stats = feature_data
  summary_stats[, NonMissingStats := MSstats:::.getNonMissingFilterStats(
    .SD, NULL)]
  summary_stats[, NumMeasuredFeature := sum(NonMissingStats),
                by = c("PROTEIN", "RUN")]
  summary_stats[, MissingPercentage := 1 - (NumMeasuredFeature / total_features)]
  summary_stats[, more50missing := MissingPercentage >= 0.5]
  summary_stats[, nonmissing_orig := LABEL == "L" & !censored]
  summary_stats[, NumImputedFeature := sum(LABEL == "L" & !nonmissing_orig),
                by = c("PROTEIN", "RUN")]
  summary_stats = summary_stats[, c("PROTEIN", "RUN",
                                    "NumMeasuredFeature", "MissingPercentage",
                                    "more50missing", "NumImputedFeature")]


  # Format feature data
  feature_data[ ,c("GROUP", "SUBJECT", "nonmissing", "n_obs",
                   "n_obs_run", "total_features", "prop_features") := NULL]

  setnames(feature_data,
           c("GROUP_ORIGINAL", "SUBJECT_ORIGINAL", "imp_mean", "imp_sd"),
           c("GROUP", "SUBJECT", "predicted", "predicted_sd"))

  feature_data[, newABUNDANCE := ifelse(!is.na(ABUNDANCE),
                                        ABUNDANCE, predicted)]

  feature_data = feature_data[, c("PROTEIN", "PEPTIDE", "TRANSITION",
                                  "FEATURE", "LABEL", "GROUP", "RUN",
                                  "SUBJECT", "FRACTION", "originalRUN",
                                  "censored", "INTENSITY", "ABUNDANCE",
                                  "newABUNDANCE", "predicted", "predicted_sd"
  )]


  # Add run level summarization
  summary_cols = c("PROTEIN_original", "RUN_original", "summarized", "errors")
  protein_data = merge(unique(protein_data[, ..summary_cols]),
                       unique(input[, c("PROTEIN", "RUN", "originalRUN",
                                        "GROUP_ORIGINAL", "SUBJECT_ORIGINAL")]
                       ),
                       by.x = c("PROTEIN_original", "RUN_original"),
                       by.y = c("PROTEIN", "RUN"),
                       all.x=TRUE, all.y=FALSE)

  protein_data = merge(protein_data, unique(summary_stats),
                       by.x = c("PROTEIN_original", "RUN_original"),
                       by.y = c("PROTEIN", "RUN"),
                       all.x=TRUE, all.y=FALSE)


  protein_data = merge(protein_data,
                       feature_data[, uniqueN(.SD),
                                    by = c("PROTEIN", "GROUP"),
                                    .SDcols = c("FEATURE", "originalRUN")],
                       by.x = c("PROTEIN_original", "GROUP_ORIGINAL"),
                       by.y=c("PROTEIN", "GROUP"),
                       all.x=TRUE, all.y=FALSE)


  setnames(protein_data,
           c("PROTEIN_original", "GROUP_ORIGINAL", "RUN_original",
             "summarized", "errors", "SUBJECT_ORIGINAL", "V1"),
           c("Protein", "GROUP", "RUN", "LogIntensities",
             "LogIntensities_sd", "SUBJECT", "TotalGroupMeasurements"))

  protein_data = protein_data[,c("RUN", "Protein", "LogIntensities",
                                 "LogIntensities_sd", "originalRUN", "GROUP",
                                 "SUBJECT", "TotalGroupMeasurements",
                                 "NumMeasuredFeature", "MissingPercentage",
                                 "more50missing", "NumImputedFeature")]


  return(list(FeatureLevelData = as.data.frame(feature_data),
              ProteinLevelData = as.data.frame(protein_data),
              SummaryMethod = "Bayesian")
  )
}

prepare_for_bayes = function(data, dpc_betas){

  data_list = remove_missing_runs(data)
  data = data_list[[1]]
  missing_data = data_list[[2]]

  format_data = format_data(data)
  priors = calculate_priors(data)
  results = flatten_input(format_data)

  input = results[[1]]
  lookup_table = results[[2]]

  stan_input = arrange_stan_data(input, priors, dpc_betas)

  return(list(lookup_table, stan_input, missing_data))
}

remove_missing_runs = function(data){

  data[, run_total:=all(is.na(ABUNDANCE)), by=.(PROTEIN, RUN)]

  remove = data[run_total != 0,][, run_total := NULL]
  keep = data[run_total == 0,][, run_total := NULL]

  return(list(keep, remove))

}

arrange_stan_data = function(input, priors, dpc_betas){

  num_runs = length(unique(input$Protein_run_idx))
  num_feat = length(unique(input$Protein_feature_idx))

  ## TODO: IDK Dummy matrix problem when the data is huge (is this a problem?)
  # obs_mat = matrix(as.numeric(unlist(input[which(!is.na(input[,3])),])),
  #                  nrow=nrow(input[which(!is.na(input[,3])),]))
  # missing_mat = matrix(as.numeric(unlist(input[which(is.na(input[,3])),])),
  #                      nrow=nrow(input[which(is.na(input[,3])),]))
  # if (all(dim(missing_mat) == c(0,0))) {
  #   missing_mat = matrix(c(0,0,0,0), nrow=1)
  # }
  ## TODO: Fix this its so bad to load mid package
  library(fastDummies)
  temp = input[, c("Protein_run_idx", "Protein_feature_idx")]
  temp[] = lapply(temp, as.character)
  X = dummy_cols(temp) %>%
    dplyr::select(-Protein_run_idx, -Protein_feature_idx) %>% as.matrix()
  missing_X = X[is.na(input$ABUNDANCE), ]

  missing_idx = which(is.na(input$ABUNDANCE))

  obs = input$ABUNDANCE
  obs[is.na(obs)] = 0

  stan_input = list(N = length(input$ABUNDANCE),
                    N_missing = sum(is.na(input$ABUNDANCE )),
                    y = obs,
                    missing_idx = missing_idx,
                    N_R = num_runs,
                    N_F = num_feat,
                    X_R = X[,1:num_runs],
                    X_F = X[,(num_runs+1):(num_runs+num_feat)],
                    X_R_missing = missing_X[,1:num_runs],
                    X_F_missing = missing_X[,(num_runs+1):(num_runs+num_feat)],
                    R_mu_priors = priors$run_priors$R_mu_prior,
                    R_sd_priors = priors$run_priors$R_sd_prior,
                    F_mu_priors = priors$feature_priors$F_mu_prior,
                    F_sd_priors = priors$feature_priors$F_sd_prior,
                    sigma_mu_priors = priors$sigma_priors$sigma_mu_prior,
                    sigma_sd_priors = priors$sigma_priors$sigma_sd_prior
  )
}

format_data = function(data){

  # Missing indicator
  data[, "Missing"] = ifelse(is.na(data[, "ABUNDANCE"]), 1., 0.)

  for (col in c("PROTEIN", "RUN", "FEATURE")){
    id_map = as.numeric(droplevels(unlist(data[, ..col])))
    data[, paste0(col, "_original")] = data[,..col]
    data[,col] = id_map
  }

  # order
  data = data[with(data, order(PROTEIN, FEATURE, RUN)), ]

  return(data)
}

flatten_input = function(data) {

  # Convert the data to a data frame
  # data = data.frame(data, stringsAsFactors = FALSE)
  # colnames(data) = c("Protein", "Run", "Feature", "Intensity", "Missing")

  # Add a list index
  # data$list_index = 1:nrow(data)

  # Create a unique protein_run identifier
  data$Protein_run = paste(data$PROTEIN, data$RUN, sep = "_")

  # Create a unique protein_feature identifier
  data$Protein_feature = paste(data$PROTEIN, data$FEATURE, sep = "_")

  # Create protein_run_idx and protein_feature_idx
  protein_run_idx = unique(data$Protein_run)
  protein_feature_idx = unique(data$Protein_feature)

  data$Protein_run_idx = match(data$Protein_run, protein_run_idx)
  data$Protein_feature_idx = match(data$Protein_feature, protein_feature_idx)

  # Return the flattened data
  flatten_data = data[, c("Protein_run_idx", "Protein_feature_idx",
                          "ABUNDANCE", "PROTEIN")]

  return(list(flatten_data, data))
}

#' Calculate priors given some data for the run, feature, and error effect
#'
#' @param data MSstats formatted dataframe
#' @param run_only logical. If TRUE, only calculate priors for the run effect
#'
#' @return list of priors for the run, feature, and error effect
calculate_priors = function(data, run_only=FALSE){
  # snag some priors
  # Run
  run_summary = data %>% group_by(PROTEIN, RUN) %>%
    summarize(run_mean = mean(ABUNDANCE, na.rm=TRUE),
              run_sd = sd(ABUNDANCE, na.rm=TRUE))
  run_priors = fitdistr(run_summary$run_mean, "normal")
  M = log(run_priors[[1]][[2]]^2 /
            sqrt(run_priors[[1]][[2]]^2 + run_priors[[2]][[2]]^2))
  S = sqrt(log(1 + (run_priors[[2]][[2]]^2 / run_priors[[1]][[2]]^2)))

  run_priors = list(R_mu_prior = c(run_priors[[1]][[1]],
                                   run_priors[[2]][[1]]),
                    R_sd_prior = c(M,S))

  # Feature
  feature_summary = data %>% merge(run_summary, by = c("PROTEIN", "RUN")) %>%
    mutate(feature_intensity = ABUNDANCE - run_mean) %>% group_by(FEATURE) %>%
    summarize(feature_mean = mean(feature_intensity, na.rm=TRUE),
              feature_sd = sd(feature_intensity, na.rm=TRUE))

  feature_priors = fitdistr(feature_summary$feature_mean, "normal")
  M = log(feature_priors[[1]][[2]]^2 /
            sqrt(feature_priors[[1]][[2]]^2 + feature_priors[[2]][[2]]^2))
  S = log(1 + (feature_priors[[2]][[2]]^2 / feature_priors[[1]][[2]]^2))

  feature_priors = list(F_mu_prior = c(feature_priors[[1]][[1]],
                                       feature_priors[[2]][[1]]),
                        F_sd_prior = c(M, S))

  # sigma
  sigma_summary = data %>%
    merge(run_summary, by = c("PROTEIN", "RUN")) %>%
    merge(feature_summary, by = "FEATURE") %>%
    mutate(error = ABUNDANCE - run_mean - feature_mean) %>%
    filter(!is.na(error)) %>%
    group_by(RUN) %>%
    summarize(sigma_sd = sd(error, na.rm=TRUE))

  sigma_priors = fitdistr(na.omit(sigma_summary$sigma_sd), "lognormal")

  M = log(sigma_priors[[1]][[2]]^2 /
            sqrt(sigma_priors[[1]][[2]]^2 + sigma_priors[[2]][[2]]^2))
  S = log(1 + (sigma_priors[[2]][[2]]^2 / sigma_priors[[1]][[2]]^2))

  sigma_priors = list(sigma_mu_prior = c(sigma_priors[[1]][[1]],
                                         sigma_priors[[2]][[1]]),
                      sigma_sd_prior = c(M, S))

  priors = list(run_priors=run_priors,
                feature_priors=feature_priors,
                sigma_priors=sigma_priors)

  return(priors)
}

extract_values = function(fit, n_proteins, n_runs){

  # Subset model into runs and sigmas
  # model_overview = summary(fit)$summary
  model_overview = fit$summary()
  run_results = model_overview[grepl("R\\[", model_overview$variable), ]
  sigma_results = model_overview[grepl("sigma\\[", model_overview$variable), "mean"]

  # Grab summaries
  summarized_values = run_results$mean

  # Grab errors
  R_errors = run_results$sd
  errors = sqrt(sigma_results$mean^2 + R_errors^2)

  return(list("summarized" = summarized_values,
              "errors" = errors))
}

recover_data = function(fit, feature_data, missing_runs){

  # parameters = c("run", "sigma", "beta0", "beta1")#"obs_mis",

  # Missing imputation
  results = fit$summary()[grepl("y_impute\\[", fit$summary()$variable),
                          c("mean", "sd")]
  imp_index = which(is.na(feature_data$ABUNDANCE))

  feature_data[, "imp_mean"] = 0
  feature_data[, "imp_sd"] = 0
  feature_data[imp_index, "imp_mean" := results[,1]]
  feature_data[imp_index, "imp_sd" := results[,2]]


  # Feature estimation
  # results = rstan::summary(fit, pars="feature")
  # results = as.data.frame(results$summary[,c("mean", "sd")])
  # results = rename(results, c("feature_mean" = "mean",
  #                             "feature_sd" = "sd"))
  # results$Protein_feature_idx = seq_len(nrow(results))
  # feature_data = merge(feature_data, results, all.x=TRUE,
  #                      by="Protein_feature_idx")
  n_proteins = nrow(unique(feature_data[, "PROTEIN"]))
  n_runs = nrow(unique(feature_data[, c("PROTEIN", "RUN")]))

  model_output = extract_values(fit, n_proteins, n_runs)

  # Run estimation
  results = as.data.frame(model_output)
  results$Protein_run_idx = seq_len(nrow(results))
  feature_data = merge(feature_data, results, all.x=TRUE,
                       by="Protein_run_idx")


  # beta0 = rstan::summary(fit, pars="beta0")
  # beta1 = rstan::summary(fit, pars="beta1")
  # run = rstan::summary(fit, pars="run")
  # feature = rstan::summary(fit, pars="feature")
  # obs_mis = rstan::summary(fit, pars="obs_mis")
  # sigma = rstan::summary(fit, pars="sigma")

  return(feature_data)
    # list(result_df = feature_data,
              # bayes_results=list(
              #   # beta = rbind(beta0$summary, beta1$summary),
              #   run = run$summary,
              #   feature = feature$summary,
              #   # missing = obs_mis$summary,
              #   sigma = sigma$summary)
}
