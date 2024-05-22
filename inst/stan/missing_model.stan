functions{
    vector merge_missing( array[] int miss_indexes , vector x_obs , vector x_miss ) {
        int N = dims(x_obs)[1];
        int N_miss = dims(x_miss)[1];
        vector[N] merged;
        merged = x_obs;
        for ( i in 1:N_miss )
            merged[ miss_indexes[i] ] = x_miss[i];
        return merged;
    }
}
data {
    int<lower=0> N;
    int<lower=0> N_missing;

    vector[N] y;
    // int missing_idx[N_missing];
    array[N_missing] int missing_idx;
    
    int<lower=0> N_R;
    int<lower=0> N_F;

    matrix[N, N_R] X_R;
    matrix[N, N_F] X_F;
    
    matrix[N_missing, N_R] X_R_missing;
    matrix[N_missing, N_F] X_F_missing;
    
    // Priors
    vector[2] R_mu_priors;
    vector[2] R_sd_priors;
    
    vector[2] F_mu_priors;
    vector[2] F_sd_priors;
    
    vector[2] sigma_mu_priors;
    vector[2] sigma_sd_priors;
}
parameters {
    real R_mu;
    real<lower=0> R_sd;
    vector[N_R] R;
    // vector[N_R] R_imp;

    real F_mu;
    real<lower=0> F_sd;
    vector[N_F] F;

    real sigma_mu;
    real<lower=0> sigma_sd;
    vector<lower=0>[N_R] sigma;
    
    vector[N_missing] y_impute;
    // vector<lower=0>[N_R] shift;
}
// transformed parameters {
// 
//     // real shift;
//     vector<lower=0>[N_missing] shift;
//     //.4 is the hard coded slope in the sim function. Swap to estimated slope irl.
//     // shift = .4*exp(sigma_mu_priors[1] + (sigma_mu_priors[2]/2));
//     // shift = .4*(sigma_mu_priors[1]);
//     shift = ;
// }
model {
    
    vector[N] y_merge;
    
    R_mu ~ normal(R_mu_priors[1], R_mu_priors[2]);
    R_sd ~ lognormal(R_sd_priors[1], R_sd_priors[2]);
    R ~ normal(R_mu, R_sd);
    // R_imp ~ normal(R_mu, R_sd);

    F_mu ~ normal(F_mu_priors[1], F_mu_priors[2]);
    F_sd ~ lognormal(F_sd_priors[1], F_sd_priors[2]);
    F ~ normal(F_mu, F_sd);

    sigma_mu ~ normal(sigma_mu_priors[1], sigma_mu_priors[2]);
    sigma_sd ~ lognormal(sigma_sd_priors[1], sigma_sd_priors[2]);
    sigma ~ lognormal(sigma_mu, sigma_sd);
    
    y_impute ~ normal((X_R_missing*R) + (X_F_missing*F) - (.4 .* (X_R_missing*sigma)), X_R_missing*sigma);
    
    y_merge = merge_missing(missing_idx, to_vector(y), y_impute);
    // y ~ normal((X_R*R) + (X_F*F), X_R*sigma);
    y_merge ~ normal((X_R*R) + (X_F*F), X_R*sigma);
    
    // y_mar ~ normal((X_R_missing*R) + (X_F_missing*F), X_R_missing*sigma);
    
    // missing_reason ~ bernoulli(mnar_prob);
    // if (missing_reason){
    //     y_missing ~ normal((X_R_missing*R) + (X_F_missing*F) - (.4*(X_R_missing*error)), X_R_missing*error);
    // } else {
    //     y_missing ~ normal((X_R_missing*R) + (X_F_missing*F), X_R_missing*error);
    // }
}
