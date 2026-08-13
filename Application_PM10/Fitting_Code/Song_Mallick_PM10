rm(list=ls(all=TRUE))
library("fda")        
library("MASS")     
library("ggplot2")
library("dplyr")
library("tidyr")
library("posterior")
library("bayesplot")
library("coda")
library("wavethresh")
library("cmdstanr")
library("snowfall")

setwd("/home/alextk4/Ap")
load("PM10_2025_day.RData")
load("coordinates_utm_full_km.RData")
load("coordinates_utm_obs_km.RData")

#-------------------------------------------------------------------------------
Sample  = PM10_2025_day %>% 
          dplyr::select(-Day)

Sample  = apply(Sample, c(1, 2), function(x) ifelse(is.na(x), NA, log(x)))
Sample  = as.data.frame(Sample)

log_sample = Sample %>%
             dplyr::select(-ACO, -IZT, -TAH)

Yobs_pred  = Sample %>%
            dplyr::select(ACO, IZT, TAH)

m_obs  = dim(log_sample)[2]
m_pred = dim(Yobs_pred)[2]

# (1) Global standardization

matriz_Y = as.matrix(log_sample)
mu_Y     = mean(matriz_Y)
sd_Y     = sd(matriz_Y)

# Standardize the entire matrix while keeping the spatial structure intact
Sample_padronizado = (matriz_Y - mu_Y) / sd_Y
Sample_padronizado = as.data.frame(Sample_padronizado)

# ==============================================================================
# DYNAMIC PADDING (MIRRORING) - identical to the proposed model
# ==============================================================================
# Now 'n' takes the value of the inflated size (n_padded) for Stan to process
n_original = dim(Sample_padronizado)[1] 

# Compute how much needs to be added on each side (half the signal)
pad_size = n_original / 2 

# The new total size will be exactly double the original (e.g., 128 -> 256)
n_padded = n_original * 2 
Y_padded = matrix(NA, nrow = n_padded, ncol = m_obs)

for(j in 1:m_obs) {
  curva = Sample_padronizado[, j]
  # Perfect mirroring to trick the periodic algorithm
  Y_padded[, j] = c(rev(curva[1:pad_size]), curva, rev(curva[(pad_size + 1):n_original]))
}

Sample_padronizado = as.data.frame(Y_padded)

n = dim(Sample_padronizado)[1]
L = log2(n)                                             # Maximum resolution level
n_l = c(2, 2^(1:(L - 1)));sum(n_l)                      # number of coefficients per level n_l 
level_id = rep(0:(L - 1), times = n_l);length(level_id) # Map each of the n positions in d_matrix to its level 


# ------------------------------------------------------------------------------
# SCALE/DETAIL INDICATOR
is_detail = c(0, rep(1, n - 1))

# Spatial prior variant: 1=GG, 2=GL, 3=LG, 4=LL
prior_type = 4

# Hyperparameters of Song & Mallick's Inverse-Gamma priors: IG(2, 1)
ig_shape = 2.0
ig_rate  = 1.0
jitter_chol = 1e-6

locations_obs  = coordinates_utm_obs_km
locations_full = coordinates_utm_full_km

dist_mat = as.matrix(dist(locations_full[, c("Easting", "Northing")]))
D_max = max(dist_mat)

# ==============================================================================
# COORDINATE RESCALING (CRITICAL FOR A FAIR COMPARISON)
#
# This is the ONLY conceptual change relative to S&M's literal specification,
# and it is deliberate. It amounts to reparameterizing the prior on phi.
#
# Quantitative justification:
#   The DGP uses Sigma = sigma_s^2 * exp(-d / phi_s) with phi_s = 0.4, i.e.,
#   true decay = 2.5 on the unit square.
#
#   WITHOUT rescaling: phi ~ IG(2,1) has mean 1 against a truth of 2.5. Since
#   P(phi > 2.5) = P(Gamma(2,1) < 0.4) ~ 0.062, the true value falls in the
#   94th percentile of the prior. The prior would pull S&M toward overly long
#   ranges -> oversmoothing -> artificially poor performance, and the
#   comparison would be won by a calibration artifact.
#
#   WITH rescaling (factor 4.89 / 1.414 = 3.458): the true decay in the new
#   units becomes 2.5 / 3.458 = 0.72, against a prior mean of 1. Practically
#   centered.
#
# Also run with rescale_coords = FALSE as a sensitivity analysis and report
# both in the manuscript.
# ==============================================================================
rescale_coords = TRUE
D_max_SM = 4.89

if (rescale_coords) {
  escala_SM      = D_max_SM / D_max
  locations_obs  = locations_obs  * escala_SM
  locations_full = locations_full * escala_SM
  D_max          = D_max * escala_SM
}

locations_obs  = as.matrix(locations_obs)
locations_full = as.matrix(locations_full)
# ------------------------------------------------------------------------------
# (2). Discrete Wavelet Transform using Daubechies
# ------------------------------------------------------------------------------
wavelet_list = vector("list", m_obs)
for(j in 1:m_obs) {
  wd_obj = wavethresh::wd(Sample_padronizado[, j], 
                          filter.number = 4, 
                          family = "DaubExPhase",
                          bc="periodic") # Keep periodic so as not to break the 2^j tree
  wavelet_list[[j]] = wd_obj
}

# ------------------------------------------------------------------------------
# (3). Function to extract and organize the coefficients into the format 
# ------------------------------------------------------------------------------
coef_empirical = function(wavelet_list, n, L) {
  m = length(wavelet_list)
  d_empirical = matrix(NA, nrow = n, ncol = m_obs)
  
  for (j in seq_len(m_obs)) {
    w = wavelet_list[[j]]
    c0 = accessC(w, level = 0)     
    d0 = accessD(w, level = 0)      
    vetor_organizado = c(c0, d0)
    
    for (i in 1:(L - 1)) {
      di = accessD(w, level = i)    
      vetor_organizado = c(vetor_organizado, di)
    }
    
    if (length(vetor_organizado) != n) {
      stop(sprintf("Incorrect length at j = %d: got %d, expected %d", j, length(vetor_organizado), n))
    }
    d_empirical[, j] = vetor_organizado
  }
  return(d_empirical)
}

d_empirical = t(coef_empirical(wavelet_list = wavelet_list, n = n, L = L))

###############################
### Monte Carlo Simulation ####
###############################

song_mallick_model="
data {
  int<lower=1> m_obs; 
  int<lower=1> m_pred;                    
  int<lower=1> n;                                                             
  matrix[m_obs, n] d_empirical;
  matrix[m_obs, 2] locations_obs;                         // matrix of observed coordinates
  matrix[(m_obs + m_pred), 2] locations_full;             // matrix of observed + predicted coordinates
  array[n] int<lower=0, upper=1> is_detail;               // 1 = detail (beta), 0 = scale (alpha)
  int<lower=1, upper=4> prior_type;                       // 1=GG, 2=GL, 3=LG, 4=LL
  real<lower=0> ig_shape;                                 // 2.0  (Song & Mallick: IG(2,1))
  real<lower=0> ig_rate;                                  // 1.0
  real<lower=0> jitter_chol;                              // 1e-6
}

transformed data {
  array[m_obs] vector[2] s_locs;
  array[m_obs + m_pred] vector[2] all_locs;
  array[n] int use_laplace;
  
  for (i in 1:m_obs) {
    s_locs[i] = to_vector(locations_obs[i]);
  }
  
  for (i in 1:(m_obs + m_pred)) {
  all_locs[i] = to_vector(locations_full[i]);
  }
  
  // Select which coefficients receive the Laplace prior (scale mixture)
  // according to the chosen variant
  for (j in 1:n) {
    if (prior_type == 1)      use_laplace[j] = 0;                    // GG
    else if (prior_type == 2) use_laplace[j] = is_detail[j];         // GL
    else if (prior_type == 3) use_laplace[j] = 1 - is_detail[j];     // LG
    else                      use_laplace[j] = 1;                    // LL
  }
}

parameters {
  real<lower=0> sigma2_obs;                                                   
  matrix[m_obs, n] z;                                                             
  vector<lower=0>[n] kappa;                               // spatial variation (sill), one per coefficient
  vector<lower=0>[n] phi;                                 // spatial decay, one per coefficient
  vector<lower=0>[n] lambda;                              // Laplace mixture parameter, one per coefficient
}

transformed parameters {
  matrix[m_obs, n] theta;                                                         
  vector<lower=0>[n] coef_scale;
  
{
  for (j in 1:n) {
  // Song & Mallick exponential covariance: kappa * exp(-phi * d).
  // gp_exponential_cov(x, sigma, length_scale) returns sigma^2 * exp(-d/length_scale),
  // so the decay phi corresponds to length_scale = 1/phi. We build the
  // CORRELATION matrix (magnitude 1) and keep the scale in coef_scale,
  // exactly as in the proposed model.
  //
  // NOTE: this convention is RECIPROCAL to the proposed model's, where phi
  // enters directly as the length scale in gp_matern32_cov. When comparing
  // estimates from the two models, convert before tabulating.
  matrix[m_obs, m_obs] R = gp_exponential_cov(s_locs, 1, inv(phi[j]));
  for (i in 1:m_obs) {R[i, i] += jitter_chol;}
  matrix[m_obs, m_obs] L_cov = cholesky_decompose(R);
  
  if (use_laplace[j] == 1) {
      coef_scale[j] = sqrt(kappa[j] / lambda[j]);                    // Laplace prior
    }
  else {
      coef_scale[j] = sqrt(kappa[j]);                                // Gaussian prior
    }
  
  theta[, j] = coef_scale[j] * (L_cov * z[, j]);
   }
 }
}
  
model {
sigma2_obs ~ inv_gamma(ig_shape, ig_rate);
kappa ~ inv_gamma(ig_shape, ig_rate);
phi ~ inv_gamma(ig_shape, ig_rate);

// Laplace mixture parameter: Eq. (5) of Song & Mallick.
// Sampled for all coefficients to keep the block rectangular; where
// use_laplace = 0 it does not enter theta and merely samples from its own
// prior. When building the parameter-count table, count only the ACTIVE ones
// (sum(use_laplace)), otherwise the S&M model is overstated.
lambda ~ inv_gamma((m_obs + 1) / 2.0, 0.25);

to_vector(z) ~ std_normal();

to_vector(d_empirical) ~ normal(to_vector(theta), sqrt(sigma2_obs));
}

generated quantities {
  
  real<lower=0> sigma_obs = sqrt(sigma2_obs);
  
  matrix[m_obs, n] log_lik; 

  for (j in 1:m_obs) {
  for (i in 1:n) {
  log_lik[j, i] = normal_lpdf(d_empirical[j, i] | theta[j,i], sigma_obs);
  }
 }
  
// 2. spatial prediction
   
  matrix[m_pred, n] theta_pred; // Coefficients at the new locations
   
  for (j in 1:n) {
    
  matrix[(m_obs + m_pred),(m_obs + m_pred) ] C_full = gp_exponential_cov(all_locs, 1, inv(phi[j]));
  
  for (k in 1:(m_obs + m_pred)){C_full[k, k] += jitter_chol;} 
  
//  Partition the matrix
  
  matrix[m_obs, m_obs]   Sigma_11 = C_full[1:m_obs, 1:m_obs];                     // Obs x Obs
  matrix[m_obs, m_pred]  Sigma_12 = C_full[1:m_obs, (m_obs + 1):(m_obs + m_pred)];               // Obs x Pred
  matrix[m_pred, m_pred] Sigma_22 = C_full[(m_obs + 1):(m_obs + m_pred), (m_obs + 1):(m_obs + m_pred)];         // Pred x Pred
      
  matrix[m_obs, m_obs] L_Sigma11 = cholesky_decompose(Sigma_11);
  
  matrix[m_obs, m_pred] K_div = mdivide_left_tri_low(L_Sigma11, Sigma_12);
  
  // Conditional covariance (correlation only)
  matrix[m_pred, m_pred] Sigma_cond = Sigma_22 - K_div' * K_div;
      
  // Enforce symmetry 
      for(r in 1:m_pred) {
        for(col_idx in (r+1):m_pred) Sigma_cond[r, col_idx] = Sigma_cond[col_idx, r];
        Sigma_cond[r, r] += jitter_chol; 
      }
      
  // Sample the predicted coefficients
  vector[m_obs] theta_obs_k = col(theta, j);
  
  // Conditional mean
  vector[m_obs] alpha_vec = mdivide_left_tri_low(L_Sigma11, theta_obs_k);
  vector[m_pred] mu_cond = K_div' * alpha_vec;
  
  matrix[m_pred, m_pred] Cov_cond_scaled = square(coef_scale[j]) * Sigma_cond;
  theta_pred[, j] = multi_normal_rng(mu_cond, Cov_cond_scaled);
  }
}
"

Model=write_stan_file(song_mallick_model)
mod=cmdstan_model(Model)

# 4) Group the data and hyperparameters into a list for Stan
data_list = list(
  m_obs                 = m_obs,
  m_pred                = m_pred,
  n                     = n,
  locations_obs         = locations_obs,
  locations_full        = locations_full,
  d_empirical           = d_empirical,
  is_detail             = is_detail,
  prior_type            = prior_type,
  ig_shape              = ig_shape,
  ig_rate               = ig_rate,
  jitter_chol           = jitter_chol
)

# Timing for the computational-cost table (Table 3 of the manuscript)
tempo_inicio = Sys.time()

fit_model = mod$sample(
  data            = data_list,
  seed            = 123,
  chains          = 4,
  parallel_chains = 1, 
  iter_warmup     = 2000,
  iter_sampling   = 1000,
  adapt_delta     = 0.99,
  max_treedepth   = 12
)

tempo_total = difftime(Sys.time(), tempo_inicio, units = "mins")
print(tempo_total)
