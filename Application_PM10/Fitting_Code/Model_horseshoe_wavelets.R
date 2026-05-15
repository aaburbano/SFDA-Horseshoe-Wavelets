rm(list=ls(all=TRUE))

library("MASS")    
library("fields")  
library("ggplot2")
library("dplyr")
library("tidyr")
library("wavethresh")
library("posterior")
library("bayesplot")
library("coda")
library("matrixcalc")
library("cmdstanr")

load("../Main_Data/PM10_2025_day.RData")
load("../Main_Data/coord.RData")

load("../Main_Data/coordinates_utm_full_km.RData")
load("../Main_Data/coordinates_utm_obs_km.RData")
#-------------------------------------------------------------------------------
Sample  = PM10_2025_day %>% 
          dplyr::select(-Day)

# Apply log transformation to the sample matrix
Sample = apply(Sample, c(1, 2), function(x) ifelse(is.na(x), NA, log(x)))
Sample = as.data.frame(Sample)

log_sample = Sample %>%
             dplyr::select(-ACO, -IZT, -TAH)

Yobs_pred = Sample %>%
            dplyr::select(ACO, IZT, TAH)

m_obs  = dim(log_sample)[2]
m_pred = dim(Yobs_pred)[2]

# (1) Global Standardization

matriz_Y = as.matrix(log_sample)
mu_Y     = mean(matriz_Y)
sd_Y     = sd(matriz_Y)

# Standardize the entire matrix while keeping the spatial structure intact
Sample_padronizado = (matriz_Y - mu_Y) / sd_Y
Sample_padronizado = as.data.frame(Sample_padronizado)

#save(Yobs_pred, file = "C:/Users/alexb/OneDrive/Artigo2/Dados_Sinteticos/Simulation_Horseshoe/Ap_PM10_2/Data_PM10/teste_2/Yobs_pred.RData")
# ==============================================================================
# CHANGE 1: DYNAMIC PADDING (MIRRORING)
# ==============================================================================
n_original = dim(Sample_padronizado)[1] 
         
# Calculate padding size for each side (half of the signal)
pad_size = n_original / 2 

# New total size will be exactly double the original (e.g., 128 -> 256)
n_padded = n_original * 2 
Y_padded = matrix(NA, nrow = n_padded, ncol = m_obs)

for(j in 1:m_obs) {
  curva = Sample_padronizado[, j]
  # Perfect mirroring to account for the periodic algorithm behavior
  Y_padded[, j] = c(rev(curva[1:pad_size]), curva, rev(curva[(pad_size + 1):n_original]))
}

Sample_padronizado = as.data.frame(Y_padded)
# ==============================================================================
# END OF CHANGE 1
# ==============================================================================

# 'n' now takes the value of the inflated size (n_padded) for Stan processing
n = dim(Sample_padronizado)[1]
L = log2(n)                                                                 # Maximum resolution level
coarse_levels = floor(log2(log(n))) + 1                                     # Coarse levels
shrunk_levels = L - coarse_levels                                           # Number of levels undergoing shrinkage
n_l = c(2, 2^(1:(L - 1))); sum(n_l)                                         # Number of coefficients per level (n_l) 

n_shrunk = sum(n_l[(coarse_levels + 1):L])                                  # Total number of coefficients undergoing shrinkage
p0 = round(sapply((coarse_levels + 1):L, function(l) max(1, 0.1 * n_l[l])), 0) # Expected number of non-zero coefficients per level (p0) 
level_id = rep(0:(L - 1), times = n_l); length(level_id)                    # Map each of the n positions in d_matrix to its level 

# Pre-compute id_l in R
max_nl = max(n_l)
id_l   = array(0, dim = c(L, max_nl))

for (l in 1:L) {
  cols = which(level_id == (l - 1))
  id_l[l, 1:length(cols)] = cols
}


locations_obs  = coordinates_utm_obs_km
locations_full = coordinates_utm_full_km

dist_mat = as.matrix(dist(locations_full[, c("Easting", "Northing")]))
D_max = max(dist_mat)
slab_scale = 2.0
slab_df = 4.0

locations_obs  = as.matrix(locations_obs)
locations_full = as.matrix(locations_full)
# ------------------------------------------------------------------------------
# (2) Discrete Wavelet Transform using Daubechies
# ------------------------------------------------------------------------------
wavelet_list = vector("list", m_obs)
for(j in 1:m_obs) {
  wd_obj = wavethresh::wd(Sample_padronizado[, j], 
                          filter.number = 4, 
                          family = "DaubExPhase",
                          bc="periodic") # Keep periodic to maintain the 2^j tree structure
  wavelet_list[[j]] = wd_obj
}

# ------------------------------------------------------------------------------
# (3) Function to extract and organize coefficients format
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
### Monte Carlo Simulation  ###
###############################

hierarchical_wavelet_model="
data {
  int<lower=1> m_obs; 
  int<lower=1> m_pred;                    
  int<lower=1> n;                                                             
  int<lower=1> L;                                                             
  int<lower=1> shrunk_levels;                                                 
  int<lower=1> coarse_levels;                                                  
  array[L] int n_l;                                                           
  int<lower=1> max_nl;                                                        
  array[L,max_nl] int id_l;                                                   
  matrix[m_obs, n] d_empirical;
  matrix[m_obs, 2] locations_obs;                         // Observed coordinates matrix 
  matrix[(m_obs + m_pred), 2] locations_full;             // Observed + Prediction coordinates matrix
  real<lower=0> D_max;                                                        
  int<lower=1> n_shrunk;                                                      
  vector<lower=0>[shrunk_levels] p0;                                          
  real<lower=0> slab_scale;                                                   
  real<lower=0> slab_df;                                                      
}

transformed data {
  array[m_obs] vector[2] s_locs;
  array[m_obs + m_pred] vector[2] all_locs;
  
  for (i in 1:m_obs) {
    s_locs[i] = to_vector(locations_obs[i]);
  }
  
  for (i in 1:(m_obs + m_pred)) {
  all_locs[i] = to_vector(locations_full[i]);
  }
}

parameters {
  real<lower=0> sigma_obs;                                                    
  matrix[m_obs, n] z;                                                             
  vector<lower=0>[coarse_levels] sigma_coarse;                                
  vector<lower=0>[L] phi;                                                     
  vector<lower=0>[shrunk_levels] aux1_global;
  vector<lower=0>[shrunk_levels] aux2_global;
  vector<lower=0>[n_shrunk] aux1_local;
  vector<lower=0>[n_shrunk] aux2_local;
  vector<lower=0>[shrunk_levels] caux;
}

transformed parameters {
  matrix[m_obs, n] theta;                                                         
  vector<lower=0>[shrunk_levels] tau_global;
  vector<lower=0>[shrunk_levels] c;
  vector<lower=0>[n] coef_scale;
  
{
  int local_idx = 1;                                                          
  
  for (l in 1:L) {
  int nl = n_l[l];
  if (nl > 0) {
  matrix[m_obs, m_obs] R = gp_matern32_cov(s_locs, 1, phi[l]);                        
  for (i in 1:m_obs) {R[i, i] += 1e-8;}                                          
  matrix[m_obs, m_obs] L_cov = cholesky_decompose(R);
  
  array[nl] int idxs = id_l[l, 1:nl];                                         
                                                                              
  if (l <= coarse_levels) {                                                   
      for (k in 1:nl) {
      int col_idx = idxs[k];
      coef_scale[col_idx] = sigma_coarse[l];
      theta[, col_idx] = sigma_coarse[l] * (L_cov * z[, col_idx]);            
        }
      }
  else {                                                                      
      int sl = l - coarse_levels;                                             
      real tau0 = (p0[sl] / (nl - p0[sl])) * (sigma_obs / sqrt(m_obs * 1.0));
      tau_global[sl] = aux1_global[sl] * sqrt(aux2_global[sl]) * tau0;        
      c[sl] = slab_scale * sqrt(caux[sl]);                                    
      
  for (k in 1:nl) {
      int col_idx = idxs[k];
      real lambda_k = aux1_local[local_idx] * sqrt(aux2_local[local_idx]);    
      
      real lambda_tilde_k = sqrt( (c[sl]^2 * square(lambda_k)) / (c[sl]^2 + tau_global[sl]^2 * square(lambda_k)) );
      coef_scale[col_idx] = tau_global[sl] * lambda_tilde_k;
      theta[, col_idx] = (tau_global[sl] * lambda_tilde_k) * (L_cov * z[, col_idx]);
      local_idx += 1;
        }
       }
     }
   }
 }
}
  
model {
sigma_obs ~ student_t(2, 0 , 1);
sigma_coarse ~ student_t(2, 0 , 1);
phi ~ lognormal(log(0.35 * D_max), 1);

to_vector(z) ~ std_normal();

aux1_global ~ std_normal();
aux2_global ~ inv_gamma(0.5, 0.5); 

aux1_local ~ std_normal();
aux2_local ~ inv_gamma(0.5, 0.5); 

caux ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);

to_vector(d_empirical) ~ normal(to_vector(theta), sigma_obs);
}

generated quantities {
  
  matrix[m_obs, n] log_lik; 

  for (j in 1:m_obs) {
  for (i in 1:n) {
  log_lik[j, i] = normal_lpdf(d_empirical[j, i] | theta[j,i], sigma_obs);
  }
 }
  
// 2. Spatial prediction
   
  matrix[m_pred, n] theta_pred; // Coefficients at new locations
   
  for (l in 1:L) {
  int nl = n_l[l];
    
  if (nl > 0) {
    
  matrix[(m_obs + m_pred),(m_obs + m_pred) ] C_full = gp_matern32_cov(all_locs, 1, phi[l]);
  
  for (k in 1:(m_obs + m_pred)){C_full[k, k] += 1e-8;} 
  
//  Matrix Partitioning
  
  matrix[m_obs, m_obs]   Sigma_11 = C_full[1:m_obs, 1:m_obs];                     // Obs x Obs
  matrix[m_obs, m_pred]  Sigma_12 = C_full[1:m_obs, (m_obs + 1):(m_obs + m_pred)];               // Obs x Pred
  matrix[m_pred, m_pred] Sigma_22 = C_full[(m_obs + 1):(m_obs + m_pred), (m_obs + 1):(m_obs + m_pred)];         // Pred x Pred
      
  matrix[m_obs, m_obs] L_Sigma11 = cholesky_decompose(Sigma_11);
  
  matrix[m_obs, m_pred] K_div = mdivide_left_tri_low(L_Sigma11, Sigma_12);
  
  // Conditional Covariance (Correlation only)
  matrix[m_pred, m_pred] Sigma_cond = Sigma_22 - K_div' * K_div;
      
  // Ensure symmetry 
      for(r in 1:m_pred) {
        for(col_idx in (r+1):m_pred) Sigma_cond[r, col_idx] = Sigma_cond[col_idx, r];
        Sigma_cond[r, r] += 1e-8; 
      }
      
  // Sampling of predicted coefficients
  array[nl] int idxs = id_l[l, 1:nl];
      
  for (k in 1:nl) {
  int current_col = idxs[k];
  vector[m_obs] theta_obs_k = col(theta, current_col);
  
  // Conditional mean
  vector[m_obs] alpha_vec = mdivide_left_tri_low(L_Sigma11, theta_obs_k);
  vector[m_pred] mu_cond = K_div' * alpha_vec;
  
  matrix[m_pred, m_pred] Cov_cond_scaled = square(coef_scale[current_col]) * Sigma_cond;
  theta_pred[, current_col] = multi_normal_rng(mu_cond, Cov_cond_scaled);
   }
  }
 }
}
"

Model=write_stan_file(hierarchical_wavelet_model)
mod=cmdstan_model(Model)

# 4) Group data and hyperparameters into a list for Stan
data_list <- list(
  m_obs                 = m_obs,
  m_pred                = m_pred,
  n                     = n,
  L                     = L,
  shrunk_levels         = shrunk_levels,
  coarse_levels         = coarse_levels,
  n_shrunk              = n_shrunk,
  n_l                   = n_l,
  max_nl                = max_nl,
  id_l                  = id_l,
  locations_obs         = locations_obs,
  locations_full        = locations_full,
  D_max                 = D_max,
  d_empirical           = d_empirical,
  p0                    = p0,
  slab_scale            = slab_scale,
  slab_df               = slab_df
)

fit_model = mod$sample(
  data            = data_list,
  seed            = 123,
  chains          = 2,
  parallel_chains = 2,
  iter_warmup     = 2000,
  iter_sampling   = 2000,
  adapt_delta     = 0.99,
  max_treedepth   = 12
)
