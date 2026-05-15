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

load(".../Simulation_Studies/Main_Data/SNR/Medium/Data_Simulation_3/locations_obs.RData")
load(".../Simulation_Studies/Main_Data/SNR/Medium/Data_Simulation_3/locations_full.RData")
load(".../Simulation_Studies/Main_Data/SNR/Medium/Data_Simulation_3/Yobs.RData")

#-------------------------------------------------------------------------------
# General parameters
#-------------------------------------------------------------------------------
m_obs  = 18
m_pred = 2
n_original = 128
n = n_original * 2
pad_size = n_original / 2
L = log2(n) # L will be 8 due to padding (256)
coarse_levels = floor(log2(log(n))) + 1
shrunk_levels = L - coarse_levels
n_l = c(2, 2^(1:(L - 1)))
max_nl = max(n_l)
n_shrunk = sum(n_l[(coarse_levels+ 1):L])
p0 = round(sapply((coarse_levels+1):L, function(l) max(1, 0.1 * n_l[l])),0)

level_id = rep(0:(L - 1), times = n_l)
id_l   = array(0, dim = c(L, max_nl))
for (l in 1:L) {
  cols = which(level_id == (l - 1))
  id_l[l, 1:length(cols)] = cols
}

dist_mat = as.matrix(dist(locations_full[, c("longitude", "latitude")]))
D_max = max(dist_mat)
slab_scale = 2.0
slab_df = 4.0

ii=1:n
jj=1:m_pred
combinations = expand.grid(jj, ii); dim(combinations)
tab = c(paste0("theta_pred[", combinations$Var1, ",", combinations$Var2,"]"))

#-------------------------------------------------------------------------------
###############################
### Monte Carlo Simulation ####
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
  matrix[(m_obs + m_pred), 2] locations_full;             // Observed + predicted coordinates matrix
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
  
// Partition the Matrix
  
  matrix[m_obs, m_obs]   Sigma_11 = C_full[1:m_obs, 1:m_obs];                                     // Obs x Obs
  matrix[m_obs, m_pred]  Sigma_12 = C_full[1:m_obs, (m_obs + 1):(m_obs + m_pred)];                // Obs x Pred
  matrix[m_pred, m_pred] Sigma_22 = C_full[(m_obs + 1):(m_obs + m_pred), (m_obs + 1):(m_obs + m_pred)]; // Pred x Pred
      
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

#-------------------------------------------------------------------------------
# Note: Global counters don't work natively inside parallel workers. 
# Removed the global counter and updated log message.

MonteCarlo=function(step){
  
  log_file = "/home/alextk4/Simu/MonteCarlo_log.txt"  # Defines the single log file name
  
  log_msg <- function(msg) {
    cat(msg, file = log_file, append = TRUE, sep = "\n")  # Writes the message to the file
    flush.console()  # Flushes the console and file immediately
  }
  #----
  Y = Yobs %>%
      filter(id == step) %>%
      dplyr::select(dplyr::starts_with("S") & !c("S9", "S19")) %>%
      as.matrix()
  
  mu_rep = mean(Y)
  sd_rep = sd(Y)
  Y_pad  = (Y - mu_rep) / sd_rep
  
  Y_padded = matrix(NA, nrow = n, ncol = m_obs)
  for(j in 1:m_obs) {
    curve = Y_pad[, j]
    Y_padded[, j] = c(rev(curve[1:pad_size]), 
                      curve, 
                      rev(curve[(pad_size + 1):n_original]))
  }
  
  # DWT with DaubExPhase
  wavelet_list = vector("list", m_obs)
  for(j in 1:m_obs) {
    wd_obj = wavethresh::wd(Y_padded[, j], 
                            filter.number = 4, 
                            family = "DaubExPhase",
                            bc="periodic") 
    wavelet_list[[j]] = wd_obj
  }
  
  # ----------------------------------------------------------------------------
  # Function to extract and organize coefficients into the required format 
  # ----------------------------------------------------------------------------
  
  coef_empirical = function(wavelet_list, n, m_obs, L) {
    d_empirical = matrix(NA, nrow = n, ncol = m_obs)
    
    for (j in seq_len(m_obs)) {
      w = wavelet_list[[j]]
      c0 = accessC(w, level = 0)     
      d0 = accessD(w, level = 0)      
      organized_vector = c(c0, d0)
      
      for (i in 1:(L - 1)) {
        di = accessD(w, level = i)    
        organized_vector = c(organized_vector, di)
      }
      
      if (length(organized_vector) != n) {
        stop(sprintf("Incorrect length at j = %d: obtained %d, expected %d", j, length(organized_vector), n))
      }
      d_empirical[, j] = organized_vector
    }
    return(d_empirical)
  }
  
  d_empirical = t(coef_empirical(wavelet_list = wavelet_list, n = n, m_obs = m_obs, L = L))
  
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
  
  ##################
  ### Parameters ###
  ##################
  #1
  draws      = fit_model$draws(c("theta", "phi", "sigma_obs", "sigma_coarse"))
  variable   = as_draws_matrix(draws)
  diagnostic = summarise_draws(variable, "mean","median", "sd", default_convergence_measures())
  
  #2
  draws2      = fit_model$draws("theta_pred")
  variable2   = as_draws_matrix(draws2)
  diagnostic2 = summarise_draws(variable2, "mean","median", "sd", default_convergence_measures())

  
  # Prepares the data.frame
  
  #1
  diag = as.data.frame(diagnostic)
  diag$id = step # Simulation identifier
  
  #2
  diag2 = as.data.frame(diagnostic2)
  diag2$id = step # Simulation identifier
  
  #3
  diag3 = as.data.frame(variable2)
  diag3$id = step # Simulation identifier
  
  #---
  write.table(diag,file="/home/alextk4/Simu/diag.txt",
              sep = ",",
              col.names=FALSE,
              row.names=FALSE ,
              append = TRUE)
  
  write.table(diag2,file="/home/alextk4/Simu/diag2.txt",
              sep = ",",
              col.names=FALSE,
              row.names=FALSE ,
              append = TRUE)
  
  write.table(diag3,file="/home/alextk4/Simu/diag3.txt",
              sep = ",",
              col.names=FALSE,
              row.names=FALSE ,
              append = TRUE)
  #-----------------------------------------------------------------------------
  
  # Curve Reconstruction
  
  coeff = diag2 %>%
          filter(variable %in% tab) %>%
          dplyr::select(mean)%>%
          as.matrix() %>%
          c()
  
  theta_hat = matrix(coeff, nrow=m_pred, ncol=n, byrow = FALSE)
  
  reorder_theta_hat = function(theta_hat) {
    if (!is.matrix(theta_hat)) stop("The 'theta_hat' object must be an (m x n) matrix.")
    m = nrow(theta_hat)
    n = ncol(theta_hat)
    L = log2(n)
    if (abs(L - round(L)) > .Machine$double.eps^0.5) stop("n is not a power of two.")
    L = as.integer(L)
    
    theta_hat_reordered = matrix(NA_real_, nrow = m, ncol = n)
    
    for (j in seq_len(m)) {
      original = theta_hat[j, ]          
      new_vector = numeric(length = n)     
      new_vector[1] = original[1]
      idx <- 2  
      
      for (i in seq(from = L - 1, to = 0, by = -1)) {
        start_i  = 1 + 2^i
        end_i    = 2^(i + 1)
        length_i = 2^i  
        
        new_vector[idx:(idx + length_i - 1)] = original[start_i:end_i]
        idx = idx + length_i
      }
      theta_hat_reordered[j, ] = new_vector
    }
    return(theta_hat_reordered)
  }
  
  theta_hat_reordered = reorder_theta_hat(theta_hat)
  
  template = wd(rep(0, n), 
                filter.number = 4, 
                family = "DaubExPhase", 
                bc = "periodic")
  
  hat_curves = matrix(0, nrow=m_pred, ncol=n)
  
  for (i in 1:m_pred) {
    wd_obj = putC(template, level = 0, v = theta_hat_reordered[i,1])
    wd_obj$D = theta_hat_reordered[i,2:n ]
    hat_curves[i, ] = wr(wd_obj)
  }
  
  Yhat_standardized = as.data.frame(t(hat_curves))
  
  # Return to original signal
  # Remove the edges mirrored at the beginning
  real_start = pad_size + 1
  real_end = pad_size + n_original
  
  Y_hat_standardized = Yhat_standardized[real_start:real_end, ]
  
  # FIXED: Unstandardize to return to the original physical scale using only the core (inner part)
  # Uses mu_rep and sd_rep defined locally at the start of the function
  Y_hat = (Y_hat_standardized * sd_rep) + mu_rep
  
  #3
  
  Y_hat$id = step # Simulation identifier
  
  write.table(Y_hat,file="/home/alextk4/Simu/Y_hat.txt",
              sep = ",",
              col.names=FALSE,
              row.names=FALSE ,
              append = TRUE)
  #-----------------------------------------------------------------------------
  
  # Logs the progress in the log file
  log_msg(paste("Finishing iteration:", step, "\n"))
}

# Define number of CPUs and total simulations beforehand
num_cpus = 6 # Example, adjust as necessary
N_sim = 300  # Replaced 'int' with a proper variable name

sfInit(parallel=TRUE, cpus=num_cpus)
sfLibrary(dplyr)
sfLibrary(cmdstanr)
sfLibrary(posterior)
sfLibrary(bayesplot)
sfLibrary(coda)
sfLibrary(wavethresh)
sfExportAll()

# Function that I want to compute multiple times using sfLapply:
sfLapply(1:N_sim, fun=MonteCarlo) 
sfStop()
