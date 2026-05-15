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

load(".../Main_Data/coords.RData")
load(".../Main_Data/SNR/.../Yobs.RData")

#-------------------------------------------------------------------------------
# General parameters
#-------------------------------------------------------------------------------
m = 20
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
id_l = array(0, dim = c(L, max_nl))
for (l in 1:L) {
  cols = which(level_id == (l - 1))
  id_l[l, 1:length(cols)] = cols
}

dist_mat = as.matrix(dist(coords[, c("longitude", "latitude")]))
D_max = max(dist_mat)
slab_scale = 2.0
slab_df = 4.0

ii=1:n
jj=1:m
combinations = expand.grid(jj, ii); dim(combinations)
tab = c(paste0("theta[", combinations$Var1, ",", combinations$Var2,"]"))

#-------------------------------------------------------------------------------
###############################
### Monte Carlo Simulation ####
###############################
hierarchical_wavelet_model="
data {
  int<lower=1> m;                                                         
  int<lower=1> n;                                                         
  int<lower=1> L;                                                         
  int<lower=1> shrunk_levels;                                             
  int<lower=1> coarse_levels;                                             
  array[L] int n_l;                                                       
  int<lower=1> max_nl;                                                    
  array[L,max_nl] int id_l;                                               
  matrix[m, n] d_empirical;                                               
  matrix[m, 2] coords;                                                    
  real<lower=0> D_max;                                                    
  int<lower=1> n_shrunk;                                                  
  vector<lower=0>[shrunk_levels] p0;                                      
  real<lower=0> slab_scale;                                               
  real<lower=0> slab_df;                                                  
}

transformed data {
  array[m] vector[2] s_locs;
  for (i in 1:m) {
    s_locs[i] = to_vector(coords[i]);
  }
}

parameters {
  real<lower=0> sigma_obs;                                                
  matrix[m, n] z;                                                         
  vector<lower=0>[coarse_levels] sigma_coarse;                            
  vector<lower=0>[L] phi;                                                 
  vector<lower=0>[shrunk_levels] aux1_global;
  vector<lower=0>[shrunk_levels] aux2_global;
  vector<lower=0>[n_shrunk] aux1_local;
  vector<lower=0>[n_shrunk] aux2_local;
  vector<lower=0>[shrunk_levels] caux;
}

transformed parameters {
  matrix[m, n] theta;                                                     
  vector<lower=0>[shrunk_levels] tau_global;
  vector<lower=0>[shrunk_levels] c;
  
  {
    int local_idx = 1;                                                      
    
    for (l in 1:L) {
      int nl = n_l[l];
      if (nl > 0) {
        matrix[m, m] R = gp_matern32_cov(s_locs, 1, phi[l]);                        
        for (i in 1:m) {R[i, i] += 1e-8;}                                           
        matrix[m, m] L_cov = cholesky_decompose(R);
        
        array[nl] int idxs = id_l[l, 1:nl];                                         
        
        if (l <= coarse_levels) {                                                   
          for (k in 1:nl) {
            int col_idx = idxs[k];
            theta[, col_idx] = sigma_coarse[l] * (L_cov * z[, col_idx]);            
          }
        }
        else {                                                                      
          int sl = l - coarse_levels;                                               
          real tau0 = (p0[sl] / (nl - p0[sl])) * (sigma_obs / sqrt(m * 1.0));
          tau_global[sl] = aux1_global[sl] * sqrt(aux2_global[sl]) * tau0;        
          c[sl] = slab_scale * sqrt(caux[sl]);                                    
          
          for (k in 1:nl) {
            int col_idx = idxs[k];
            real lambda_k = aux1_local[local_idx] * sqrt(aux2_local[local_idx]);    
            
            real lambda_tilde_k = sqrt( (c[sl]^2 * square(lambda_k)) / (c[sl]^2 + tau_global[sl]^2 * square(lambda_k)) );
            
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
"
Model=write_stan_file(hierarchical_wavelet_model)
mod=cmdstan_model(Model)

#-------------------------------------------------------------------------------
# Note: Global counters don't work natively inside parallel workers. 
# Removed the global counter and updated the log message strategy.
counter=0
MonteCarlo = function(step){
  
  log_file = "/home/alextk4/Simu/MonteCarlo_log.txt"  # Defines the single log file name
  
  log_msg <- function(msg) {
    cat(msg, file = log_file, append = TRUE, sep = "\n")  # Writes the message to the file
    flush.console()  # Flushes the console and file immediately
  }
  
  #----
  Y = Yobs %>%
      filter(id == step) %>%
      dplyr::select(dplyr::starts_with("S"))%>%
      as.matrix()
  
  mu_rep = mean(Y)
  sd_rep = sd(Y)
  Y_pad  = (Y - mu_rep) / sd_rep
  
  Y_padded = matrix(NA, nrow = n, ncol = m)
  for(j in 1:m) {
    curve = Y_pad[, j]
    Y_padded[, j] = c(rev(curve[1:pad_size]), 
                      curve, 
                      rev(curve[(pad_size + 1):n_original]))
  }
  
  # DWT with DaubExPhase
  wavelet_list = vector("list", m)
  for(j in 1:m) {
    wd_obj = wavethresh::wd(Y_padded[, j], 
                            filter.number = 4, 
                            family = "DaubExPhase",
                            bc="periodic") 
    wavelet_list[[j]] = wd_obj
  }
  
  # ----------------------------------------------------------------------------
  # Function to extract and organize coefficients into the required format 
  # ----------------------------------------------------------------------------
  coef_empirical = function(wavelet_list, n, m, L) {
    d_empirical = matrix(NA, nrow = n, ncol = m)
    
    for (j in seq_len(m)) {
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
  
  d_empirical = t(coef_empirical(wavelet_list = wavelet_list, n = n, m = m, L = L))
  
  data_list <- list(
    m                     = m,
    n                     = n,
    L                     = L,
    shrunk_levels         = shrunk_levels,
    coarse_levels         = coarse_levels,
    n_shrunk              = n_shrunk,
    n_l                   = n_l,
    max_nl                = max_nl,
    id_l                  = id_l,
    coords                = coords,
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
    parallel_chains = 2, # NOTE: Be mindful of CPU allocation when using sfLapply + parallel_chains
    iter_warmup     = 2000,
    iter_sampling   = 2000,
    adapt_delta     = 0.99,
    max_treedepth   = 12
  )
  
  ##################
  ### Parameters ###
  ##################
  
  # 1
  draws      = fit_model$draws(c("phi", "sigma_obs", "sigma_coarse"))
  variable   = as_draws_matrix(draws)
  diagnostic = summarise_draws(variable, "mean","median", "sd", default_convergence_measures())
  
  # 2
  draws2      = fit_model$draws("theta")
  variable2   = as_draws_matrix(draws2)
  diagnostic2 = summarise_draws(variable2, "mean","median", "sd", default_convergence_measures())
  
  # Prepares the data.frame
  
  # 1
  diag = as.data.frame(diagnostic)
  diag$id = step # Simulation identifier
  
  # 2
  diag2 = as.data.frame(diagnostic2)
  diag2$id = step # Simulation identifier
  
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
  
  #-----------------------------------------------------------------------------
  # Curve Reconstruction
  #-----------------------------------------------------------------------------
  coeff = diag2 %>%
          filter(variable %in% tab) %>%
          dplyr::select(mean)%>%
          as.matrix() %>%
          c()
  
  theta_hat = matrix(coeff, nrow=m, ncol=n, byrow = FALSE)
  
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
  
  hat_curves = matrix(0, nrow=m, ncol=n)
  
  for (i in 1:m) {
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
  
  # 3
  Y_hat$id = step # Simulation identifier
  
  write.table(Y_hat,file="/home/alextk4/Simu/Y_hat.txt",
              sep = ",",
              col.names=FALSE,
              row.names=FALSE ,
              append = TRUE)
              
  #-----------------------------------------------------------------------------
  counter <<- counter + 1
  # Logs the progress in the log file
  log_msg(paste("Finishing iteration:", step, "\n"))
}

# Define number of CPUs and total simulations beforehand
num_cpus = 6 # Example, adjust as necessary
N_sim = 300  

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
