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

#-------------------------------------------------------------------------------
Sample = PM10_2025_day %>% 
         dplyr::select(-Day)

Sample = apply(Sample, c(1, 2), function(x) ifelse(is.na(x), NA, log(x)))


# (1) Global Standardization

matriz_Y = as.matrix(Sample)
mu_Y     = mean(matriz_Y)
sd_Y     = sd(matriz_Y)

# Standardize the entire matrix while keeping the spatial structure intact
Sample_padronizado = (matriz_Y - mu_Y) / sd_Y
Sample_padronizado = as.data.frame(Sample_padronizado)

# ==============================================================================
# START OF MODIFICATION 1: DYNAMIC PADDING (MIRRORING)
# ==============================================================================
n_original = dim(Sample_padronizado)[1] 
m = dim(Sample_padronizado)[2]           

# Calculate the amount to add on each side (half of the signal)
pad_size = n_original / 2 

# The new total size will be exactly double the original (e.g., 128 -> 256)
n_padded = n_original * 2 
Y_padded = matrix(NA, nrow = n_padded, ncol = m)

for(j in 1:m) {
  curva = Sample_padronizado[, j]
  # Perfect mirroring to handle periodic boundary conditions for the algorithm
  Y_padded[, j] = c(rev(curva[1:pad_size]), curva, rev(curva[(pad_size + 1):n_original]))
}

Sample_padronizado = as.data.frame(Y_padded)
# ==============================================================================
# END OF MODIFICATION 1
# ==============================================================================

# 'n' now takes the value of the inflated size (n_padded) for Stan processing
n = dim(Sample_padronizado)[1]
L = log2(n)                                                                 # Maximum resolution level
coarse_levels = floor(log2(log(n))) + 1                                     # Coarse levels
shrunk_levels = L - coarse_levels                                           # Number of levels that will undergo shrinkage
n_l = c(2, 2^(1:(L - 1))); sum(n_l)                                         # Number of coefficients per level (n_l) 

n_shrunk = sum(n_l[(coarse_levels + 1):L])                                   # Total number of coefficients undergoing shrinkage
p0 = round(sapply((coarse_levels + 1):L, function(l) max(1, 0.1 * n_l[l])), 0) # Expected number of non-zero coefficients per level (p0) 
level_id = rep(0:(L - 1), times = n_l); length(level_id)                     # Map each of the 'n' positions in d_matrix to its level 

# Pre-compute id_l in R
max_nl = max(n_l)
id_l   = array(0, dim = c(L, max_nl))

for (l in 1:L) {
  cols = which(level_id == (l - 1))
  id_l[l, 1:length(cols)] = cols
}

# Calculate the Euclidean distance matrix
locations = coord %>%
            dplyr::select(Lon, Lat)

dist_mat = as.matrix(dist(locations[, c("Lon", "Lat")]))
D_max = max(dist_mat)
slab_scale = 2.0
slab_df = 4.0

# ------------------------------------------------------------------------------
# (2) Discrete Wavelet Transform using Daubechies
# ------------------------------------------------------------------------------
wavelet_list = vector("list", m)
for(j in 1:m) {
  wd_obj = wavethresh::wd(Sample_padronizado[, j], 
                          filter.number = 4, 
                          family = "DaubExPhase",
                          bc = "periodic") # Keep "periodic" to maintain the 2^j tree structure
  wavelet_list[[j]] = wd_obj
}

# ------------------------------------------------------------------------------
# (3) Function to extract and organize coefficients format
# ------------------------------------------------------------------------------
coef_empirical = function(wavelet_list, n, L) {
  m = length(wavelet_list)
  d_empirical = matrix(NA, nrow = n, ncol = m)
  
  for (j in seq_len(m)) {
    w = wavelet_list[[j]]
    c0 = accessC(w, level = 0)     
    d0 = accessD(w, level = 0)      
    vetor_organizado = c(c0, d0)
    
    for (i in 1:(L - 1)) {
      di = accessD(w, level = i)    
      vetor_organizado = c(vetor_organizado, di)
    }
    
    if (length(vetor_organizado) != n) {
      stop(sprintf("Incorrect length at j = %d: obtained %d, expected %d", j, length(vetor_organizado), n))
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
  int<lower=1> m;                                                              
  int<lower=1> n;                                                              
  int<lower=1> L;                                                              
  int<lower=1> shrunk_levels;                                                  
  int<lower=1> coarse_levels;                                                  
  array[L] int n_l;                                                            
  int<lower=1> max_nl;                                                         
  array[L,max_nl] int id_l;                                                    
  matrix[m, n] d_empirical;                                                    
  matrix[m, 2] locations;                                                         
  real<lower=0> D_max;                                                         
  int<lower=1> n_shrunk;                                                       
  vector<lower=0>[shrunk_levels] p0;                                           
  real<lower=0> slab_scale;                                                    
  real<lower=0> slab_df;                                                       
}

transformed data {
  array[m] vector[2] s_locs;
  for (i in 1:m) {
    s_locs[i] = to_vector(locations[i]);
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

# 4) Group data and hyperparameters into a list for Stan
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
  locations             = locations,
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

#-------------------------------------------------------------------------------
fit_model$diagnostic_summary()

sumario = fit_model$summary()

pior_rhat = sumario %>% 
  arrange(desc(rhat)) %>% 
  head(15)
print(pior_rhat[, c("variable", "rhat")])

pior_ess = sumario %>% 
  arrange(ess_bulk) %>% 
  head(15)
print(pior_ess[, c("variable", "ess_bulk")])
#-------------------------------------------------------------------------------

draws_1 = fit_model$draws(c("theta", "phi", "sigma_obs", "sigma_coarse"))
variable = as_draws_matrix(draws_1)
diagnostic_pm10 = summarise_draws(variable, "mean","median", "sd", default_convergence_measures())

ii=1:n
jj=1:m
combinations = expand.grid(jj, ii);dim(combinations)

tab1 =c(paste0("theta[", combinations$Var1, ",", combinations$Var2,"]"))

coeff = diagnostic_pm10 %>%
  filter(variable %in% tab1) %>%
  dplyr::select(mean)%>%
  as.matrix() %>%
  c()

theta_hat = matrix(coeff, nrow=m, ncol=n, byrow = FALSE)

# ------------------------------------------------------------------------------
# THETA MATRIX RECONSTRUCTION
# ------------------------------------------------------------------------------
reorder_theta_hat = function(theta_hat) {
  if (!is.matrix(theta_hat)) stop("'theta_hat' must be a matrix (m x n).")
  m = nrow(theta_hat)
  n = ncol(theta_hat)
  L = log2(n)
  if (abs(L - round(L)) > .Machine$double.eps^0.5) stop("n is not a power of two.")
  L = as.integer(L)
  
  theta_hat_reordered = matrix(NA_real_, nrow = m, ncol = n)
  
  for (j in seq_len(m)) {
    original = theta_hat[j, ]          
    novo     = numeric(length = n)     
    novo[1] = original[1]
    idx <- 2  
    
    for (i in seq(from = L - 1, to = 0, by = -1)) {
      inicio_i = 1 + 2^i
      final_i  = 2^(i + 1)
      comprimento_i = 2^i  
      
      novo[idx:(idx + comprimento_i - 1)] = original[inicio_i:final_i]
      idx = idx + comprimento_i
    }
    theta_hat_reordered[j, ] = novo
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

Y_hat_padronizado_padded = as.data.frame(t(hat_curves)) # 256x20 Matrix (or 2N x m)

# ==============================================================================
# START OF MODIFICATION 3: TRIMMING AND RETURN TO ORIGINAL SIGNAL
# ==============================================================================
# Remove the boundaries mirrored at the beginning to retain only the clean original signal
inicio_real = pad_size + 1
fim_real = pad_size + n_original

Y_hat_padronizado_recortado = Y_hat_padronizado_padded[inicio_real:fim_real, ]

# Destandardize to return to the original physical scale using only the core signal
Y_hat_pm10 = (Y_hat_padronizado_recortado * sd_Y) + mu_Y
# ==============================================================================
# END OF MODIFICATION 3
# ==============================================================================

# Plotting
t_original = 1:n_original # Ensure the day axis matches the original curve size (e.g., 128)

mat_smooth = Y_hat_pm10 %>%
  mutate(Day = t_original) %>%
  relocate(Day, .before = V1)

Curves = mat_smooth %>%
  pivot_longer(!Day, names_to = "Local", values_to = "Val")

x11()
smooth = ggplot(Curves, aes(x=Day, y=Val,group=factor(Local)))+geom_line()+
  labs(x = "Day",
       y = "Value") +
  theme_light()+
  scale_x_continuous(expand = expansion(add = 2),breaks=seq(0, n_original,by=16),limits=c(0, n_original))+
  scale_y_continuous(expand = expansion(add = 0.2),breaks=seq(1,6,by=1),limits=c(1, 6))+
  theme(legend.position ="none", 
        axis.title=element_text(size = 30), 
        axis.text = element_text(size=15), 
        plot.subtitle = element_text(size=25))
smooth
