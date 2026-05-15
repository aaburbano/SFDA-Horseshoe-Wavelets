# Load required libraries
rm(list=ls(all=TRUE))
library("MASS")
library("ggplot2")
library("ggrepel")
library("ggpubr")
library("dplyr")
library("tidyr")
library("fda") 

load(".../Main_Data/coords.RData")

# ------------------------------------------------------------------------------
# 1. Grid Definition and Parameters
# ------------------------------------------------------------------------------
m = 20        # Number of spatial locations
n = 128       # Time points
t = 1:n

# Calculate the Euclidean distance matrix
dist_mat = as.matrix(dist(coords[, c("longitude", "latitude")])); max(dist_mat); min(dist_mat)

# Spatial covariance for the scalar data
sigma2_s = 0.5
phi_s = 0.4   # Moderate range

# Covariance Matrix (Exponential Function)
Sigma = sigma2_s * exp(-dist_mat / phi_s)

# ------------------------------------------------------------------------------
# 2. Construction of the Deterministic Trimodal Trend
# ------------------------------------------------------------------------------

# Auxiliary function for Gaussian peaks (Trimodal)
# Peaks at: t=26 , t=64, t=109
get_trimodal_trend <- function(t) {
  p1 = 1.8 * exp(-0.008 * (t - 26)^2)  
  p2 = 0.8 * exp(-0.003 * (t - 64)^2)
  p3 = 1.0 * exp(-0.008 * (t - 109)^2)
  return(p1 + p2 + p3)
}

# Logarithmic base (Functional intercept)
mu_base = 3.1
mu_t = mu_base + get_trimodal_trend(t)

# ------------------------------------------------------------------------------
# 3. Curve Generation via Spatial Basis Expansion
# ------------------------------------------------------------------------------
nbasis = 8 
norder = 4  
range_val = range(t)
basis_obj = create.bspline.basis(rangeval = range_val, nbasis = nbasis, norder = norder)

# Evaluate the basis on the temporal grid to obtain the design matrix (Phi)
# Dimension: n (time) x nbasis (K)
Phi = eval.basis(t, basis_obj)

# Coefficient Matrix: rows = locations (m), columns = basis functions (K)
# Each column k is a spatially correlated random field.

set.seed(123)

Coefs = matrix(NA, nrow = m, ncol = nbasis)

for(k in 1:nbasis) {
  Coefs[, k] = mvrnorm(1, mu = rep(0, m), Sigma = Sigma)
}

# ------------------------------------------------------------------------------
# 4. GENERATION OF OBSERVED SAMPLES (Yobs)
# ------------------------------------------------------------------------------
Z_latente_fixo = Phi %*% t(Coefs)

# Adding the mean trend to all curves
for(i in 1:m){ 
  Z_latente_fixo[, i] = Z_latente_fixo[, i] + mu_t
}

Ytrue = Z_latente_fixo
sigma_noise = 0.6  # sigma2 = 0.36

int = 300
Yobs = c()
Yobs_list = vector("list", int) 
for (l in 1:int) {
  # Variable noise for each iteration
  Noise = matrix(rnorm(n * m, mean = 0, sd = sigma_noise), nrow = n, ncol = m)
  
  # Add noise 
  Y_matrix = Z_latente_fixo + Noise
  
  # Prepare the data.frame
  df_temp = as.data.frame(Y_matrix)
  df_temp$id = l # Simulation identifier
  
  # Store in the list
  Yobs_list[[l]] = df_temp
}

# Consolidate everything at the end
Yobs = do.call(rbind, Yobs_list)

colnames(Yobs) = c(paste0("S", 1:m), "id")

# Calculate the temporal variance of each curve
var_sinal_curvas = apply(Ytrue, 2, var) 
snr_por_curva = var_sinal_curvas / (sigma_noise^2)
SNR_medio = mean(snr_por_curva)

#-----
any(Yobs < 0, na.rm = TRUE)

dados_apenas = Yobs[, 1:20] 

# Show the global Minimum and Maximum
range(dados_apenas, na.rm = TRUE)

# If the minimum is > 0, you are safe.

# ------------------------------------------------------------------------------
# Analysis
# ------------------------------------------------------------------------------

Sample = Yobs %>%
         filter(id == 7) %>% 
         dplyr::select(S1:S20)

par(mfrow = c(2, 1))
matplot(t, Ytrue, type = 'l', lty = 1)
matplot(t, Sample, type = 'l', lty = 1)

par(mfrow = c(1, 2))
matplot(t, Ytrue, type = 'l', lty = 1)
matplot(t, Sample, type = 'l', lty = 1)

# --- Find pairs of locations ---

dist_long = as.data.frame(as.table(dist_mat)) %>%
            rename(id1 = Var1, id2 = Var2, distance = Freq) %>%
            filter(as.numeric(id1) < as.numeric(id2)) 

# Find the closest pair
pair_nearby = dist_long %>% 
              arrange(distance) %>% 
              slice(1)

id_nearby1  = as.numeric(as.character(pair_nearby$id1))
id_nearby2  = as.numeric(as.character(pair_nearby$id2))

# Find the farthest pair
pair_distant = dist_long %>% 
               arrange(desc(distance)) %>% 
               slice(1)

id_distant1  = as.numeric(as.character(pair_distant$id1))
id_distant2  = as.numeric(as.character(pair_distant$id2))

par(mfrow = c(1, 2))
matplot(t, Ytrue[,c(7, 15)], type = 'l', lty = 1)
matplot(t, Ytrue[,c(6, 14)], type = 'l', lty = 1)

par(mfrow = c(1, 2))
matplot(t, Ytrue[,c(2, 16, 18)], type = 'l', lty = 1)
matplot(t, Ytrue[,c(4, 14, 17)], type = 'l', lty = 1)

par(mfrow = c(1, 2))
matplot(t, Ytrue[,c(6)], type = 'l', lty = 1)
matplot(t, Ytrue[,c(14,20)], type = 'l', lty = 1)

Tab = coords%>%
      mutate(name = rep(1:20))

Sites = ggplot(Tab, aes(x = longitude, y = latitude)) +
        geom_point(size = 3) +
        xlim(0,1) + ylim(0,1) +
        geom_text_repel(aes(label = name), size = 7, max.overlaps = Inf) +
        labs( x = "Coordinate 1",
              y = "Coordinate 2") +
        theme_light() +
        theme( axis.title = element_text(size = 30),
               axis.text = element_text(size = 20),
               plot.subtitle = element_text(size = 25))


# --- Visualization of the generated smooth curves ---

Ytrue_df = as.data.frame(Ytrue)
t_df = as.data.frame(t)
data_plot = bind_cols(t_df, Ytrue_df)
tab_plot = data_plot %>%
           pivot_longer(!t, names_to = "Local", values_to = "Val")

true = ggplot(tab_plot, aes(x = t, y = Val, group = factor(Local), color = factor(Local))) +
       geom_line() +
       labs( x = "Domain",
             y = "Value") +
       theme_light()+
       scale_x_continuous(expand = expansion(add = 2),breaks=seq(0,128,by=16),limits=c(0, 128))+
       scale_y_continuous(expand = expansion(add = 0.2),breaks=seq(1,6,by=1),limits=c(1, 6))+
       theme(legend.position ="none",axis.title=element_text(size = 27),axis.text = element_text(size = 20))

x11()
ggarrange(Sites, true, 
          labels = c("a )","b )"),
          font.label = list(size = 20),
          nrow = 1, ncol = 2,
          label.x = 0.0001,    # move horizontally (further right if increased)
          label.y = 0.95)      # move vertically (higher up if increased)
