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

setwd("C:/Users/alexb/OneDrive/Artigo2/Dados_Sinteticos/Simulation_Horseshoe/Ap_PM10_2/Data_PM10")
load("PM10_2025_day.RData")

Sample  = PM10_2025_day %>% 
          dplyr::select(-Day)

Sample = apply(Sample, c(1, 2), function(x) ifelse(is.na(x), NA, log(x)))
Sample = as.data.frame(Sample)

log_sample = Sample %>%
             dplyr::select(-ACO, -IZT, -TAH)

log_pred_sample_day = Sample %>%
                      dplyr::select(ACO, IZT, TAH) %>%
                      mutate(Day =1:128) %>%
                      relocate(Day,.before = ACO)

# Global Standardization properties (to rebuild scale)
matriz_Y = as.matrix(log_sample)
mu_Y     = mean(matriz_Y)
sd_Y     = sd(matriz_Y)

n_original = 128 
pad_size = 64

setwd("C:/Users/alexb/OneDrive/Artigo2/Dados_Sinteticos/Simulation_Horseshoe/Ap_PM10_2/Prediction/teste_2/mcmc/Resul")
load("diagnostic_pred.RData")
n = 256
m = 3

ii=1:n
jj=1:m
combinations = expand.grid(jj, ii);dim(combinations)

tab1 =c(paste0("theta_pred[", combinations$Var1, ",", combinations$Var2,"]"))

coeff = diagnostic_pred %>%
        dplyr::filter(variable %in% tab1) %>%
        dplyr::select(mean)%>%
        as.matrix() %>%
        c()

theta_hat = matrix(coeff, nrow=m, ncol=n, byrow = FALSE)

# ------------------------------------------------------------------------------
# RECONSTRUCTION OF THE THETA MATRIX
# ------------------------------------------------------------------------------
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

# Remove the mirrored edges added at the beginning
inicio_real = pad_size + 1
fim_real = pad_size + n_original

Y_hat_padronizado_recortado = Y_hat_padronizado_padded[inicio_real:fim_real, ]

# De-standardize to return to the original physical scale using only the core (cropped) region
Y_hat_pm10 = (Y_hat_padronizado_recortado * sd_Y) + mu_Y

# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------
t_original = 1:n_original # Ensures the days axis has the size of the original curve (e.g., 128)

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

save(log_pred_sample_day, file = ".../log_pred_sample_day.RData")

#-------------------------------------------------------------------------------
# Prediction Data Preparation
#-------------------------------------------------------------------------------
rm(list=ls(all=TRUE))

library("MASS")    
library("fields")  
library("ggplot2")
library("dplyr")
library("tidyr")

setwd("C:/Users/alexb/OneDrive/Artigo2/Dados_Sinteticos/Simulation_Horseshoe/Ap_PM10_2/Data_PM10")
load("Data_PM10_2025.RData")
load("PM10_2025_day.RData")

setwd("C:/Users/alexb/OneDrive/Artigo2/Dados_Sinteticos/Simulation_Horseshoe/Ap_PM10_2/Prediction/teste_2/mcmc")
load("log_pred_sample_day.RData")

data = Data_PM10_2025 %>%
       dplyr::select(FECHA, HORA, ACO ,IZT, TAH)

num_miss = data %>%
           dplyr::select(ACO,IZT,TAH)%>%
           summarise_all(funs(sum(is.na(.))))

tb1 = data %>%
      dplyr::select(-FECHA, - HORA) %>%
      mutate(Day = rep(1:128, each= 24)) %>%
      relocate(Day,.before = FAC) %>%
      group_by(Day)%>%
      summarise(across(FAC:UIZ, ~ median(.x, na.rm = TRUE)))

tb2 = tb1 %>%
      dplyr::select(-Day) 

tb3 = apply(tb2, c(1, 2), function(x) ifelse(is.na(x), NA, log(x)))

tb4  = PM10_2025_day %>% 
       dplyr::select(-Day)

tb5 = as.data.frame(apply(tb4, c(1, 2), function(x) ifelse(is.na(x), NA, log(x))))

tb6 = tb5 %>%
      dplyr::select(ACO)

Pred_logPM10_2025_day = bind_cols(tb1[,1], tb3, tb6)

save(Pred_logPM10_2025_day,file = "C:/Users/alexb/OneDrive/Artigo2/Dados_Sinteticos/Simulation_Horseshoe/Ap_PM10_2/Data_PM10/Pred_logPM10_2025_day.RData")
