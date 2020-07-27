library(rstan)
library(ggplot2)
library(parallel)

# ALL files should be in the same working directory. The user should set this directory under the 'path' varialbe
options(max.print=999999)
path <- "" ###### USER INPUT HERE #####
setwd(path)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

####################################################################
###################             DATA             ###################
####################################################################


##########################################
# Function for creating parameter vector #
##########################################

create.params <- function(input){
  with( as.list(input),{
    
    # List with names of all possible compartments
    comps <- list("Lungs" = "Lungs" ,"Liver"="Liver", "Spleen"="Spleen", "Kidneys"="Kidneys", "Heart"="Heart",  "Brain"="Brain",
                  "Uterus"="Uterus", "Skeleton"="Skeleton", "Skin"="Skin", "Soft"="Soft") # List with names of all possible compartments
    
    ### Density of tissues/organs
    d_tissue <- 1 #g/ml
    d_skeleton <- 1.92 #g/ml
    d_adipose <- 0.940 #g/ml
    
    Q_total <- (1.54*weight^0.75)*60 # Total Cardiac Output (ml/h)
    
    Total_Blood <- 0.06*weight+0.77 # Total blood volume (ml)
    
    #Arterial blood volume
    Vart <- 0.15*Total_Blood #(ml)
    
    #Veins blood volume
    Vven <-0.64*Total_Blood #(ml)
    
    
    
    #Tissue weight fraction 
    Tissue_fractions <- c(0.5, 3.66, 0.2, 0.73, 0.33, 0.57, 0.011, 10, 19.03, NA)/100 # % of BW. Na values refers to the volume of the rest organs(RoB)
    
    #Regional blood flow fraction
    Regional_flow_fractions <- c(100, 17.4, 1.22, 14.1, 4.9, 2, 1.1, 12.2, 5.8, NA)/100 # % of total cardiac output
    
    #Capillary volume fractions (fractions of tissue volume)
    Capillary_fractions <- c(0.36, 0.21, 0.22, 0.16, 0.26, 0.03, 0.04, 0.04, 0.02, 0.04) # fraction of tissue volume
    # Where NA, it is the same value as Rest of Body due to luck od data for uterus and adipose
    
    #Macrophage content (# per gram of tissue)
    Macrophage_fraction <- c(0.04, 0.1, 0.3, 0.02, 0.02, 0.04, 0.04, 0.04, 0.02, 0.04) 
    Macrophage_fraction_blood <- 0.01
    #li <- c(2.69e06, 2.72e07, 2.08e08, 9.90e04, 7.60e04,3.06e05)
    
    CLE_ua <- 0
    CLE_muc <- 0.01181567
    
    W_tis <- rep(0,length(comps))
    V_tis <- rep(0,length(comps))
    V_cap <- rep(0,length(comps))
    W_macro <- rep(0,length(comps))  #one more for blood compartment
    Q <- rep(0,length(comps))
    
    
    for (i in 1:(length(comps)-1)) {
      control <- comps[i]
      
      Tissue_fractions[i] <- ifelse(is.na(control), NA, Tissue_fractions[i])
      Regional_flow_fractions[i] <- ifelse(is.na(control), NA, Regional_flow_fractions[i])
      Capillary_fractions[i] <- ifelse(is.na(control), NA, Capillary_fractions[i])
      Macrophage_fraction[i] <- ifelse(is.na(control), NA, Macrophage_fraction[i])
      
      ### Calculation of tissue weights  
      W_tis[i] <- weight*Tissue_fractions[i]
      W_macro[i] <- W_tis[i]*Macrophage_fraction[i]
      
      ###Calculation of tissue volumes
      
      if (i==9){
        V_tis[i] <- W_tis[i]/d_skeleton
      } else{
        V_tis[i] <- W_tis[i]/d_tissue 
      }
      
      ###Calculation of capillary volumes
      V_cap[i] <- V_tis[i]*Capillary_fractions[i]
      
      
      ###Calculation of regional blood flows
      Q[i] <- Q_total*Regional_flow_fractions[i]
    }
    
    ### Calculations for "Rest of Body" compartment
    W_tis[10] <- weight - sum(W_tis[1:(length(W_tis)-1)], na.rm = TRUE)
    V_tis[10] <- W_tis[10]/d_adipose     #(considering that the density of the rest tissues is 1 g/ml)
    Q[10] <- Q_total - sum(Q[2:(length(Q)-1)],na.rm = TRUE)
    V_cap[10] <- V_tis[10]*Capillary_fractions[10]
    # V_cap[9] <- Total_Blood - Vven - Vart - sum(V_cap[1:(length(V_cap)-1)], na.rm = TRUE) #this is problematic because it produces negative number
    W_macro[10] <- W_tis[10]*Macrophage_fraction[10]
    #Capillary_fractions[1] <- V_cap[1]/V_tis[1]
    
    Macro_blood <- Macrophage_fraction_blood*Total_Blood
    
    
    return(c("Vlu_tis"= V_tis[1], "Vli_tis"= V_tis[2], "Vspl_tis"= V_tis[3], "Vki_tis"= V_tis[4],"Vht_tis"= V_tis[5],
             "Vbr_tis"= V_tis[6], "Vut_tis"= V_tis[7], "Vskel_tis"= V_tis[8], "Vskin_tis"= V_tis[9], "Vsoft_tis"= V_tis[10],
             "Vven"=Vven, "Vart"=Vart, "V_blood"=Total_Blood,
             
             "Vlu_cap"= V_cap[1], "Vli_cap"= V_cap[2], "Vspl_cap"= V_cap[3], "Vki_cap"= V_cap[4], "Vht_cap"= V_cap[5],
             "Vbr_cap"= V_cap[6], "Vut_cap"= V_cap[7], "Vskel_cap"= V_cap[8], "Vskin_cap"= V_cap[9], "Vsoft_cap"= V_cap[10],
             
             "Q_lu"= Q[1], "Q_li"= Q[2], "Q_spl"= Q[3], "Q_ki"= Q[4], "Q_ht"= Q[5], "Q_br"= Q[6], "Q_ut"= Q[7], "Q_skel"= Q[8],
             "Q_skin"= Q[9], "Q_soft"= Q[10],  "Q_total"=Q_total, 
             
             "Wm_blood" = Macro_blood, "Wm_lu"= W_macro[1], "Wm_li"= W_macro[2], "Wm_spl"= W_macro[3], "Wm_ki"= W_macro[4], 
             "Wm_ht"= W_macro[5], "Wm_br"= W_macro[6], "Wm_ut"= W_macro[7], "Wm_skel"= W_macro[8], "Wm_skin"= W_macro[9],
             "Wm_soft"= W_macro[10],
             
             "Inhaled.vol.rate" = Inhaled.vol.rate, "dep.ua" = dep.ua, "dep.tb" = dep.tb, "dep.al" = dep.al, 
             "exposure.concentration" = exposure.concentration, "inh.duration" = exposure.time,
             
             "CLE_ua" = CLE_ua, "CLE_muc" = CLE_muc))
    
  }) 
}

#load data from xlsx files
biodist_data <- openxlsx::read.xlsx("biodist_data.xlsx", sheet=1,colNames = TRUE,rowNames = TRUE)
input.data <- openxlsx::read.xlsx("biodist_data.xlsx", sheet=3,colNames = TRUE,rowNames = TRUE)

feces <- openxlsx::read.xlsx("feces.xlsx", sheet=2,colNames = TRUE)[,2]*input.data[1,5]
urine<- openxlsx::read.xlsx("urine.xlsx", sheet=2,colNames = TRUE)[,2]*input.data[7,5]
excreta_time <- c( 3.5, 7, 10.5, 14, 17.5, 21, 24.5, 27)*24
biodist_time  <- c(2, 6, 26, (7*24+2), (28*24+2))


# the following numbers derive from normalisation of the MPPD numbers 0.12 and 0.52
tb.depo = 0.1875 
al.depo = 0.8125
ua.depo = 0
group1 <-list("exposure.concentration" = input.data[3,1], "exposure.time" = 2,
              "weight" = 277, "dep.ua" = 0, "dep.tb" = tb.depo*input.data[2,1], "dep.al" =  al.depo*input.data[2,1], 
              "Inhaled.vol.rate" = input.data[10,1])

params.group1 <- create.params(group1)

group2 <-list("exposure.concentration" = input.data[3,2], "exposure.time" = 2,
              "weight" = 277, "dep.ua" = 0, "dep.tb" =  tb.depo*input.data[2,2], "dep.al" = al.depo*input.data[2,2], 
              "Inhaled.vol.rate" = input.data[10,2])
params.group2 <- create.params(group2)

group3 <-list("exposure.concentration" = input.data[3,3], "exposure.time" = 2,
              "weight" = 277, "dep.ua" = 0, "dep.tb" =  tb.depo*input.data[2,3], "dep.al" = al.depo*input.data[2,3], 
              "Inhaled.vol.rate" = input.data[10,3])
params.group3 <- create.params(group3)

group4 <-list("exposure.concentration" = input.data[3,4], "exposure.time" = 2,
              "weight" = 277, "dep.ua" = 0, "dep.tb" =  tb.depo*input.data[2,4], "dep.al" = al.depo*input.data[2,4],
              "Inhaled.vol.rate" = input.data[10,4])
params.group4 <- create.params(group4)

group5 <-list("exposure.concentration" = input.data[3,5], "exposure.time" = 2,
              "weight" = 277, "dep.ua" = 0, "dep.tb" =  tb.depo*input.data[2,5], "dep.al" = al.depo*input.data[2,5], 
              "Inhaled.vol.rate" = input.data[10,5])
params.group5 <- create.params(group5)

params = matrix(c(params.group1, params.group2, params.group3,params.group4, params.group5), ncol=5)

x_fast <- 0.7
CLE_ur<- 0.7 
uptake <- 1e-03
k_alpc_tb <- 3e-03
k_lu_al <-  1e-02
k_al_lu <- 0.9
uptake_al <-  4
uptake_skel <- 27e-03
P_lu <- 400
P_ki <- 1.8
P_br <- 0.005
P_ut <- 9
P_skin <- 0.2
P_soft <- 0.08
P_li <- 0.35

sd <- 5

x_fast_std <- x_fast/sd
CLE_ur_std <- CLE_ur/sd
uptake_std <-uptake/sd #microg TiO2/g PCs
k_alpc_tb_std <- k_alpc_tb/sd
k_lu_al_std <-  k_lu_al/sd
k_al_lu_std <- k_al_lu/sd
uptake_al_std <- uptake_al/sd
uptake_skel_std <- uptake_skel/sd
P_lu_std <- P_lu/sd
P_ki_std <- P_ki/sd
P_br_std <- P_br/sd
P_ut_std <- P_ut/sd
P_skin_std <- P_skin/sd
P_soft_std <- P_soft/sd
P_li_std <- P_li/sd



eta_tr <- c(x_fast, CLE_ur, uptake, k_alpc_tb, k_lu_al, k_al_lu, uptake_al, uptake_skel, P_lu,
            P_ki, P_br, P_ut, P_skin, P_soft, P_li)

eta_tr_std <- c(x_fast_std, CLE_ur_std, uptake_std, k_alpc_tb_std, k_lu_al_std, k_al_lu_std,  uptake_al_std,
                 uptake_skel_std,
                P_lu_std, P_ki_std, P_br_std, P_ut_std, P_skin_std, P_soft_std, P_li_std)


#############################################################################################
set.seed(3843)
# function form 2 with an argument named `chain_id`
initf2 <- function(chain_id = 1,eta_tr, eta_tr_std) {
  theta <- rep(NA, length(eta_tr))
  eta_std <- rep(NA, length(eta_tr))
  eta <- rep(NA, length(eta_tr))
  for (i in 1:length(eta_tr)){
    eta_std[i] <- sqrt(log(((eta_tr_std[i]^2)/(eta_tr[i])^2)+1));
    eta[i] <- log(((eta_tr[i])^2)/sqrt((eta_tr_std[i]^2)+(eta_tr[i])^2));
    theta[i] <- rnorm(1,eta[i], eta_std[i]) 
  }
  
  sigma <- abs(rnorm(1,1.2,0.1))
  list("theta[1]" = theta[1],"theta[2]" = theta[2],"theta[3]" = theta[3],"theta[4]" = theta[4],"theta[5]" = theta[5],"theta[6]" = theta[6],
       "theta[7]" = theta[7], "theta[8]" = theta[8],"theta[9]" = theta[9],"theta[10]" = theta[10],"theta[11]" = theta[11],"theta[12]" = theta[12],
       "theta[13]" = theta[13],"theta[14]" = theta[14], "theta[15]" = theta[15])
}

# generate a list of lists to specify initial values
n_chains <- 4
inits <- lapply(1:n_chains, function(id) initf2(chain_id = id, eta_tr, eta_tr_std ))

DataList <- list( eta_tr = eta_tr,
                  eta_tr_std = eta_tr_std,
                  N_diff = 39,
                  N_compart = 13,
                  N_group = 5,
                  N_param = length(eta_tr),
                  params = params,
                  feces = feces,
                  urine = urine,
                  excreta_time = excreta_time,
                  biodist_time =  biodist_time ,
                  biodist_data = biodist_data[2:14],
                  t_init = 0,
                  m0 = rep(0,39),
                  rel_tol = 1e-06, 
                  abs_tol = 1e-06,
                  max_num_steps = 1e05)

tic = proc.time()

fit <- stan(file = 'fitting_stan.stan', data = DataList, iter = 1000, warmup=400, chains=4)# , init = inits)#,
#control = list(adapt_delta = 0.9))
options(max.print=5.5E5) 

#check_hmc_diagnostics(fit)
#check_divergences(fit)
#check_treedepth(fit)
#check_energy(fit)
#if I had set  control=list(max_treedepth=15) the the correct command would be check_treedepth(fit, 15) 

#library(shinystan)
#launch_shinystan(fit)
#print(fit)


#clock<-proc.time() - tic
#print(clock)
#pairs(fit, pars=c("par"))
#traceplot(fit,c("theta[2]"));
#stan_dens(fit,("theta[1]"),separate_chains = TRUE)
#exp_fit <- extract(fit)
#mean(exp_fit$theta[,1,2])