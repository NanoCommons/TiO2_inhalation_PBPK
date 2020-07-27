functions{
  real []  pbpk(real t,
                real[] M,
                real[] theta,
                real[] rdata,
                int[] idata) {
    
    
    real dMdt[39] ;
    
    real Vtis[10]; real Vven; real Vart; real V_blood; real Vcap[10]; real Wm_blood; real Wm[10];
    real Q[10]; real QTotal; 
    
    real Inhaled_vol_rate; real dep_ua; real dep_al;  real dep_tb; real external_conc; real inhal_duration;
    
    real x_fast; real P;
    real CLE_ur; real CLE_hep;   real CLE_ua; real CLE_muc; real CLE_muc_cor;
    real k_ab0; real k_de; real k_ab0_spl; real uptake; real Wm_al;
    real k_ua_br; real k_alpc_tb; real k_lu_al; real k_al_lu;  real uptake_al;
    real uptake_skel ;real uptake_spl ;
    
    
    real Clu_tis; real  Vlu_tis; real  Q_lu; real Clu_cap; real  Vlu_cap;
    real Cli_tis; real Vli_tis; real  Q_li; real Cli_cap; real Vli_cap;
    real Cspl_tis; real Vspl_tis; real  Q_spl; real Cspl_cap; real Vspl_cap;
    real Cki_tis; real Vki_tis; real  Q_ki; real Cki_cap; real Vki_cap;
    real Cht_tis; real Vht_tis; real  Q_ht; real Cht_cap; real Vht_cap;
    real Cbr_tis; real Vbr_tis; real  Q_br; real Cbr_cap; real Vbr_cap;
    real Cut_tis; real Vut_tis; real  Q_ut; real Cut_cap; real Vut_cap;
    real Cskel_tis; real Vskel_tis; real  Q_skel; real Cskel_cap; real Vskel_cap;
    real Cskin_tis; real Vskin_tis; real  Q_skin; real Cskin_cap; real Vskin_cap;
    real Csoft_tis; real Vsoft_tis; real  Q_soft;real Csoft_cap; real Vsoft_cap;
    real Cart; real Cven;
    
    real kluab;  real kliab; real ksplab; real kkiab; real khtab;
    real kbrab; real kutab; real kskelab; real kskinab; real ksoftab;
    real kbloodab; real kalab;
    
    real Wm_lu; real Wm_li; real Wm_spl; real Wm_ki;
    real Wm_ht; real Wm_br; real Wm_ut; real Wm_skel;
    real Wm_skin; real Wm_soft;
    
    real P_lu; real P_li; real P_spl; real P_ki;
    real P_ht; real P_br; real P_ut; real P_skel; real P_skin; real P_soft;real P_skel_soft_ht; real P_li_spl;
    
    real x_lu; real x_li; real x_spl; real x_ki;
    real x_ht; real x_br; real x_ut; real x_skel; real x_skin; real x_soft;
    real x_slow;
    
    Vlu_tis = rdata[1];
    Vli_tis = rdata[2];
    Vspl_tis = rdata[3];
    Vki_tis = rdata[4];
    Vht_tis = rdata[5];
    Vbr_tis = rdata[6];
    Vut_tis = rdata[7];
    Vskel_tis = rdata[8];
    Vskin_tis = rdata[9];
    Vsoft_tis = rdata[10];
    Vven = rdata[11];
    Vart = rdata[12];
    V_blood = rdata[13];
    
    Vlu_cap = rdata[14];
    Vli_cap = rdata[15];
    Vspl_cap = rdata[16];
    Vki_cap = rdata[17];
    Vht_cap = rdata[18];
    Vbr_cap = rdata[19];
    Vut_cap = rdata[20];
    Vskel_cap = rdata[21];
    Vskin_cap = rdata[22];
    Vsoft_cap = rdata[23];
    
    Q_lu = rdata[24];
    Q_li = rdata[25];
    Q_spl = rdata[26];
    Q_ki = rdata[27];
    Q_ht = rdata[28];
    Q_br = rdata[29];
    Q_ut = rdata[30];
    Q_skel = rdata[31];
    Q_skin = rdata[32];
    Q_soft = rdata[33];
    QTotal = rdata[34];
    
    Wm_blood = rdata[35];
    Wm_lu = rdata[36];
    Wm_li = rdata[37];
    Wm_spl = rdata[38];
    Wm_ki = rdata[39];
    Wm_ht = rdata[40];
    Wm_br = rdata[41];
    Wm_ut = rdata[42];
    Wm_skel = rdata[43];
    Wm_skin = rdata[44];
    Wm_soft = rdata[45];
    
    Inhaled_vol_rate = rdata[46];
    dep_ua = rdata[47];
    dep_tb = rdata[48];
    dep_al = rdata[49];
    external_conc = rdata[50];
    inhal_duration = rdata[51];
    
    CLE_ua = rdata[52];
    CLE_muc = rdata[53];
    k_ua_br = 0;
    Wm_al = 0.05;
    CLE_muc_cor = 1; 
    k_de =  6.7e-19;
    CLE_hep =  1.52e-03;   
    k_ab0 = 1.45; 
    k_ab0_spl = 0.5; 
    
    ///////////params///////////
      x_fast = exp(theta[1]);
      CLE_ur = exp(theta[2]);   
      uptake = exp(theta[3]); 
      k_alpc_tb = exp(theta[4]);  
      k_lu_al = exp(theta[5]); 
      k_al_lu = exp(theta[6]);
      uptake_al = exp(theta[7]);
      uptake_skel = exp(theta[8]);
      uptake_spl = uptake;
      
      P_skel_soft_ht = exp(theta[14]);
      P_li_spl = exp(theta[15]);
      
      P_lu = exp(theta[9]);
      P_li = P_li_spl;
      P_spl = P_li_spl;
      P_ki = exp(theta[10]); 
      P_ht = P_skel_soft_ht;
      P_br = exp(theta[11]);
      P_ut = exp(theta[12]);
      P_skel = P_skel_soft_ht;
      P_skin = exp(theta[13]);
      P_soft = P_skel_soft_ht;
      
      x_lu = x_fast;
      x_li = x_fast;
      x_spl = x_fast;
      x_ki = x_fast;
      x_ht =  x_fast;
      x_br =  x_fast;
      x_ut = x_fast;
      x_skel = x_fast;
      x_skin = x_fast;
      x_soft = x_fast;
      
      ///////////////
        // ODEs system //
        ///////////////
        // concentrations in tissues
      Clu_tis  =  M[6]/Vlu_tis;
      Clu_cap =  M[5]/Vlu_cap;
      Cli_tis  =  M[9]/Vli_tis;
      Cli_cap =  M[8]/Vli_cap;
      Cspl_tis  =  M[12]/Vspl_tis;
      Cspl_cap  =  M[11]/Vspl_cap;
      Cki_tis  =  M[15]/Vki_tis;
      Cki_cap  =  M[14]/Vki_cap;
      Cht_tis   =   M[18]/Vht_tis;
      Cht_cap  =  M[17]/Vht_cap;
      Cbr_tis  =  M[21]/Vbr_tis;
      Cbr_cap  =  M[20]/Vbr_cap;
      Cut_tis  =  M[24]/Vut_tis;
      Cut_cap  =  M[23]/Vut_cap;
      Cskel_tis  =  M[27]/Vskel_tis;
      Cskel_cap  =  M[26]/Vskel_cap;
      Cskin_tis  =  M[30]/Vskin_tis;
      Cskin_cap  =  M[29]/Vskin_cap;
      Csoft_tis  =  M[33]/Vsoft_tis;
      Csoft_cap  =  M[32]/Vsoft_cap;
      
      Cart = M[35]/(0.19*V_blood);
      Cven = M[36]/(0.81*V_blood);
      
      // Uptake rates by phagocytizing cells
      kluab = k_ab0*(1-(M[7]/(Wm_lu*uptake)));
      kliab = k_ab0*(1-(M[10]/(Wm_li*uptake)));
      ksplab = k_ab0_spl*(1-(M[13]/(Wm_spl*uptake_spl)));
      kkiab = k_ab0*(1-(M[16]/(Wm_ki*uptake)));
      khtab = k_ab0*(1-(M[19]/(Wm_ht*uptake)));
      kbrab = k_ab0*(1-(M[22]/(Wm_br*uptake)));
      kutab = k_ab0*(1-(M[25]/(Wm_ut*uptake)));
      kskelab = k_ab0*(1-(M[28]/(Wm_skel*uptake_skel)));
      kskinab = k_ab0*(1-(M[31]/(Wm_skin*uptake)));
      ksoftab = k_ab0*(1-(M[34]/(Wm_soft*uptake)));
      kbloodab = k_ab0*(1-(M[37]/(Wm_blood*uptake)));
      kalab = k_ab0*(1-(M[4]/(Wm_al*uptake_al)));
      
      //Upper airways
      if (t < inhal_duration){
        dMdt[1] = -(k_ua_br*M[1]) - (CLE_ua * M[1]) + (Inhaled_vol_rate * external_conc * dep_ua);
        
        //Trancheobronchial region
        dMdt[2] = (k_alpc_tb * M[4]) - (CLE_muc* CLE_muc_cor * M[2])+ (Inhaled_vol_rate * external_conc * dep_tb);
        
        //Alveolar region
        dMdt[3] = (k_lu_al * M[6]) - (k_al_lu * M[3]) - (M[3]*kalab - M[4]*k_de)+ (Inhaled_vol_rate * external_conc * dep_al);
      }
      else {
        dMdt[1] = -(k_ua_br*M[1]) - (CLE_ua * M[1]);
        
        //Trancheobronchial region
        dMdt[2] = (k_alpc_tb * M[4]) - (CLE_muc* CLE_muc_cor * M[2]);
        
        //Alveolar region
        dMdt[3] = (k_lu_al * M[6]) - (k_al_lu * M[3]) - (M[3]*kalab - M[4]*k_de);
      }
      dMdt[4] = (M[3]*kalab - M[4]*k_de) - (k_alpc_tb * M[4]);
      
      //Lungs
      dMdt[5]  =  (Cven*Q_lu)- (Clu_cap*Q_lu) + (x_lu*Q_lu)*Clu_tis/P_lu - (x_lu*Q_lu)*Clu_cap ; 
      dMdt[6]  =  - (x_lu*Q_lu)*Clu_tis/P_lu + (x_lu*Q_lu)*Clu_cap - (Vlu_tis*Clu_tis*kluab - M[7]*k_de) -
        (k_lu_al * M[6]) + (k_al_lu * M[3]) ; //Lung interstitium
      dMdt[7]  = (Vlu_tis*Clu_tis*kluab - M[7]*k_de); 
      
      //Liver
      dMdt[8]  = (Cart*Q_li) + (Cspl_cap*Q_spl) - (Cli_cap*(Q_li+Q_spl)) + (x_li*Q_li)*Cli_tis/P_li - (x_li*Q_li)*Cli_cap;  //capillary
      dMdt[9] = - (x_li*Q_li)*Cli_tis/P_li + (x_li*Q_li)*Cli_cap -(Vli_tis*Cli_tis*kliab - M[10]*k_de) - (CLE_hep*M[9]) ; //tissue
      dMdt[10]  = (Vli_tis*Cli_tis*kliab - M[10]*k_de) ; //Phagocytized
      
      //Spleen
      dMdt[11] = (Cart*Q_spl) - (Cspl_cap*Q_spl) +  (x_spl*Q_spl)*Cspl_tis/P_spl - (x_spl*Q_spl)*Cspl_cap ;  //capillary
      dMdt[12] =  - (x_spl*Q_spl)*Cspl_tis/P_spl + (x_spl*Q_spl)*Cspl_cap - (Vspl_tis*Cspl_tis*ksplab - M[13]*k_de)  ; //tissue
      dMdt[13]  = (Vspl_tis*Cspl_tis*ksplab - M[13]*k_de) ; //seq
      
      //Kidneys
      dMdt[14]  = (Cart*Q_ki) - (Cki_cap*Q_ki) +  (x_ki*Q_ki)*Cki_tis/P_ki - (x_ki*Q_ki)*Cki_cap- (CLE_ur*M[14]) ;  //capillary
      dMdt[15] =  - (x_ki*Q_ki)*Cki_tis/P_ki + (x_ki*Q_ki)*Cki_cap - (Vki_tis*Cki_tis*kkiab - M[16]*k_de)  ; //tissue
      dMdt[16]  = (Vki_tis*Cki_tis*kkiab - M[16]*k_de) ; //seq
      
      //Heart
      dMdt[17]  = (Cart*Q_ht) - (Cht_cap*Q_ht) +  (x_ht*Q_ht)*Cht_tis/P_ht - (x_ht*Q_ht)*Cht_cap ;  //capillary
      dMdt[18] =  - (x_ht*Q_ht)*Cht_tis/P_ht + (x_ht*Q_ht)*Cht_cap - (Vht_tis*Cht_tis*khtab - M[19]*k_de)  ; //tissue
      dMdt[19]  = (Vht_tis*Cht_tis*khtab - M[19]*k_de) ; //seq
      
      //Brain
      dMdt[20]  = (Cart*Q_br) - (Cbr_cap*Q_br) +  (x_br*Q_br)*Cbr_tis/P_br - (x_br*Q_br)*Cbr_cap ;  //capillary
      dMdt[21] =  - (x_br*Q_br)*Cbr_tis/P_br + (x_br*Q_br)*Cbr_cap - (Vbr_tis*Cbr_tis*kbrab - M[22]*k_de) + (k_ua_br*M[1])  ; //tissue
      dMdt[22]  =  (Vbr_tis*Cbr_tis*kbrab - M[22]*k_de) ; //seq
      
      //Uterus
      dMdt[23]  = (Cart*Q_ut) - (Cut_cap*Q_ut) +  (x_ut*Q_ut)*Cut_tis/P_ut - (x_ut*Q_ut)*Cut_cap ;  //capillary
      dMdt[24] =  - (x_ut*Q_ut)*Cut_tis/P_ut+ (x_ut*Q_ut)*Cut_cap - (Vut_tis*Cut_tis*kutab - M[25]*k_de)  ; //tissue
      dMdt[25]  = (Vut_tis*Cut_tis*kutab - M[25]*k_de) ; //seq
      
      //Skeleton
      dMdt[26]  = (Cart*Q_skel) - (Cskel_cap*Q_skel) +  (x_skel*Q_skel)*Cskel_tis/P_skel - (x_skel*Q_skel)*Cskel_cap ;  //capillary
      dMdt[27] =  - (x_skel*Q_skel)*Cskel_tis/P_skel + (x_skel*Q_skel)*Cskel_cap - (Vskel_tis*Cskel_tis*kskelab - M[28]*k_de)  ; //tissue
      dMdt[28]  = (Vskel_tis*Cskel_tis*kskelab - M[28]*k_de) ; //seq
      
      //Skin
      dMdt[29]  = (Cart*Q_skin) - (Cskin_cap*Q_skin) +  (x_skin*Q_skin)*Cskin_tis/P_skin - (x_skin*Q_skin)*Cskin_cap ;  //capillary
      dMdt[30] =  - (x_skin*Q_skin)*Cskin_tis/P_skin + (x_skin*Q_skin)*Cskin_cap - (Vskin_tis*Cskin_tis*kskinab - M[31]*k_de)  ; //tissue
      dMdt[31]  = (Vskin_tis*Cskin_tis*kskinab - M[31]*k_de) ; //seq
      
      //Soft
      dMdt[32]  = (Cart*Q_soft) - (Csoft_cap*Q_soft) +  (x_soft*Q_soft)*Csoft_tis/P_soft - (x_soft*Q_soft)*Csoft_cap ;  //capillary
      dMdt[33] =  - (x_soft*Q_soft)*Csoft_tis/P_soft + (x_soft*Q_soft)*Csoft_cap - (Vsoft_tis*Csoft_tis*ksoftab - M[34]*k_de)  ; //tissue
      dMdt[34]  = (Vsoft_tis*Csoft_tis*ksoftab - M[34]*k_de) ; //seq
      
      //Blood
      dMdt[35] = (Clu_cap*Q_lu) - Cart* (Q_li + Q_spl+ Q_ki + Q_ht + Q_br + Q_ut + Q_skel + Q_skin + Q_soft) -  
        (0.19*V_blood*Cart*kbloodab - 0.19* M[37]*k_de); 
      dMdt[36] = (Cli_cap*(Q_li+Q_spl)) +  (Cki_cap*Q_ki) + (Cht_cap*Q_ht) + (Cbr_cap*Q_br) + (Cut_cap*Q_ut) + (Cskel_cap*Q_skel) +
        (Cskin_cap*Q_skin) + (Csoft_cap*Q_soft) - (Cven*Q_lu)  -  (0.81*V_blood*Cven*kbloodab - 0.81* M[37]*k_de);
      dMdt[37] =  ((0.19*V_blood*Cart + 0.81*V_blood*Cven)*kbloodab - M[37]*k_de) ;
      
      //Feces
      dMdt[38] = (CLE_hep*M[9]) + (CLE_ua * M[1]) + (CLE_muc* CLE_muc_cor * M[2])  ;
      
      //Urine
      dMdt[39] = (CLE_ur*M[14]) ;
      
      return dMdt;
  }
}

//////////////////////////////////////////////////////////////////////////
  data{
    int<lower=0> N_param;                 // Number of parameters to be estimated
    int<lower=0> N_compart;               //number of observed compartments
    int<lower=0> N_diff;              // number of differential equations
    int<lower=0> N_group;              // number of subjects or groups
    real  biodist_time[N_group];
    real  excreta_time[8];
    real  feces[8];
    real  urine[8];
    real  biodist_data[N_group,N_compart];
    real m0[N_diff];           // Initial concentration in compartments
    real t_init;                  // Initial time
    real eta_tr[N_param];
    real eta_tr_std[N_param];
    real  params[53,N_group];      // Matrix containing the individual parameters
    real rel_tol;
    real abs_tol;
    real max_num_steps;
  }

////////////////////////////////////////////////////////////////////////////
  transformed data {
    real rdata[0];
    int idata[0];
    //vector[N_param]  eta_tr_std ;
    vector[N_param]  eta_std ;
    vector[N_param]  eta ;
    vector [N_param] H;                //covariance matrix
    
    for (i in 1:N_param){
      //eta_tr_std[i] = eta_tr[i];
      eta_std[i]= sqrt(log(((eta_tr_std[i]^2)/(eta_tr[i])^2)+1));
      eta[i]=log(((eta_tr[i])^2)/sqrt((eta_tr_std[i]^2)+(eta_tr[i])^2));
      H[i] =  eta_std[i];
    }
    
    
  // print(integrate_ode_bdf(pbpk, m0, t_init, biodist_time,to_array_1d(eta[:]), params[:,5], idata,
                     //       rel_tol, abs_tol, max_num_steps))
    
  }
//////////////////////////////////////////////////////////////////
  
  parameters{
    
    real<lower=0>  sigma[3];
    vector [N_param] theta_tr;
  }

////////////////////////////////////////////////////////////////////
  
  model{
    real ode_res[N_group,N_diff, N_group];
    real excreta[8,N_diff];
    real ode_hat[N_group,N_compart];
    real feces_hat[8];
    real urine_hat[8];

    
    
    //priors
    sigma[:] ~ normal(0,1);
  
    theta_tr[:] ~normal(eta[:],H[:]);
    
    
    //likelihood
    for(i in 1:N_group){
      ode_res[:,:,i] = integrate_ode_bdf(pbpk, m0, t_init, biodist_time,to_array_1d(theta_tr[:]), params[:,i], idata,
                                         rel_tol, abs_tol, max_num_steps);
      
    } 
    excreta[:,:] = integrate_ode_bdf(pbpk, m0, t_init, excreta_time,to_array_1d(theta_tr[:]), params[:,5], idata,
                                     rel_tol, abs_tol, max_num_steps);
    
    for(i in 1:N_group){
      //Amount in Trachea
      //ode_hat[i,1] = ode_res[i,2,i];
      
      //Amount in BALF
      ode_hat[i,1] = ode_res[i,3,i];
      
      //Amount in BALC
      ode_hat[i,2] = ode_res[i,4,i];
      
      //Amount in lavaged lungs
      ode_hat[i,3] = ode_res[i,6,i] + ode_res[i,7,i];
      
      //Amount in liver
      ode_hat[i,4] = ode_res[i,9,i] + ode_res[i,10,i];
      
      //Amount in spleen
      ode_hat[i,5] = ode_res[i,12,i] + ode_res[i,13,i];
      
      //Amount in kidneys
      ode_hat[i,6] = ode_res[i,15,i] + ode_res[i,16,i];
      
      //Amount in heart
      ode_hat[i,7] = ode_res[i,18,i] + ode_res[i,19,i];
      
      //Amount in brain
      ode_hat[i,8] = ode_res[i,21,i] + ode_res[i,22,i];
      
      //Amount in uterus
      ode_hat[i,9] = ode_res[i,24,i] + ode_res[i,25,i];
      
      //Amount in skeleton
      ode_hat[i,10] = ode_res[i,27,i] + ode_res[i,28,i];
      
      //Amount in skin
      ode_hat[i,11] = ode_res[i,30,i] + ode_res[i,31,i];
      
      //Amount in soft tissue
      ode_hat[i,12] = ode_res[i,33,i] + ode_res[i,34,i];
      
      //Amount in blood
      ode_hat[i,13] = ode_res[i,35,i] + ode_res[i,36,i] + ode_res[i,37,i];
      
    }
    feces_hat[:] = excreta[:,38];
    urine_hat[:] = excreta[:,39];
    
    //error grouping based on mass magnitude
    for (j in 1:3){
      to_vector(biodist_data[:,j]) ~ normal( to_vector(ode_hat[:,j]),sigma[1]);
    }
    to_vector(biodist_data[:,4]) ~ normal( to_vector(ode_hat[:,4]),sigma[2]);
    to_vector(biodist_data[:,5]) ~ normal( to_vector(ode_hat[:,5]),sigma[3]);
    to_vector(biodist_data[:,6]) ~ normal( to_vector(ode_hat[:,6]),sigma[2]);
   
     for (j in 7:9){
      to_vector(biodist_data[:,j]) ~ normal( to_vector(ode_hat[:,j]),sigma[3]);
    }
     for (j in 10:13){
      to_vector(biodist_data[:,j]) ~ normal( to_vector(ode_hat[:,j]),sigma[2]);
     }

    to_vector(feces[:]) ~ normal(to_vector(feces_hat[:]),sigma[1]);
    to_vector(urine[:]) ~ normal(to_vector(urine_hat[:]),sigma[1]);
    
  }

generated quantities{
  vector [N_param] theta;

  theta = exp(theta_tr);
}
