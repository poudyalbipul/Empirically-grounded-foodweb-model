
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 28 12:03:10 2025

@author: bipul
"""

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import Foodweb_model_qss as fm

beta = fm.beta

keywords = ['mean_length_list', 'mean_trophic_list', 'max_stabil_list',
            'total_biomass_list'
            ] 

def simulate_foodweb(beta_values, beta_index, n_spec=150, beta = beta):
    
    dict_list = {key: [] for key in keywords}
    #initialize a list
    
    for i in range(50):  # number of realizations
        species_id = fm.generate_species(n_spec, random=True,e = 0.85, B0=1e-6, K =1e3, beta = beta)
        #initialize list
        keywords_2 = ['length_list', 'trophic_list', 'stabil_list', 'biomass_list'] 

        dict_runs = {key: [] for key in keywords_2}
        
        for beta_value in beta_values:
            beta[beta_index] = beta_value
            species_id = fm.change_temperature(species_id, beta=beta)
            mu, A = fm.compute_LV_param(species_id, beta=beta)            
            mu_eff, A_eff = fm.reduce_LV(mu,A)
            N_0 = np.full(n_spec-1, 1)            
            #N_0[0] = species_id["K"]
            #N_0_log = np.log(N_0)
            t_eval = np.linspace(0, 2000, 51)
            sol = solve_ivp(fm.LV_model, [0, 2000], 
                            t_eval= t_eval, 
                            y0=N_0, method="LSODA", args=(mu_eff, A_eff))  
            
            #compute number of survivors and total population
            survivors = sol.y[:, -1] > 1e-3            
            biomass = sum(num for num in sol.y[0:,-1] if num > 1e-3)            

            #compute stability
            A = A_eff #/ np.linalg.norm(A, ord='fro')      #normalize A
            ind_s = np.where(sol.y[:,-1] > 1e-3)[0]
            A_surv = A_eff[ind_s[:, np.newaxis], ind_s]
            J = -np.diag(sol.y[ind_s,-1])@A_surv       
            J_norm = J 
            #np.corrcoef(off_diag_values, off_diag_values_T)[0,1]
            eigenvalues = np.linalg.eigvals(J_norm).real
            #stability = np.nanmax(eigenvalues) 
                        
            n_full = fm.compute_links(species_id).shape[0]
            surv_full = np.zeros(n_full, dtype = bool)
            surv_full[1:] = survivors
            surv_full[0] = True
            
            trophic_level = fm.compute_trophic_level(species_id, surv_full)                            
            #s = fm.compute_trophic_level(species_id,survivors)
            mean_trophic_level = np.nanmean(trophic_level)
            
         
            
            if sum(survivors) == 0 or sum(survivors) == n_spec-1:
                for key in keywords_2:
                    dict_runs[key].append(np.nan)
               
            else:
                dict_runs['length_list'].append(sum(survivors))            
                dict_runs['trophic_list'].append(mean_trophic_level)                                    
                dict_runs['biomass_list'].append(biomass)     
                dict_runs['stabil_list'].append(
                    np.nan if eigenvalues.size == 0 else np.nanmax(eigenvalues))
        
        for k1, k2 in zip(keywords, keywords_2):
            dict_list[k1].append(dict_runs[k2]) 
      
    #proportion_formed = 1- (np.isnan(dict_list['mean_length_list']).sum(axis=0)/len(dict_list['mean_length_list']))
    
    #return community metrics
    return (np.nanmedian(dict_list['mean_length_list'], axis=0),            
            np.nanmean(dict_list['mean_trophic_list'], axis=0),        
            np.nanmedian(dict_list['max_stabil_list'], axis = 0),             
            np.nanmedian(dict_list['total_biomass_list'], axis = 0)          
            )

# Initialize parameters
#Change parameter values from 1/2 to 2 times the reference value
beta_00 = np.sort(np.geomspace(beta[0, 0] / 2, 2 * beta[0, 0], 11))
beta_01 = np.sort(np.geomspace(beta[0, 1] / 2, 2 * beta[0, 1], 11))
beta_10 = np.sort(np.geomspace(beta[1, 0] / 2, 2 * beta[1, 0], 11))
beta_11 = np.sort(np.geomspace(beta[1, 1] / 2 ,2 * beta[1, 1], 11))
beta_20 = np.sort(np.geomspace(beta[2, 0] / 2, 2 * beta[2, 0], 11))
beta_21 = np.sort(np.geomspace(beta[2, 1] / 2, 2 * beta[2, 1], 11))
beta_30 = np.sort(np.geomspace(beta[3, 0] / 2, 2 * beta[3, 0], 11))
beta_31 = np.sort(np.geomspace(beta[3, 1] / 2, 2 * beta[3, 1], 11))
beta_40 = np.sort(np.geomspace(beta[4, 0] / 2, 2 * beta[4, 0], 11))
beta_50 = np.sort(np.geomspace(beta[5, 0] / 2, 2 * beta[5, 0], 11))
beta_60 = np.sort(np.geomspace(beta[6, 0] / 2, 2 * beta[6, 0], 11))
beta_61 = np.sort(np.geomspace(beta[6, 1] / 2, 2 * beta[6, 1], 11))
beta_42 = np.sort(np.geomspace(beta[4, 2] / 2, 2 * beta[4, 2], 11))

#beta_40 = np.sort(np.linspace(0.5,4, 20))

#Change parameter values up to aquatic-terrestrial reference value.
"""beta_00 = np.sort(np.linspace(beta[0, 0] - 2.59, beta[0, 0] + 2.59, 10))
beta_01 = np.sort(np.linspace(beta[0, 1] - 0.13667, beta[0, 1] + 0.13667, 10))
beta_10 = np.sort(np.linspace(beta[1, 0] - 0.9865, beta[1, 0] + 0.9865, 10))
beta_11 = np.sort(np.linspace(beta[1, 1] - 0.0001, beta[1, 1] + 0.0001, 10))
beta_20 = np.sort(np.linspace(beta[2, 0] - 1.0883, beta[2, 0] + 1.0883, 10))
beta_21 = np.sort(np.linspace(beta[2, 1] - 0.03114, beta[2, 1] + 0.03114, 10))
beta_30 = np.sort(np.linspace(beta[3, 0] - 0, beta[3, 0] + 0, 10))
beta_31 = np.sort(np.linspace(beta[3, 1] - 0, beta[3, 1] + 0, 10))
beta_50 = np.sort(np.linspace(beta[5, 0] - 7.47, beta[5, 0] + 7.47, 10))
beta_60 = np.sort(np.linspace(beta[6, 0] - 0.49, beta[6, 0] + 0.49, 10))
beta_61 = np.sort(np.linspace(beta[6, 1] - 0.0495, beta[6,1] + 0.0495, 10))
"""
#Change parameter values up to +- 0.5 of reference value.
"""beta_00 = np.sort(np.linspace(beta[0, 0] - 0.5, beta[0, 0] + 0.5, 10))
beta_01 = np.sort(np.linspace(beta[0, 1] - 0.5, beta[0, 1] + 0.5, 10))
beta_10 = np.sort(np.linspace(beta[1, 0] - 0.5, beta[1, 0] + 0.5, 10))
beta_11 = np.sort(np.linspace(beta[1, 1] - 0.5, beta[1, 1] + 0.5, 10))
beta_20 = np.sort(np.linspace(beta[2, 0] - 0.5, beta[2, 0] + 0.5, 10))
beta_21 = np.sort(np.linspace(beta[2, 1] - 0.5, beta[2, 1] + 0.5, 10))
beta_30 = np.sort(np.linspace(beta[3, 0] - 0.5, beta[3, 0] + 0.5, 10))
beta_31 = np.sort(np.linspace(beta[3, 1] - 0.5, beta[3, 1] + 0.5, 10))
beta_40 = np.sort(np.linspace(0, 0.5, 20))
beta_50 = np.sort(np.linspace(beta[5, 0] - 0.5, beta[5, 0] + 0.5, 10))
beta_60 = np.sort(np.linspace(beta[6, 0] - 0.5, beta[6, 0] + 0.5, 10))
beta_61 = np.sort(np.linspace(beta[6, 1] - 0.5, beta[6, 1] + 0.5, 10))"""


# Run simulations
results = simulate_foodweb(beta_60, beta_index=(6, 0))
 
#reset beta values for next simulation run
beta[0,0] = -3.1503  
beta[0, 1] = 0.234519  
beta[1, 0] = 1.6645 
beta[1, 1] = 0.09433
beta[2, 0] = -1.1012   
beta[2, 1] = -0.066577
beta[3, 0] = 0.9637
beta[3, 1] = -0.04
beta[4,0] = 1
beta[5, 0] = -11.83
beta[6, 0] = 0.54
beta[6, 1] = 0.005 
beta[4,2] = 0.85

#Save in an excel file

datapm = pd.DataFrame({
                         "change":fm.normalize_x(beta_60), "diversity":results[0], 
                         "biomass":results[3], "mean_tl": results[1],"stability": results[2]
                         })
    
    
excel_path = "betaai_sensitivity.xlsx"
datapm.to_excel(excel_path, index=False)

##########################################     
