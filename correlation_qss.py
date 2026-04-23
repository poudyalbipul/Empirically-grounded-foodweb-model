#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 14 14:29:00 2025

@author: bipul
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import solve_ivp
import itertools
from sklearn.metrics.pairwise import cosine_similarity
import Foodweb_model_qss as fm
import networkx as nx
from scipy.stats import entropy


beta = fm.beta

def simulate_communities(n_times, n_spec, beta = beta ):
    
    keylist = ['population_list', 'link_density_list', 'diversity_list', 'stability_list', 'mtl_list']
    
    dict_list = {key: [] for key in keylist}
    
   
    for i in range(n_times):
        species_id = fm.generate_species(n_spec, e = 0.85, random=True, B0=1e-6, K = 1e3, beta= beta )
        species_id = fm.change_temperature(species_id, beta=beta)
        mu, A = fm.compute_LV_param(species_id, beta=beta)     
        mu_eff, A_eff = fm.reduce_LV(mu,A)
        N_0 = np.full(n_spec-1, 1)
        t_eval = np.linspace(0, 2000, 51)
        sol = solve_ivp(fm.LV_model, [0, 2000], 
                        t_eval= t_eval, 
                        y0=N_0, method="LSODA", args=(mu_eff, A_eff))  
        
        # compute total population and species richness
        survivors = sol.y[:, -1] > 1e-3            
        #m_i = np.array(species_id["m_i"])
        survivors_nonbasal = sol.y[:,-1] > 1e-3
        #biomass = sum(num for num in sol.y[1:,-1] if num > 1)     
        #bodysize = m_i[1:][survivors_nonbasal] * sol.y[1:,-1][survivors_nonbasal]
 
        
        if np.any(survivors_nonbasal):
            # At least one non-basal species survived
       
            surviving_pops = sol.y[0:,-1][survivors_nonbasal]
            diversity = np.sum(survivors)
            
            nonbasal_population = np.nansum(surviving_pops) 
            basal_population = sol.y[0,-1]
                        
            ind_s = np.where(sol.y[:,-1] > 1e-3)[0]
            A_surv = A_eff[np.ix_(ind_s, ind_s)]
            J_norm = -np.diag(sol.y[:,-1][ind_s])@A_surv
    
            eigenvalues = np.linalg.eigvals(J_norm).real
            stability = np.nanmax(eigenvalues)
        
        else:
            # No non-basal species survived - food web collapsed         
            diversity = np.nan
            stability = np.nan
            nonbasal_population = np.nan
            basal_population = np.nan
        
        #compute number of links and link density
        adjacency_matrix = fm.compute_links(species_id).T[1:,1:][survivors][:,survivors]
        link_density = np.sum(adjacency_matrix)/(diversity)
        
        n_full = fm.compute_links(species_id).shape[0]
        surv_full = np.zeros(n_full, dtype = bool)
        surv_full[1:] = survivors
        surv_full[0] = True
        
        trophic_level = fm.compute_trophic_level(species_id, surv_full)                            
        mean_trophic_level = np.nanmean(trophic_level)
   
            
        if sum(survivors) == 1 or sum(survivors) == n_spec:
            for key in keylist:
                dict_list[key].append(np.nan)
        else:
            dict_list['population_list'].append(nonbasal_population + basal_population)
            dict_list['link_density_list'].append(link_density)
            dict_list['diversity_list'].append(diversity)
            dict_list['stability_list'].append(stability)
            dict_list['mtl_list'].append(mean_trophic_level)
           
    return (dict_list['diversity_list'], dict_list['population_list'], 
            dict_list['link_density_list'], dict_list['stability_list'],     
            dict_list['mtl_list']) 
          

Results = simulate_communities(500, 150)



simulation_data = pd.DataFrame({"Diversity": Results[0], "Total_biomass": Results[1],                                                                
                                "Link_density": Results[2] , "Stability":Results[3], 
                                "Mean_trophic_level": Results[4]                                                                                                                                                 
                               })




excel_path = "community_outcomes_trial.xlsx"
simulation_data.to_excel(excel_path, index=False)