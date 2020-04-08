# Assembly-emergence-via-triplet-STDP
Code for the paper "Autonomous emergence of connectivity assemblies via spike triplet interactions"


This file explains how to use the code for the publication:
Lisandro Montangie, Christoph Miehl, Julijana Gjorgjieva (2020) PLoS Computational Biology


Each figure can be reproduced by generating the respective data via "Triplet_plasticity_code_etam" (or "..._taue", or "..._external_input") and the data is processed and visualized with "K_means_clustering_code" and "Graph_measures_code".

The code of the simulation is based on the code from the publication: Ravid Tannenbaum, Burak (2016) PLoS Comput Biol 12(8):e1005056.

The graph measures are from the Brain Connectivity Toolbox (brain-connectivity-toolbox.net): Rubinov M, Sporns O (2010) NeuroImage 52:1059-69. 


We provide a detailed explanation of how to generate each of the figures:

Fig. 6:
Use the code "Triplet_plasticity_code_etam" to generate weight matrices for different values of the modulation parameter etam. This will generate a file 'data_triplet_STDP.mat' which can be used to visualize the results. The k-means clustering and the fraction of bidirectional connections can be performed with running "K_means_clustering_code".

Fig. 7:
Use the code "Triplet_plasticity_code_etam" to generate weight matrices for different values of the modulation parameter etam. This will generatea file 'data_triplet_STDP.mat' which can be used to visualize the results. To calculate the graph measures, run "Graph_measures_code". 

Fig. 8:
Same procedure as Fig. 6 and 7. To get the respective cases of only cross-cov, no loops, cross-cov + loops comment the respective terms in the code "Triplet_plasticity_code_etam". For details of which lines to exclude in each respective case see Fig. 5 and text in main manuscript.

Fig. 9:
Use the code "Triplet_plasticity_code_taue" to generate weight matrices for different values of the first membrane time constant taue. Then use "K_means_clustering_code" and "Graph_measures_code" to visualize the data. Note to change the value of the x-axis to taue in these two scripts".

Fig. 10:
Use the code "Triplet_plasticity_code_external_input" to generate weight matrices for different values of external input correlations. Then use "K_means_clustering_code" and "Graph_measures_code" to visualize the data.

Fig. 11: 
Same as Fig. 8. To get the respective cases, change the delta values in the "Triplet_plasticity_code_etam" script. 

Suppl. Fig. 1:
Numerical integration of Eq. 2 of the main text with N=12 and etam=13.

Suppl. Fig. 2: 
Same as Fig. 8.

Suppl. Fig. 3:
Use the code "Triplet_plasticity_code_etam" but instead of looping over etam, change the external input rate 'mu0' (Suppl. Fig. 3A) or include a standard deviation term for the external input rate (Suppl. Fig. 3B).

Suppl. Fig. 4: 
Use the code "Triplet_plasticity_code_etam" but instead of looping over etam change the parameter 'inh_mult' and visualize it via "Graph_measures_code".

Suppl. Fig. 5: 
Use the code "Triplet_plasticity_code_etam" and save the firing rate parameter 'firing_rate' through the whole simulation. 
