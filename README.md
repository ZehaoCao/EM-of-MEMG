# EM-of-MEMG
Credible Joint Chance-Constrained Low-Carbon Energy Management for Multi-energy Microgrids
This document is intended to help you understand our codes, and we will describe each file on the following section.
1.	The MATLAB code “main_MEMG.m”
This file defines most of the parameters that appear in the manuscript. Detailed information can be found in the program comments, such as the range of values for the parameters and the choice of methods: 
2.	The MATLAB code “initialize_parameter.m”
This file introduces the operation parameters of the MEMG, such as the power limit of the CHP units and so on. Detailed information can be found in the code comments.
3.	The MATLAB code “LDPF_parameter.m”
This file is forced to generated the power network parameters of the MEMG, such as the Resistance matrix X and connectivity matrix B, etc. Detailed information can be found in the code comments.
4.	The MATLAB code “case33mg.m”
This file introduces the modified IEEE 33-bus distribution system. Detailed information can be found in the code comments.
5.	The MATLAB code “DERs_org_data.m”
This file introduces the original forecast power output of the DERs within MEMG. Detailed information can be found in the code comments.
6.	The MATLAB code “data_treat_clean.m”
This file introduces the uncertain forecast error samples belonging to two different sample sets. Notably, “clean” means these two sets do not contain extreme samples. Detailed information can be found in the code comments.
7.	The MATLAB code “data_treat_dirty.m”
This file introduces the uncertain forecast error samples belonging to a sample set that contains about 10% extreme samples. Detailed information can be found in the code comments.
8.	The MATLAB code “out_sample_Test.m”
This file is used for out-of-sample test of the optimal decisions, including testing for DRJCC violation and calculating out-of-sample operation costs. Detailed information can be found in the code comments.
9.	The MATLAB code “Model_C_DRObased_MEMG.m”
This file introduces the energy management model based on C-DRO, which proposed in our paper. Detailed information can be found in the code comments.
10.	The MATLAB code “Model_W_DRObased_MEMG.m”
This file introduces the energy management model based on W-DRO, which is used in the comparisons. Detailed information can be found in the code comments.
The MATLAB code “Model_RObased_MEMG.m”
This file introduces the energy management model based on RO, which is used in the comparisons. Detailed information can be found in the code comments.
11.	The MATLAB code “Model_NoUncertain_MEMG”
This file introduces the energy management model based on DO, which is used in the comparisons. Detailed information can be found in the code comments.
12.	The MATLAB code “sample_pruning.m”
This file introduces the sample-pruning algorithm proposed in our paper. Detailed information can be found in the code comments.
13.	The Excel document “Elia&DSO_DERs_data.xlsx”
This file introduces the dataset used in the case study, which is obtained from Elia. The relevant source is attached below:
https://www.elia.be/en/grid-data/power-generation/wind-power-generation

