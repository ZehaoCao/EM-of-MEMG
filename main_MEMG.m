%EM of MEMG: scripted by Zehao Cao
clear all; rehash toolboxcache;yalmip('clear');
%% Initialization
n_T = 24; % dispactch periods
dt = 1;% dispactch interval
B = 10;% number of K-means
Method = 4;% 1:DO   2:RO   3:W-DR0   4:C-DR0
choice = 1;% 1:data set from Elia&DSO  0:data set from GMM
for_extr = 0;% 0:clean samples 1:samples with outliers
T_pruning = 0.1;%percentage of sample-pruning: [0,1]
%% Parameters of DRO
N_Omega_normal = 2000; %number of Trainning samples
N_extreme = 200;%number of extreme samples
N_Beta_normal = 2000; %number of Test samples
N_Sub = 200;%size of representative sub-set
rho = 1e-3;% radius of Wasserstein sphere
u_max = 0.1;% confidence interval of first-order moment: [0,1]
s_max = 4;% Max forecast errors
risk_factor = 0.05; %risk coefficient: [0,1]
%% Processing
switch for_extr
%%===============================clean samples===============================%%
    case 0 %clean samples
data_treat_clean%introduce the original data:non outliers→generate the representative sub-set
initialize_parameter;%introduce the parameters of MEMG
    switch Method
        case 1
        Model_NoUncertain_MEMG;%DO-method
        case 2
        Model_RObased_MEMG;%RO-method
        case 3
        Model_W_DRObased_MEMG;%W-DRO method
        case 4
        Model_C_DRObased_MEMG;%C-DRO method
    end
out_sample_Test;%Calculate the out-of-sample cost&verify DRJCCs

%%========================samples with outliers==============================%%
    case 1 %samples with outliers
        Intra = 1;%difine number of interations
        delta_pr = 0;%difine initial violation
        data_treat_dirty%introduce the original data:with outliers
        
      while delta_pr <= risk_factor && Intra <= 8
        initialize_parameter;%introduce the parameters of MEMG
        sample_pruning;%generate the representative sub-set
          
        switch Method
            case 3
                Model_W_DRObased_MEMG;%W-DRO method
            case 4
                Model_C_DRObased_MEMG;%C-DRO method
          end
        out_sample_Test;%Calculate the out-of-sample cost&verify DRJCCs
 
        Intra = Intra + 1;
      end
end
