%Prepare "clean" sample set: scripted by Zehao Cao
%% Initialization
disp("preparing data set...");
[PW_max,PW_min,PV_max,PV_min] = DERs_org_data;
xi_max_tran = s_max;%define boundry of samples generation
if choice == 0%select generation method 1:GMM  2:from Elia
%% Generation of wind power forecast error samples (generated using a Gaussian Mixture Model)
rng(123,'twister');%set seed
Wdf_normal = zeros(N_Omega_normal,n_T);%define original sampleset
Wdr = zeros(N_Omega_normal,n_T);%define sampleset for out of sample test
for t = 1:n_T
%generate samples ↓↓↓↓↓
    Wdf_normal(1:0.3*N_Omega_normal,t) = xi_max_tran*normrnd(-0.1,0.1,[1,0.3*N_Omega_normal]);
    Wdf_normal(0.3*N_Omega_normal+1:0.8*N_Omega_normal,t) = xi_max_tran*normrnd(0.1,0.15,[1,0.5*N_Omega_normal]);
    Wdf_normal(0.8*N_Omega_normal+1:N_Omega_normal,t) = xi_max_tran*normrnd(0.1,0.2,[1,0.2*N_Omega_normal]);
    Wdf_normal(Wdf_normal <= -xi_max_tran) = -xi_max_tran; 
    Wdf_normal(Wdf_normal >= -0.2+xi_max_tran) = xi_max_tran;
    Wdr(1:0.3*N_Omega_normal,t) = xi_max_tran*normrnd(-0.1,0.1,[1,0.3*N_Omega_normal]);
    Wdr(0.3*N_Omega_normal+1:0.8*N_Omega_normal,t) = xi_max_tran*normrnd(0.1,0.15,[1,0.5*N_Omega_normal]);
    Wdr(0.8*N_Omega_normal+1:N_Omega_normal,t) = xi_max_tran*normrnd(0.1,0.2,[1,0.2*N_Omega_normal]);
    Wdr(Wdr <= -xi_max_tran) = -xi_max_tran; 
    Wdr(Wdr >= -0.2+xi_max_tran) = xi_max_tran;
end
Wdf = Wdf_normal;
Wdf0 = Wdf;
Wdf_mean = mean(Wdf);%calculate mean of samples
else
%% Generation of wind power forecast error samples (generated using actual data from Elia)
WT_erros_org = readmatrix('Elia&DSO_DERs_data.xlsx','sheet','Feuil1', 'Range','U2:U48001');%get original samples
Wdf = reshape(WT_erros_org(1:N_Omega_normal*n_T),[n_T,N_Omega_normal]);%generate 2000 trainning samples 
Wdf = Wdf'*s_max;%scaling samples
Wdr = Wdf;%generate 2000 test samples, for Wdf_S<<Wdf, so let wdr==wdf
Wdf_mean = mean(Wdf);%calculate mean of samples
end
%% Representative sub set generation
rng(12,'twister');%set seed
idx_Wdf = randi(size(Wdf,1),N_Sub,1);
Wdf_S = Wdf(idx_Wdf(1:N_Sub),:);
    
disp("*data set has been defined*");
