%Prepare "unclean" sample set: scripted by Zehao Cao
%% Initialization
disp("preparing data set...");
[PW_max,PW_min,PV_max,PV_min] = DERs_org_data;
xi_max_tran = s_max;%define boundry of samples generation
%% Generation of wind power forecast error samples (generated using actual data from Elia) with outliers
WT_erros_org = readmatrix('Elia&DSO_DERs_data.xlsx','sheet','Feuil1', 'Range','U2:U48001');%get original samples
Wdf = reshape(WT_erros_org(1:N_Omega_normal*n_T),[n_T,N_Omega_normal]);%generate 2000 trainning samples 
Wdf = Wdf'*s_max;%scaling samples
rng(123,'twister');%set seed
Wdf_normal = zeros(N_extreme,n_T);
for t = 1:n_T
Wdf_normal(1:0.4*N_extreme,t) = xi_max_tran*normrnd(-0.1,0.3,[1,0.4*N_extreme]);
Wdf_normal(0.4*N_extreme+1:N_extreme,t) = xi_max_tran*normrnd(-0.15,0.3,[1,0.6*N_extreme]);
%Wdf_normal(1:N_extreme,t) = -0.5+ -0.5*xi_max_tran + 2*0.5*xi_max_tran*rand(1,N_extreme);
Wdf_normal(Wdf_normal <= -xi_max_tran) = 0.15 - xi_max_tran; 
Wdf_normal(Wdf_normal >= xi_max_tran) = xi_max_tran;
end
WDF_D = zeros(N_extreme+N_Omega_normal, n_T);
A_index = 1;
B_index = 1;
C_index = 1;
while B_index <= size(Wdf_normal, 1)
    WDF_D(C_index:C_index+9, :) = Wdf(A_index:A_index+9, :); 
    A_index = A_index + 10;
    C_index = C_index + 10;
    WDF_D(C_index, :) = Wdf_normal(B_index, :);
    B_index = B_index + 1;
    C_index = C_index + 1;
end
WDF_D(C_index:end, :) = Wdf(A_index:end, :);

%Wdf = [Wdf;Wdf_normal];%insert outlies
Wdf = WDF_D;%insert outlies
%random_order = randperm(size(Wdf,1));
%Wdf = Wdf(random_order, :);%randperm Wdf
Wdr = Wdf;%generate 2000 test samples, for Wdf_S<<Wdf, so let wdr==wdf
Wdf_mean = mean(Wdf);%calculate mean of samples

disp("*data set has been defined*");
