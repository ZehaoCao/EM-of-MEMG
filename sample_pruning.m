% Sample-pruning: scripted by Zehao Cao
%% Sample-pruning
V_test = V_0*1e3;X_test = X*1e3;R_test = R*1e3;% to Kohms, to KV, for avoiding Numerical problems
if for_extr == 1 && Intra >= 2
        for n = 1:size(Wdf,1)%in-sample data size
        Del0_chp{n} = zeros(n_T,r_chp);%define chp violation
        Del0_gt{n} = zeros(n_T,r_gt);%define gt violation
        Del0_eb{n} = zeros(n_T,r_eb);%define eb violation
        Del0_bs{n} = zeros(n_T,r_bs);%define bs violation
        Del0_car{n} = zeros(1,r_car);%define car violation
        Del0_vol{n} = zeros(n_T,4);%define vol violation,it is sufficient to test nodes 14-17 for the rest of the nodes have no possibility of crossing the boundary
        end
%%========================test with in-sample data==============================%%
        for n = 1:size(Wdf,1)
    % calculate chp violation
    for n_chp = 1:N_chp
    eval([' Del0_chp{n}(:,2*n_chp-1) = Pmax_chp_',num2str(n_chp),' - (P_chp_',num2str(n_chp),'-diag(Wdf(n,:))*af_chp_',num2str(n_chp),');']); 
    eval([' Del0_chp{n}(:,2*n_chp) = -Pmin_chp_',num2str(n_chp),' + (P_chp_',num2str(n_chp),'-diag(Wdf(n,:))*af_chp_',num2str(n_chp),');']); 
    end
    % calculate chp violation
    for n_gt = 1:N_gt
    eval([' Del0_gt{n}(:,2*n_gt-1) = Pmax_GT_',num2str(n_gt),' - (P_gt_',num2str(n_gt),'-diag(Wdf(n,:))*af_gt_',num2str(n_gt),');']); 
    eval([' Del0_gt{n}(:,2*n_gt) = -Pmin_GT_',num2str(n_gt),' + (P_gt_',num2str(n_gt),'-diag(Wdf(n,:))*af_gt_',num2str(n_gt),');']); 
    end
    % calculate eb violation
    for n_eb = 1:N_eb
    eval([' Del0_eb{n}(:,2*n_eb-1) = Pmax_EB_',num2str(n_eb),' - (P_eb_',num2str(n_eb),'+ diag(Wdf(n,:))*af_eb_',num2str(n_eb),');']); 
    eval([' Del0_eb{n}(:,2*n_eb) = -Pmin_EB_',num2str(n_eb),' + (P_eb_',num2str(n_eb),'+ diag(Wdf(n,:))*af_eb_',num2str(n_eb),');']); 
    end    
    % calculate bs violation
    for n_bs = 1:N_bs
    eval([' Del0_bs{n}(:,4*n_bs-3) = Pmax_BSD_',num2str(n_bs),' - (P_bsd_',num2str(n_bs),'-diag(Wdf(n,:))*af_bsd_',num2str(n_bs),');']); 
    eval([' Del0_bs{n}(:,4*n_bs-2) = -Pmin_BSD_',num2str(n_bs),' + (P_bsd_',num2str(n_bs),'-diag(Wdf(n,:))*af_bsd_',num2str(n_bs),');']); 
    eval([' Del0_bs{n}(:,4*n_bs-1) = Pmax_BSC_',num2str(n_bs),' - (P_bsc_',num2str(n_bs),'+ diag(Wdf(n,:))*af_bsc_',num2str(n_bs),');']); 
    eval([' Del0_bs{n}(:,4*n_bs) = -Pmin_BSC_',num2str(n_bs),' + (P_bsc_',num2str(n_bs),'+ diag(Wdf(n,:))*af_bsc_',num2str(n_bs),');']); 
    end    
    % calculate car violation
    Del0_car{n}(1,1) = (1-lambda_e)*e_max-value(  repmat(e_chp,1,n_T)*(P_chp_1 + P_chp_2 + P_chp_3)+ repmat(e_gt,1,n_T)*(P_gt_1 + P_gt_2) + repmat(e_grid,1,n_T)*P_buy...
        + repmat(e_gas,1,n_T)*G_L - repmat(e_p2g,1,n_T)*(P_p2g)+ Wdf(n,:)*( - e_chp*(af_chp_1+af_chp_2+af_chp_3) + e_gt*(- af_gt_1 - af_gt_2)...
        - e_p2g*af_p2g ));
    % calculate vol violation
     Del0_vol{n}(:,:) = -1*transpose((value(V_bus(:,15:18)')*V_test - (Bra_B(14:17,:)*((R_test*Bra_B*Conf_WT'+ X_test*Bra_B*Conf_WT'*0.75)*Wdf(n,:)*1e3))...
         - (Bra_B(14:17,:)*(- R_test*Bra_B*value(af_nbus_VP(:,:)') - X_test*Bra_B*value(afQ_nbus(:,:)'))*diag(Wdf(n,:)*1e3)) - ones(4, 1)*1.06*V_test^2)/V_test);
    end
    Del_sum = zeros(6,size(Wdf,1));%define violation amplitude
    eva = zeros(size(Wdf,1),1);
     for n = 1:size(Wdf,1)
        Del_sum(1,n) = sum(Del0_chp{n}(Del0_chp{n}<-1e-3),'all')/(Pmax_chp_2(1) - Pmin_chp_2(1));
        Del_sum(2,n) = sum(Del0_gt{n}(Del0_gt{n}<-1e-3),'all')/(Pmax_GT_1(1)-Pmin_GT_1(1));
        Del_sum(3,n) = sum(Del0_eb{n}(Del0_eb{n}<-1e-3),'all')/(Pmax_EB(1)-Pmin_EB(1));
        Del_sum(4,n) = sum(Del0_bs{n}(Del0_bs{n}<-1e-3),'all')/(Pmax_BS_C(1)-Pmin_BS_C(1));
        Del_sum(5,n) = sum(Del0_car{n}(Del0_car{n}<-1e-3),'all')/(e_max*(1-lambda_e));
        Del_sum(6,n) = sum(Del0_vol{n}(Del0_vol{n}<-1e-3),'all');
        eva(n) = sum(Del_sum(:,n));% calculate the total violation amplitude
     end 
    [eva_amp,idx_eva_Wdf] = sortrows(eva); 
    eva_amp = eva_amp(eva_amp <= -3e-2);
    idx_sub_Wdf = idx_eva_Wdf( 1:round(T_pruning*(size(eva_amp))) );%get T% index of Wdf
    Wdf(idx_sub_Wdf,:) = [];%remove T% samples from Wdf
end

%% Representative sub set generation
method_sub = 0;%select method of sub sample set generation 1:random  2:cluster
tic
if method_sub == 0
rng(123,'twister');%set seed
%fix initial centers of clustering
%if np == 1
    %if data == 1
        %load('ini_kmeans_dirty.mat');
    %else
        %load('ini_kmeans_clean.mat');
    %end
    %initial_centroids = Wdf_C{2};
%end
for b = 1:B
    Wdf_C{b} = zeros(N_Sub,n_T);%define empty set to store samples
end
Wdf_S = zeros(N_Sub,n_T);%define representative sub set
sumD = zeros(B,N_Sub);
%B clustering is performed, then the data set that best represents the original sample is selected
for b = 1:B
    %if np ==1
    %[~,Wdf_C{b},sumD(b,:)] = kmeans(Wdf,N_Sub,'Start', initial_centroids);
    %else
    [~,Wdf_C{b},sumD(b,:)] = kmeans(Wdf,N_Sub);
    %end
end
[~,idx_sumD_min] = min(sum(sumD'));
Wdf_S = Wdf_C{idx_sumD_min};
%Wdf_S = Wdf(1:N_Sub,:);
toc
else
    rng(212,'twister');%set seed
    idx_Wdf = randi(size(Wdf,1),N_Sub,1);
    Wdf_S = Wdf(idx_Wdf(1:N_Sub),:);
end