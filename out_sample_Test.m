%Out-of-sample test: scripted by Zehao Cao
%% Calculate the out-of-sample cost
Cost_out = zeros(1,N_Omega_normal);% predefine the out-of-sample cost
if Method == 1
%get the "value" of optimal decision
for n_chp = 1:N_chp
    eval(['P_chp_',num2str(n_chp),'= value(P_chp_',num2str(n_chp),');']);
    eval(['af_chp_',num2str(n_chp),'= zeros(n_T,1);']);
end
for n_bs = 1:N_bs
    eval(['P_bsc_',num2str(n_bs),'= value(P_bsc_',num2str(n_bs),');']);
    eval(['P_bsd_',num2str(n_bs),'= value(P_bsd_',num2str(n_bs),');']);
    eval(['af_bsc_',num2str(n_bs),'= value(af_bsc_',num2str(n_bs),');']);
    eval(['af_bsd_',num2str(n_bs),'= value(af_bsd_',num2str(n_bs),');']);
end
for n_gt = 1:N_gt
    eval(['P_gt_',num2str(n_gt),'= value(P_gt_',num2str(n_gt),');']);
    eval(['af_gt_',num2str(n_gt),'= value(af_gt_',num2str(n_gt),');']);
end
for n_eb = 1:N_eb
    eval(['P_eb_',num2str(n_eb),'= value(P_eb_',num2str(n_eb),');']);
    eval(['af_eb_',num2str(n_eb),'= zeros(n_T,1);']);
end
P_p2g = value(P_p2g);af_p2g = value(af_p2g); 
P_buy = value(P_buy);P_sell = value(P_sell);G_buy = value(G_buy);

else
%get the "value" of optimal decision
for n_chp = 1:N_chp
    eval(['P_chp_',num2str(n_chp),'= value(P_chp_',num2str(n_chp),');']);
    eval(['af_chp_',num2str(n_chp),'= value(af_chp_',num2str(n_chp),');']);
end
for n_bs = 1:N_bs
    eval(['P_bsc_',num2str(n_bs),'= value(P_bsc_',num2str(n_bs),');']);
    eval(['P_bsd_',num2str(n_bs),'= value(P_bsd_',num2str(n_bs),');']);
    eval(['af_bsc_',num2str(n_bs),'= value(af_bsc_',num2str(n_bs),');']);
    eval(['af_bsd_',num2str(n_bs),'= value(af_bsd_',num2str(n_bs),');']);
end
for n_gt = 1:N_gt
    eval(['P_gt_',num2str(n_gt),'= value(P_gt_',num2str(n_gt),');']);
    eval(['af_gt_',num2str(n_gt),'= value(af_gt_',num2str(n_gt),');']);
end
for n_eb = 1:N_eb
    eval(['P_eb_',num2str(n_eb),'= value(P_eb_',num2str(n_eb),');']);
    eval(['af_eb_',num2str(n_eb),'= value(af_eb_',num2str(n_eb),');']);
end
P_p2g = value(P_p2g);af_p2g = value(af_p2g); 
P_buy = value(P_buy);P_sell = value(P_sell);G_buy = value(G_buy);
end
%calculate out-of-sample cost 
for i = 1:N_Omega_normal
Cost_out(i) = C_buy_E'*P_buy - C_sell_E(1:n_T)'*P_sell + C_chp_1'*P_chp_1 + C_chp_2'*P_chp_2 + C_chp_3'*P_chp_3...
    - C_BS_ch'*(P_bsc_1 + P_bsc_2) + C_BS_dch'*(P_bsd_1 + P_bsd_2) + C_P2G'*P_p2g + C_EB'*(P_eb_1 + P_eb_2) + C_Gas'*G_buy;%发电成本+购气成本
Cost_out(i) = Cost_out(i) + c_e_car'*(e_chp*(P_chp_1 + P_chp_2 + P_chp_3) + e_gt*(P_gt_1 + P_gt_2) + e_gas*G_L(1:n_T) + e_grid*P_buy - e_p2g*P_p2g); %碳成本
Cost_out(i) = Cost_out(i) - C_chp_1'*(diag(Wdr(i,1:n_T))*af_chp_1) - C_chp_2'*(diag(Wdr(i,1:n_T))*af_chp_2) - C_chp_3'*(diag(Wdr(i,1:n_T))*af_chp_3)...
    - C_BS_ch'*(diag(Wdr(i,1:n_T))*(af_bsc_1+af_bsc_2)) - C_BS_dch'*(diag(Wdr(i,1:n_T))*(af_bsd_1+af_bsd_2)) + ...
   C_P2G'*(diag(Wdr(i,1:n_T))*af_p2g) + C_EB'*(diag(Wdr(i,1:n_T))*(af_eb_1+af_eb_2));
Cost_out(i) = Cost_out(i) + c_e_car(1:n_T)'*(-e_chp*(diag(Wdr(i,1:n_T))*af_chp_1 + diag(Wdr(i,1:n_T))*af_chp_2 + diag(Wdr(i,1:n_T))*af_chp_3)...
    - e_gt*diag(Wdr(i,1:n_T))*(af_gt_1+af_gt_2) -e_p2g*diag(Wdr(i,1:n_T))*af_p2g);
end
%% Verify the violation of DRJCCs
r_chp = 2*N_chp;r_bs = 4*N_bs;r_eb = 2*N_eb;r_gt = 2*N_gt;r_car = 1;r_vol = 2;%set of DRJCCs
for n = 1:size(Wdr,1)
    Del_chp{n} = zeros(n_T,r_chp);%define chp violation
    Del_gt{n} = zeros(n_T,r_gt);%define gt violation
    Del_eb{n} = zeros(n_T,r_eb);%define eb violation
    Del_bs{n} = zeros(n_T,r_bs);%define bs violation
    Del_car{n} = zeros(1,r_car);%define car violation
    Del_vol{n} = zeros(n_T,4);%define vol violation,it is sufficient to test nodes 14-17 for the rest of the nodes have no possibility of crossing the boundary
end
Pmax_EB_1 = Pmax_EB;Pmax_EB_2 = Pmax_EB;%parameters standardization   
Pmin_EB_1 = Pmin_EB;Pmin_EB_2 = Pmin_EB;
Pmax_BSC_1 = Pmax_BS_C;Pmax_BSC_2 = Pmax_BS_C;Pmax_BSD_1 = Pmax_BS_D;Pmax_BSD_2 = Pmax_BS_D;
Pmin_BSC_1 = Pmin_BS_C;Pmin_BSC_2 = Pmin_BS_C;Pmin_BSD_1 = Pmin_BS_D;Pmin_BSD_2 = Pmin_BS_D;
violate = zeros(n_T,6);%define violation matrix
for n = 1:size(Wdr,1)
    % calculate chp violation
    for n_chp = 1:N_chp
    eval([' Del_chp{n}(:,2*n_chp-1) = Pmax_chp_',num2str(n_chp),' - (P_chp_',num2str(n_chp),'-diag(Wdr(n,:))*af_chp_',num2str(n_chp),');']); 
    eval([' Del_chp{n}(:,2*n_chp) = -Pmin_chp_',num2str(n_chp),' + (P_chp_',num2str(n_chp),'-diag(Wdr(n,:))*af_chp_',num2str(n_chp),');']); 
    end
    % calculate chp violation
    for n_gt = 1:N_gt
    eval([' Del_gt{n}(:,2*n_gt-1) = Pmax_GT_',num2str(n_gt),' - (P_gt_',num2str(n_gt),'-diag(Wdr(n,:))*af_gt_',num2str(n_gt),');']); 
    eval([' Del_gt{n}(:,2*n_gt) = -Pmin_GT_',num2str(n_gt),' + (P_gt_',num2str(n_gt),'-diag(Wdr(n,:))*af_gt_',num2str(n_gt),');']); 
    end
    % calculate eb violation
    for n_eb = 1:N_eb
    eval([' Del_eb{n}(:,2*n_eb-1) = Pmax_EB_',num2str(n_eb),' - (P_eb_',num2str(n_eb),'+ diag(Wdr(n,:))*af_eb_',num2str(n_eb),');']); 
    eval([' Del_eb{n}(:,2*n_eb) = -Pmin_EB_',num2str(n_eb),' + (P_eb_',num2str(n_eb),'+ diag(Wdr(n,:))*af_eb_',num2str(n_eb),');']); 
    end    
    % calculate bs violation
    for n_bs = 1:N_bs
    eval([' Del_bs{n}(:,4*n_bs-3) = Pmax_BSD_',num2str(n_bs),' - (P_bsd_',num2str(n_bs),'-diag(Wdr(n,:))*af_bsd_',num2str(n_bs),');']); 
    eval([' Del_bs{n}(:,4*n_bs-2) = -Pmin_BSD_',num2str(n_bs),' + (P_bsd_',num2str(n_bs),'-diag(Wdr(n,:))*af_bsd_',num2str(n_bs),');']); 
    eval([' Del_bs{n}(:,4*n_bs-1) = Pmax_BSC_',num2str(n_bs),' - (P_bsc_',num2str(n_bs),'+ diag(Wdr(n,:))*af_bsc_',num2str(n_bs),');']); 
    eval([' Del_bs{n}(:,4*n_bs) = -Pmin_BSC_',num2str(n_bs),' + (P_bsc_',num2str(n_bs),'+ diag(Wdr(n,:))*af_bsc_',num2str(n_bs),');']); 
    end    
    % calculate car violation
    Del_car{n}(1,1) = (1-lambda_e)*e_max-value(  repmat(e_chp,1,n_T)*(P_chp_1 + P_chp_2 + P_chp_3)+ repmat(e_gt,1,n_T)*(P_gt_1 + P_gt_2) + repmat(e_grid,1,n_T)*P_buy...
        + repmat(e_gas,1,n_T)*G_L - repmat(e_p2g,1,n_T)*(P_p2g)+ Wdr(n,:)*( - e_chp*(af_chp_1+af_chp_2+af_chp_3) + e_gt*(- af_gt_1 - af_gt_2)...
        - e_p2g*af_p2g ));
    % calculate vol violation
     Del_vol{n}(:,:) = -1*transpose((value(V_bus(:,15:18)')*V_0 - (Bra_B(14:17,:)*((R*Bra_B*Conf_WT'+ X*Bra_B*Conf_WT'*0.75)*Wdr(n,:)*1e3))...
         - (Bra_B(14:17,:)*(- R*Bra_B*value(af_nbus_VP(:,:)') - X*Bra_B*value(afQ_nbus(:,:)'))*diag(Wdr(n,:)*1e3)) - ones(4, 1)*1.06*V_0^2)/V_0);
end
for n = 1:size(Wdr,1)
        [idx_chp,~] = find([Del_chp{n}] < -1e-2);%Determine if there are DRJCCs that are not satisfied in t
        [idx_gt,~] = find([Del_gt{n}] < -1e-2);
        [idx_bs,~] = find([Del_bs{n}] < -1e-2);
        [idx_eb,~] = find([Del_eb{n}] < -1e-2);
        idx_car = find([Del_car{n}] < -1e-2);
        [idx_vol,~] = find([Del_vol{n}] < -1e-2);
            if isempty(idx_chp) ~= 1
                idx_chp = unique(idx_chp);%Delete duplicate value
                sgn_chp = zeros(n_T,1);sgn_chp(idx_chp) = 1;
                violate(:,1) = violate(:,1) + sgn_chp;%Add number of violation
                sgn_chp = zeros(n_T,1);
            end
            if isempty(idx_gt) ~= 1
                idx_gt = unique(idx_gt);%Delete duplicate value
                sgn_gt= zeros(n_T,1);sgn_gt(idx_gt) = 1;
                violate(:,2) = violate(:,2) + sgn_gt;%Add number of violation
                sgn_gt = zeros(n_T,1);
             end
            if isempty(idx_bs) ~= 1
                idx_bs = unique(idx_bs);%Delete duplicate value
                sgn_bs = zeros(n_T,1);sgn_bs(idx_bs) = 1;
                violate(:,3) = violate(:,3) + sgn_bs;%Add number of violation
                sgn_bs = zeros(n_T,1);
            end
            if isempty(idx_eb) ~= 1
                idx_eb = unique(idx_eb);%Delete duplicate value
                sgn_eb = zeros(n_T,1);sgn_eb(idx_eb) = 1;
                violate(:,4) = violate(:,4) + sgn_eb;%Add number of violation
                sgn_eb = zeros(n_T,1);
            end
            if isempty(idx_car) ~= 1
                violate(1,5) = violate(1,5)+ 1;
            end
            if isempty(idx_vol) ~= 1
                idx_vol = unique(idx_vol);%Delete duplicate value
                sgn_vol = zeros(n_T,1);sgn_vol(idx_vol) = 1;
                violate(:,6) = violate(:,6) + sgn_vol;%Add number of violation
                sgn_vol = zeros(n_T,1);
            end
end
%% Present the out-of-sample test results
Cost_out_mean = mean(Cost_out);% out-of-sample cost
Cost_out_median = median(Cost_out);
COST = sprintf('The out-of-sample cost of MEMG is  %d',Cost_out_mean);
   disp(COST);
if for_extr == 0   
violate_pro = (sum(sum(violate)))/(size(Wdr,1)*5*n_T+1);
vio = sprintf('The violation probability of DRJCCs is %d %%',round(violate_pro*100,3));
disp(vio);
elseif for_extr == 1
violate_pro_max = max([mean(violate(:,1:4))/size(Wdr,1) violate(1,5)/(size(Wdr,1)*n_T) mean(violate(:,6))/size(Wdr,1)]);
delta_pr = violate_pro_max;
vio_max = sprintf('The violation probability of DRJCCs is %d %%',round(violate_pro_max*100,3));
    disp(vio_max);    
end

