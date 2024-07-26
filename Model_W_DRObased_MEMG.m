% W-DRO method: scripted by Zehao Cao
yalmip('clear');
%% Support set
suppH_tran = [eye(1);-eye(1)];
supph_tran = [repmat(s_max,1,1);repmat(s_max,1,1)];
Wdf_max = repmat(s_max,1,n_T);%forecast errors max positive deviation
Wdf_min = -repmat(s_max,1,n_T);%forecast errors max negative deviation
Wdf_up = max(Wdf_S);%first-order moment max positive deviation
Wdf_down = min(Wdf_S);%first-order moment max negative deviation
%% Loading power network data
mpc = loadcase('case33mg');%modified IEEE 33-BUS
node = mpc.bus;
P_D = repmat(node(:,3),1,n_T);P_D = P_D*diag(P_L(1:n_T));%active load matrix
Q_D = repmat(node(:,4),1,n_T);Q_D = Q_D*diag(P_L(1:n_T));%reactive load matrix
%% Define dicision variables
for n_chp = 1 : N_chp
eval(['P_chp_',num2str(n_chp),'= sdpvar(n_T,1);']);% P&Q of CHP
eval(['af_chp_',num2str(n_chp),'= sdpvar(n_T,1);']);
eval(['Q_chp_',num2str(n_chp),'= sdpvar(n_T,1);']);
eval(['afQ_chp_',num2str(n_chp),'= sdpvar(n_T,1);']);% adjust factor of CHP
end
for n_gt = 1 : N_gt
eval(['P_gt_',num2str(n_gt),'= sdpvar(n_T,1);']);
eval(['af_gt_',num2str(n_gt),'= sdpvar(n_T,1);']);
end
for n_eb = 1 : N_eb
eval(['P_eb_',num2str(n_eb),'= sdpvar(n_T,1);']);
eval(['af_eb_',num2str(n_eb),'= sdpvar(n_T,1);']);
end
for n_bs = 1 : N_bs
eval(['P_bsc_',num2str(n_bs),'= sdpvar(n_T,1);']);
eval(['P_bsd_',num2str(n_bs),'= sdpvar(n_T,1);']);
eval(['af_bsc_',num2str(n_bs),'= sdpvar(n_T,1);']);
eval(['af_bsd_',num2str(n_bs),'= sdpvar(n_T,1);']);
eval(['af_bs_',num2str(n_bs),'= sdpvar(n_T,1);']);
eval(['E_bs_',num2str(n_bs),'= sdpvar(n_T,1);']);
end 
P_p2g = sdpvar(n_T,1);
af_p2g = sdpvar(n_T,1);
P_buy = sdpvar(n_T,1);%elec purchase
P_sell = sdpvar(n_T,1);%elec sell
G_buy = sdpvar(n_T,1);%gaspurchase
P_flow_out = sdpvar(n_T,n_bus);Q_flow_out = sdpvar(n_T,n_bus);%powerflow of lines "in"
P_flow_in = sdpvar(n_T,n_bus);Q_flow_in = sdpvar(n_T,n_bus);%powerflow of lines "out"
V_bus = sdpvar(n_T,n_bus);%nodal voltage
af_nbus_VP = sdpvar(n_T,n_bus-1);%nodal ajust factor
af_nbus_VN = sdpvar(n_T,n_bus-1);
afQ_nbus = sdpvar(n_T,n_bus-1);
%% ===============================Define General Constraints=============================== %%
con_General = [];
disp("Modeling MEMG as W-DRO...");
V_0 = V_0*1e3;X = X*1e3;R = R*1e3;% to Kohms, to KV, for avoiding Numerical problems
for t = 1:n_T
    %power balance lim
    con_General = con_General + [ P_buy(t) - P_sell(t) + P_chp_1(t) + P_chp_2(t) + P_chp_3(t) + PW_1(t) + PW_2(t) + PW_3(t) + (n_bus-1)*PV(t)...
        + P_gt_1(t) + P_gt_2(t) == sum(P_D(:,t)) + P_eb_1(t) + P_eb_2(t) + P_bsc_1(t) + P_bsc_2(t) - P_bsd_1(t) - P_bsd_2(t)];
    con_General = con_General + [Q_chp_1(t) + Q_chp_2(t) + Q_chp_3(t) + 0.75*PW_1(t) + 0.75*PW_2(t) + 0.75*PW_3(t)...
        == sum(Q_D(:,t))];    
    %active power flow lim
    P_dev = [P_chp_1(t); P_chp_2(t); P_chp_3(t);-P_bsc_1(t);P_bsd_1(t);-P_bsc_2(t);P_bsd_2(t);-P_eb_1(t);-P_eb_2(t);P_gt_1(t);P_gt_2(t);-P_p2g(t)];%decision vectors stacking
    P_DER = [PW_1(t); PW_2(t); PW_3(t);PV(t)];%DERs' vectors stacking
    for n = 1:n_bus
        if n == 1
        con_General = con_General + [ P_flow_in(t,n) == P_buy(t) - P_sell(t) ]+[ P_flow_out(t,n) == P_flow_in(t,n)];%power flow of root
        elseif n == 3||n == 19
        con_General = con_General + [P_flow_in(t,3) + P_flow_in(t,19) == P_flow_out(t,2)]+...
            [P_flow_in(t,n)+ Conf(1:12,n)'*P_dev + Conf(13:16,n)'*P_DER == P_flow_out(t,n) + P_D(n,t)];%power flow of forked lines
        elseif n == 4||n == 23
        con_General = con_General + [P_flow_in(t,4) + P_flow_in(t,23) == P_flow_out(t,3)]+...
            [P_flow_in(t,n)+ Conf(1:12,n)'*P_dev + Conf(13:16,n)'*P_DER == P_flow_out(t,n) + P_D(n,t)];%power flow of forked lines
        elseif n == 7||n == 26
        con_General = con_General + [P_flow_in(t,7) + P_flow_in(t,26) == P_flow_out(t,6)]+...
            [P_flow_in(t,n)+ Conf(1:12,n)'*P_dev + Conf(13:16,n)'*P_DER == P_flow_out(t,n) + P_D(n,t)];%power flow of forked lines
        elseif n == 18||n == 22||n == 33
        con_General = con_General + [P_flow_in(t,n) == P_flow_out(t,n-1)] + [P_flow_in(t,n)+ Conf(1:12,n)'*P_dev + Conf(13:16,n)'*P_DER == P_D(n,t)] + [P_flow_out(t,n) == 0];%power flow of other lines     
        else    
        con_General = con_General + [P_flow_in(t,n) == P_flow_out(t,n-1)] + [P_flow_in(t,n)+ Conf(1:12,n)'*P_dev + Conf(13:16,n)'*P_DER == P_flow_out(t,n) + P_D(n,t)];%power flow of other lines 
        end
    end
    %reactive power flow lim
    Q_dev = [Q_chp_1(t); Q_chp_2(t); Q_chp_3(t);];%decision vectors stacking
    Q_DER = [0.75*PW_1(t); 0.75*PW_2(t); 0.75*PW_3(t)];%DERs' vectors stacking
    for n = 1:n_bus
        if n == 1
        con_General = con_General + [ Q_flow_in(t,n) == 0 ]+[ Q_flow_out(t,n) == 0];%power flow of root
        elseif n == 3||n == 19
        con_General = con_General + [Q_flow_in(t,3) + Q_flow_in(t,19) == Q_flow_out(t,2)]+...
            [Q_flow_in(t,n)+ Conf_Q(1:3,n)'*Q_dev + Conf_Q(4:6,n)'*Q_DER == Q_flow_out(t,n) + Q_D(n,t)];%power flow of forked lines
        elseif n == 4||n == 23
        con_General = con_General + [Q_flow_in(t,4) + Q_flow_in(t,23) == Q_flow_out(t,3)]+...
            [Q_flow_in(t,n)+ Conf_Q(1:3,n)'*Q_dev + Conf_Q(4:6,n)'*Q_DER == Q_flow_out(t,n) + Q_D(n,t)];%power flow of forked lines
        elseif n == 7||n == 26
        con_General = con_General + [Q_flow_in(t,7) + Q_flow_in(t,26) == Q_flow_out(t,6)]+...
            [Q_flow_in(t,n)+ Conf_Q(1:3,n)'*Q_dev + Conf_Q(4:6,n)'*Q_DER == Q_flow_out(t,n) + Q_D(n,t)];%power flow of forked lines
        elseif n == 18||n == 22||n == 33
        con_General = con_General + [Q_flow_in(t,n) == Q_flow_out(t,n-1)]+ [Q_flow_in(t,n)+ Conf_Q(1:3,n)'*Q_dev + Conf_Q(4:6,n)'*Q_DER == Q_D(n,t)] + [Q_flow_out(t,n) == 0];%power flow of other lines    
        else
        con_General = con_General + [Q_flow_in(t,n) == Q_flow_out(t,n-1)]+ [Q_flow_in(t,n)+ Conf_Q(1:3,n)'*Q_dev + Conf_Q(4:6,n)'*Q_DER == Q_flow_out(t,n) + Q_D(n,t)];%power flow of other lines 
        end
    end
    %nodal voltage lim
     for n = 1:n_bus
         if n == 1
         con_General = con_General + [ V_bus(t,n) == V_0 ];%vlotage of root node
        elseif n == 19
        con_General = con_General + [V_bus(t,n) == V_bus(t,2)-(R(n-1,n-1)*P_flow_in(t,n)*1e3+X(n-1,n-1)*Q_flow_in(t,n)*1e3)/V_0];%voltage of forked nodes
        elseif n == 23
        con_General = con_General + [V_bus(t,n) == V_bus(t,3)-(R(n-1,n-1)*P_flow_in(t,n)*1e3+X(n-1,n-1)*Q_flow_in(t,n)*1e3)/V_0];%voltage of forked nodes
        elseif n == 26
        con_General = con_General + [V_bus(t,n) == V_bus(t,6)-(R(n-1,n-1)*P_flow_in(t,n)*1e3+X(n-1,n-1)*Q_flow_in(t,n)*1e3)/V_0];%voltage of forked nodes
        else
        con_General = con_General + [V_bus(t,n) == V_bus(t,n-1)-(R(n-1,n-1)*P_flow_in(t,n)*1e3+X(n-1,n-1)*Q_flow_in(t,n)*1e3)/V_0];%voltage of other nodes
        end
     end   
    %ensuer variables non-negative
    for n_chp = 1 : N_chp
    eval([' con_General = con_General + [P_chp_',num2str(n_chp),'(t) >= 0];']);
    eval([' con_General = con_General + [0 <= af_chp_',num2str(n_chp),'(t) <= 1];']);
    eval([' con_General = con_General + [0 <= afQ_chp_',num2str(n_chp),'(t) <= 1];']);
    end
    for n_gt = 1 : N_gt
    eval([' con_General = con_General + [P_gt_',num2str(n_gt),'(t) >= 0];']);
    eval([' con_General = con_General + [0 <= af_gt_',num2str(n_gt),'(t) <= 1];']);
    end
    for n_eb = 1 : N_eb
    eval([' con_General = con_General + [P_eb_',num2str(n_eb),'(t) >= 0];']);
    eval([' con_General = con_General + [0 <= af_eb_',num2str(n_eb),'(t) <= 1];']);
    end
    for n_bs = 1 : N_bs
    eval([' con_General = con_General + [P_bsc_',num2str(n_bs),'(t) >= 0];']);
    eval([' con_General = con_General + [0 <= af_bsc_',num2str(n_bs),'(t) <= 1];']);
    eval([' con_General = con_General + [P_bsd_',num2str(n_bs),'(t) >= 0];']);
    eval([' con_General = con_General + [0 <= af_bsd_',num2str(n_bs),'(t) <= 1];']);
    eval([' con_General = con_General + [0 <= af_bs_',num2str(n_bs),'(t) <= 1];']);
    eval([' con_General = con_General + [E_bs_',num2str(n_bs),'(t) >= 0];']);
    end
    con_General = con_General + [P_p2g(t) >= 0] + [0 <= af_p2g(t) <= 1];
    con_General = con_General + [700>= G_buy(t) >= 0];
    con_General = con_General + [P_buy(t) >= 0] + [P_sell(t) >= 0];
    %LDR lim, all t to 1
    con_General = con_General + [af_chp_1(t)  + af_chp_2(t)+ af_chp_3(t) + af_eb_1(t) + af_eb_2(t)...
        + af_gt_1(t) + af_gt_2(t) + af_p2g(t) + af_bsc_1(t) + af_bsc_2(t) + af_bsd_1(t) + af_bsd_2(t) == 1];
    con_General = con_General + [afQ_chp_1(t) + afQ_chp_2(t)+ afQ_chp_3(t) == 1];
        
    %heat balance lim
    con_General = con_General + [eta_EB*(af_eb_1(t) + af_eb_2(t)) - H_E*(af_chp_1(t) + af_chp_2(t) + af_chp_3) == 0 ];
    con_General = con_General + [H_E*(P_chp_1(t) + P_chp_2(t) + P_chp_3(t)) + eta_EB*(P_eb_1(t) + P_eb_2(t)) == H_L(t,t)];
end
%% ===============================Define Robust Constraints=============================== %%
con_Robust = [];
for t = 1:n_T
    %ramping cons of CHP
    if t>=2% from t=2
            for n_chp = 1 : N_chp
            eval([' con_Robust = con_Robust + [ P_chp_',num2str(n_chp),'(t)  - P_chp_',num2str(n_chp),'(t-1) - af_chp_',num2str(n_chp),'(t)*Wdf_min(t) + af_chp_',num2str(n_chp),...
                '(t-1)*Wdf_max(t-1) <= ramp_chp_',num2str(n_chp),'];']);
            eval([' con_Robust = con_Robust + [ P_chp_',num2str(n_chp),'(t)  - P_chp_',num2str(n_chp),'(t-1) - af_chp_',num2str(n_chp),'(t)*Wdf_max(t) + af_chp_',num2str(n_chp),...
                '(t-1)*Wdf_min(t-1) >= -down_chp_',num2str(n_chp),'];']);
            end
    end
    %capacity lim of BS
    if t == 1
        for n_bs = 1 : N_bs
        eval(['con_Robust = con_Robust + [ E_bs_',num2str(n_bs),'(t)  == 0];']);
        eval(['con_Robust = con_Robust + [ af_bs_',num2str(n_bs),'(t)  == 0];']);
        eval(['con_Robust = con_Robust + [ af_bsd_',num2str(n_bs),'(t)  == 0] + [ P_bsd_',num2str(n_bs),'(t)  == 0];']);
        end 
    end
    if t>=2
        for n_bs = 1 : N_bs
        eval(['con_Robust = con_Robust + [ E_bs_',num2str(n_bs),'(t) ==  E_bs_',num2str(n_bs),'(t-1) + eta_BS*P_bsc_',num2str(n_bs),'(t)*3600 - 3600*P_bsd_',num2str(n_bs),'(t)/eta_BS];']);
        eval(['con_Robust = con_Robust + [ E_bs_',num2str(n_bs),'(t) <= E_BS_max(t)];']);
        eval(['con_Robust = con_Robust + [ E_bs_',num2str(n_bs),'(t) >= 0];']);
        eval(['con_Robust = con_Robust + [ af_bs_',num2str(n_bs),'(t) ==  af_bs_',num2str(n_bs),'(t-1) + eta_BS*af_bsc_',num2str(n_bs),'(t) - af_bsd_',num2str(n_bs),'(t)/eta_BS];']);        
        end 
    end         
    %gas balance lim
    con_Robust = con_Robust + [G_L(t) <= eta_P2G*P2Gf*P_p2g(t) + eta_P2G*P2Gf*af_p2g(t)*Wdf_min(t) - ( a_GT_1*P_gt_1(t)+b_GT_1)- ( a_GT_2*P_gt_2(t)+b_GT_2)...
        + a_GT_1*af_gt_1(t)*Wdf_min(t) + a_GT_2*af_gt_2(t)*Wdf_min(t) + G_buy(t) <= 1.1*G_L(t)];
    
    %power lim of P2G
    con_Robust = con_Robust + [P_p2g(t) + af_p2g(t)*Wdf_max(t) <= Pmax_P2G(t)];
    con_Robust = con_Robust + [P_p2g(t) + af_p2g(t)*Wdf_min(t) >= Pmin_P2G(t)];
end
%% =================================Define Chance Constraints==========================================%%
con_Chance = [];con_Test = [];
non_bus = n_bus - 1;%number of non-root nodes

    %Defining auxiliary variables of DRJCCs↓↓↓ for CHP,EB,GT and so on
    r_chp = 2*N_chp;
    delta_chp = sdpvar(n_T,r_chp);eposi_chp = sdpvar(n_T,r_chp);
    gamma_chp = sdpvar(N_Sub,n_T,r_chp);Zpos_chp = sdpvar(n_T,r_chp);
    Zneg_chp = sdpvar(n_T,r_chp);    
    A_chp = blkdiag([0 -1;0 1],[0 -1;0 1],[0 -1;0 1]);
    B_chp = blkdiag([1 0;-1 0],[1 0;-1 0],[1 0;-1 0]);
    b_chp =[];zeta_chp = [];
    for n_chp = 1:N_chp
    eval(['b_chp = [b_chp -Pmax_chp_',num2str(n_chp),'(1) Pmin_chp_',num2str(n_chp),'(1)];']);
    eval(['zeta_chp = [zeta_chp P_chp_',num2str(n_chp),' af_chp_',num2str(n_chp),'];']);
    end
    deltaQ_chp = sdpvar(n_T,r_chp);eposiQ_chp = sdpvar(n_T,r_chp);
    gammaQ_chp = sdpvar(N_Sub,n_T,r_chp);ZposQ_chp = sdpvar(N_Sub,n_T,r_chp);
    ZnegQ_chp = sdpvar(N_Sub,n_T,r_chp);    
    AQ_chp = blkdiag([0 -1;0 1],[0 -1;0 1],[0 -1;0 1]);
    BQ_chp = blkdiag([1 0;-1 0],[1 0;-1 0],[1 0;-1 0]);
    bQ_chp =[];zetaQ_chp = [];
    for n_chp = 1:N_chp
    eval(['bQ_chp = [bQ_chp -Qmax_chp_',num2str(n_chp),'(1) Qmin_chp_',num2str(n_chp),'(1)];']);
    eval(['zetaQ_chp = [zetaQ_chp Q_chp_',num2str(n_chp),' afQ_chp_',num2str(n_chp),'];']);
    end
    r_bs = 4*N_bs;
    delta_bs = sdpvar(n_T,r_bs);eposi_bs = sdpvar(n_T,r_bs);
    gamma_bs = sdpvar(N_Sub,n_T,r_bs);Zpos_bs = sdpvar(N_Sub,n_T,r_bs);
    Zneg_bs = sdpvar(N_Sub,n_T,r_bs);    
    A_bs = blkdiag([0 1;0 -1],[0 -1;0 1],[0 1;0 -1],[0 -1;0 1]);
    B_bs = blkdiag([1 0;-1 0],[1 0;-1 0],[1 0;-1 0],[1 0;-1 0]);
    b_bs =[];zeta_bs = [];
    for n_bs = 1:N_bs
    b_bs = [b_bs -Pmax_BS_C(1) Pmin_BS_C(1) -Pmax_BS_D(1) Pmin_BS_D(1)];
    eval(['zeta_bs = [zeta_bs P_bsc_',num2str(n_bs),' af_bsc_',num2str(n_bs),' P_bsd_',num2str(n_bs),' af_bsd_',num2str(n_bs),'];']);
    end    
    r_eb = 2*N_eb;
    delta_eb = sdpvar(n_T,r_eb);eposi_eb = sdpvar(n_T,r_eb);
    gamma_eb = sdpvar(N_Sub,n_T,r_eb);Zpos_eb = sdpvar(N_Sub,n_T,r_eb);
    Zneg_eb = sdpvar(N_Sub,n_T,r_eb);    
    A_eb = blkdiag([0 1;0 -1],[0 1;0 -1]);
    B_eb = blkdiag([1 0;-1 0],[1 0;-1 0]);
    b_eb =[];zeta_eb = [];
    for n_eb = 1:N_eb
    b_eb = [b_eb -Pmax_EB(1) Pmin_EB(1)];
    eval(['zeta_eb = [zeta_eb P_eb_',num2str(n_eb),' af_eb_',num2str(n_eb),'];']);
    end    
    r_gt = 2*N_gt;r_gtR = N_gt;
    delta_gt = sdpvar(n_T,r_gt);eposi_gt = sdpvar(n_T,r_gt);
    gamma_gt = sdpvar(N_Sub,n_T,r_gt);Zpos_gt = sdpvar(N_Sub,n_T,r_gt);
    rno = (9.863e-2)*(1e6)^(rho);
    Zneg_gt = sdpvar(N_Sub,n_T,r_gt);    
    A_gt = blkdiag([0 -1;0 1],[0 -1;0 1]);
    B_gt = blkdiag([1 0;-1 0],[1 0;-1 0]);
    b_gt =[]; zeta_gt = [];
    for n_gt = 1:N_gt
    eval(['b_gt = [b_gt -Pmax_GT_',num2str(n_gt),'(1) Pmin_GT_',num2str(n_gt),'(1)];']);
    eval(['zeta_gt = [zeta_gt P_gt_',num2str(n_gt),' af_gt_',num2str(n_gt),'];']);
    end
    delta_gtR = sdpvar(n_T,r_gt);eposi_gtR = sdpvar(n_T,r_gt);
    gamma_gtR = sdpvar(N_Sub,n_T,r_gt);Zpos_gtR = sdpvar(N_Sub,n_T,r_gt);
    Zneg_gtR = sdpvar(N_Sub,n_T,r_gt);    
    A_gtR = blkdiag([-1 0;0 1]);
    b_gtR =[];zeta_gtR = [];
    for n_gt = 1:N_gt
    eval(['b_gtR = [b_gtR Pmax_GT_',num2str(n_gt),'(1)*0.15 Pmax_GT_',num2str(n_gt),'(1)*0.15];']);
    eval(['zeta_gtR = [zeta_gtR af_gt_',num2str(n_gt),'];']);
    end
    r_car = 1;
    delta_car = sdpvar(n_T,r_car);
    gamma_car = sdpvar(N_Sub,n_T,r_car);
    eposi_car = sdpvar(n_T,r_car);
    Zpos_car = sdpvar(N_Sub,n_T,r_car);
    Zneg_car = sdpvar(N_Sub,n_T,r_car);
    zeta_car = [P_chp_1 af_chp_1 P_chp_2 af_chp_2 P_chp_3 af_chp_3  P_gt_1 af_gt_1 P_gt_2 af_gt_2  P_p2g af_p2g P_buy];
    rf =0.55*ones(n_T,1);
    A_car = [0 -e_chp 0 -e_chp 0 -e_chp 0 -e_gt 0 -e_gt 0 -e_p2g 0];
    B_car = [e_chp 0 e_chp 0 e_chp 0 e_gt 0 e_gt 0 -e_p2g 0 e_grid];
    b_car = [-(1-lambda_e)*e_max*1.01];
    r_vol = 2;
    delta_vol = sdpvar(n_T,r_vol);
    gamma_vol = sdpvar(n_T,N_Sub,r_vol);
    eposi_vol = sdpvar(n_T,r_vol);
    Zpos_vol = sdpvar(N_Sub,n_T,r_vol);
    Zneg_vol = sdpvar(N_Sub,n_T,r_vol);
    
    %Chance cons of CHPs ↓↓↓↓↓ POWER LIM   
    for t = 1:n_T
        tran_xi = repmat(risk_factor,n_T,1)/r_chp;
        for r = 1:r_chp
        con_Chance = con_Chance + [delta_chp(t,r)*tran_xi(t) + eposi_chp(t,r)*rho + 1e-9 + sum(gamma_chp(:,t,r))/N_Sub <= 0];
        con_Chance = con_Chance + [Zpos_chp(t,r) >= 0] + [Zneg_chp(t,r) >= 0] + [eposi_chp(t,r) >= 0];
        con_Chance = con_Chance + [eposi_chp(t,r) >= norm(A_chp(r,:)*zeta_chp(t,:)'-Zpos_chp(t,r)+Zneg_chp(t,r),inf)];
        con_Chance = con_Chance + [eposi_chp(t,r) >= norm(-Zpos_chp(t,r)+Zneg_chp(t,r),inf)]; 
        con_Chance = con_Chance + [gamma_chp(:,t,r) >= repmat((B_chp(r,:)*zeta_chp(t,:)'+ b_chp(r)) - delta_chp(t,r)+ Zpos_chp(t,r)*Wdf_max(t) - Zneg_chp(t,r)*Wdf_min(t),N_Sub,1)...
                 + diag(Wdf_S(:,t))*( repmat(A_chp(r,:)*zeta_chp(t,:)' - Zpos_chp(t,r) + Zneg_chp(t,r),N_Sub,1) )];
        con_Chance = con_Chance + [gamma_chp(:,t,r) >= repmat(Zpos_chp(t,r)*Wdf_max(t) - Zneg_chp(t,r)*Wdf_min(t),N_Sub,1)...
                 + diag(Wdf_S(:,t))*(repmat(- Zpos_chp(t,r) + Zneg_chp(t,r),N_Sub,1) )];
        end
    
        for r = 1:r_chp
            con_Chance = con_Chance + [deltaQ_chp(t,r)*tran_xi(t) + eposiQ_chp(t,r)*rho+ 1e-9 + sum(gammaQ_chp(:,t,r))/N_Sub <= 0];
            con_Chance = con_Chance + [ZposQ_chp(:,t,r) >= 0] + [ZnegQ_chp(:,t,r) >= 0] + [eposiQ_chp(t,r) >= 0];
            con_Chance = con_Chance + [repmat(eposiQ_chp(t,r),N_Sub,1) >= -ZposQ_chp(:,t,r)+ZnegQ_chp(:,t,r)];
            con_Chance = con_Chance + [repmat(eposiQ_chp(t,r),N_Sub,1) >= -ZposQ_chp(:,t,r)+ZnegQ_chp(:,t,r)]; 
            con_Chance = con_Chance + [gammaQ_chp(:,t,r) >= repmat((BQ_chp(r,:)*zetaQ_chp(t,:)'+ bQ_chp(r)) - deltaQ_chp(t,r),N_Sub,1) + ZposQ_chp(:,t,r)*Wdf_max(t) - ZnegQ_chp(:,t,r)*Wdf_min(t)...
                 + diag(Wdf_S(:,t))*( repmat(AQ_chp(r,:)*zetaQ_chp(t,:)',N_Sub,1) - ZposQ_chp(:,t,r) + ZnegQ_chp(:,t,r))];
            con_Chance = con_Chance + [gammaQ_chp(:,t,r) >= ZposQ_chp(:,t,r)*Wdf_max(t) - ZnegQ_chp(:,t,r)*Wdf_min(t)...
                 + diag(Wdf_S(:,t))*( - ZposQ_chp(:,t,r) + ZnegQ_chp(:,t,r))];
        end

        %Chance cons of BSs ↓↓↓↓↓POWER LIM    
        tran_xi = repmat(risk_factor,n_T,1)/r_bs;
        for r = 1:r_bs
        con_Chance = con_Chance + [delta_bs(t,r)*tran_xi(t) + eposi_bs(t,r)*rho+ 1e-9 + sum(gamma_bs(:,t,r))/N_Sub <= 0];
        con_Chance = con_Chance + [Zpos_bs(t,r) >= 0] + [Zneg_bs(t,r) >= 0] + [eposi_bs(t,r) >= 0];
        con_Chance = con_Chance + [eposi_bs(t,r) >= norm(A_bs(r,:)*zeta_bs(t,:)'-Zpos_bs(t,r)+Zneg_bs(t,r),inf)];
        con_Chance = con_Chance + [eposi_bs(t,r) >= norm(-Zpos_bs(t,r)+Zneg_bs(t,r),inf)]; 
        con_Chance = con_Chance + [gamma_bs(:,t,r) >= repmat((B_bs(r,:)*zeta_bs(t,:)'+ b_bs(r)) - delta_bs(t,r)+ Zpos_bs(t,r)*Wdf_max(t) - Zneg_bs(t,r)*Wdf_min(t),N_Sub,1)...
                 + diag(Wdf_S(:,t))*( repmat(A_bs(r,:)*zeta_bs(t,:)' - Zpos_bs(t,r) + Zneg_bs(t,r),N_Sub,1) )];
        con_Chance = con_Chance + [gamma_bs(:,t,r) >= repmat(Zpos_bs(t,r)*Wdf_max(t) - Zneg_bs(t,r)*Wdf_min(t),N_Sub,1)...
                 + diag(Wdf_S(:,t))*(repmat(- Zpos_bs(t,r) + Zneg_bs(t,r),N_Sub,1) )];
        end

        %Chance cons of EBs ↓↓↓↓↓POWER LIM
        tran_xi = repmat(risk_factor,n_T,1)/r_eb;
        for r = 1:r_eb
        con_Chance = con_Chance + [delta_eb(t,r)*tran_xi(t) + eposi_eb(t,r)*rho+ 1e-9 + sum(gamma_eb(:,t,r))/N_Sub <= 0];
        con_Chance = con_Chance + [Zpos_eb(t,r) >= 0] + [Zneg_eb(t,r) >= 0] + [eposi_eb(t,r) >= 0];
        con_Chance = con_Chance + [eposi_eb(t,r) >= norm(A_eb(r,:)*zeta_eb(t,:)'-Zpos_eb(t,r)+Zneg_eb(t,r),inf)];
        con_Chance = con_Chance + [eposi_eb(t,r) >= norm(-Zpos_eb(t,r)+Zneg_eb(t,r),inf)]; 
        con_Chance = con_Chance + [gamma_eb(:,t,r) >= repmat((B_eb(r,:)*zeta_eb(t,:)'+ b_eb(r)) - delta_eb(t,r)+ Zpos_eb(t,r)*Wdf_max(t) - Zneg_eb(t,r)*Wdf_min(t),N_Sub,1)...
                 + diag(Wdf_S(:,t))*( repmat(A_eb(r,:)*zeta_eb(t,:)' - Zpos_eb(t,r) + Zneg_eb(t,r),N_Sub,1) )];
        con_Chance = con_Chance + [gamma_eb(:,t,r) >= repmat(Zpos_eb(t,r)*Wdf_max(t) - Zneg_eb(t,r)*Wdf_min(t),N_Sub,1)...
                 + diag(Wdf_S(:,t))*(repmat(- Zpos_eb(t,r) + Zneg_eb(t,r),N_Sub,1) )];
        end

        %Chance cons of GTs ↓↓↓↓↓POWER LIM
        tran_xi = repmat(risk_factor,n_T,1)/r_gt;
        for r = 1:r_gt
        con_Chance = con_Chance + [delta_gt(t,r)*tran_xi(t) + eposi_gt(t,r)*rho+ 1e-9 + sum(gamma_gt(:,t,r))/N_Sub <= 0];
        con_Chance = con_Chance + [Zpos_gt(t,r) >= 0] + [Zneg_gt(t,r) >= 0] + [eposi_gt(t,r) >= 0];
        con_Chance = con_Chance + [eposi_gt(t,r) >= norm(A_gt(r,:)*zeta_gt(t,:)'-Zpos_gt(t,r)+Zneg_gt(t,r),inf)];
        con_Chance = con_Chance + [eposi_gt(t,r) >= norm(-Zpos_gt(t,r)+Zneg_gt(t,r),inf)]; 
        con_Chance = con_Chance + [gamma_gt(:,t,r) >= repmat((B_gt(r,:)*zeta_gt(t,:)'+ b_gt(r)) - delta_gt(t,r)+ Zpos_gt(t,r)*Wdf_max(t) - Zneg_gt(t,r)*Wdf_min(t),N_Sub,1)...
                 + diag(Wdf_S(:,t))*( repmat(A_gt(r,:)*zeta_gt(t,:)' - Zpos_gt(t,r) + Zneg_gt(t,r),N_Sub,1) )];
        con_Chance = con_Chance + [gamma_gt(:,t,r) >= repmat(Zpos_gt(t,r)*Wdf_max(t) - Zneg_gt(t,r)*Wdf_min(t),N_Sub,1)...
                 + diag(Wdf_S(:,t))*(repmat(- Zpos_gt(t,r) + Zneg_gt(t,r),N_Sub,1) )];
        
        end
        
        for r = 1:r_gtR
        con_Chance = con_Chance + [delta_gtR(t,r)*tran_xi(t) + eposi_gtR(t,r)*rho+ 1e-9 + sum(gamma_gtR(:,t,r))/N_Sub <= 0];
        con_Chance = con_Chance + [Zpos_gtR(t,r) >= 0] + [Zneg_gtR(t,r) >= 0] + [eposi_gtR(t,r) >= 0];
        con_Chance = con_Chance + [eposi_gtR(t,r) >= norm(A_gtR(r,:)*zeta_gtR(t,:)'-Zpos_gtR(t,r)+Zneg_gtR(t,r),inf)];
        con_Chance = con_Chance + [eposi_gtR(t,r) >= norm(-Zpos_gtR(t,r)+Zneg_gtR(t,r),inf)]; 
        con_Chance = con_Chance + [gamma_gtR(:,t,r) >= repmat(-b_gtR(r) - delta_gtR(t,r)+ Zpos_gtR(t,r)*Wdf_max(t) - Zneg_gtR(t,r)*Wdf_min(t),N_Sub,1)...
                 + diag(Wdf_S(:,t))*( repmat(A_gtR(r,:)*zeta_gtR(t,:)' - Zpos_gtR(t,r) + Zneg_gtR(t,r),N_Sub,1) )];
        con_Chance = con_Chance + [gamma_gtR(:,t,r) >= repmat(Zpos_gtR(t,r)*Wdf_max(t) - Zneg_gtR(t,r)*Wdf_min(t),N_Sub,1)...
                 + diag(Wdf_S(:,t))*(repmat(- Zpos_gtR(t,r) + Zneg_gtR(t,r),N_Sub,1) )];
        
        end
        %Chance cons of Carbon Emissions ↓↓↓↓↓EMISSION LIM
        tran_xi = repmat(risk_factor,n_T,1)/r_car;
        for r = 1:r_car
            con_Chance = con_Chance + [delta_car(t,r)*tran_xi(t) + eposi_car(t,r)*rho^rf(t)+ 1e-9 + sum(gamma_car(:,t,r))/N_Sub <= 0];
            con_Chance = con_Chance + [Zpos_car(t,r) >= 0] + [Zneg_car(t,r) >= 0] + [eposi_car(t,r) >= 0];
            con_Chance = con_Chance + [eposi_car(t,r) >= norm(A_car(r,:)*zeta_car(t,:)'-Zpos_car(t,r)+Zneg_car(t,r),inf)];
            con_Chance = con_Chance + [eposi_car(t,r) >= norm(-Zpos_car(t,r)+Zneg_car(t,r),inf)]; 
            con_Chance = con_Chance + [gamma_car(:,t,r) >= repmat((B_car(r,:)*zeta_car(t,:)'+ b_car(r)) - delta_car(t,r)+ Zpos_car(t,r)*Wdf_max(t) - Zneg_car(t,r)*Wdf_min(t),N_Sub,1)...
                 + diag(Wdf_S(:,t))*( repmat(A_car(r,:)*zeta_car(t,:)' - Zpos_car(t,r) + Zneg_car(t,r),N_Sub,1) )];
            con_Chance = con_Chance + [gamma_car(:,t,r) >= repmat(Zpos_car(t,r)*Wdf_max(t) - Zneg_car(t,r)*Wdf_min(t),N_Sub,1)...
                 + diag(Wdf_S(:,t))*(repmat(- Zpos_car(t,r) + Zneg_car(t,r),N_Sub,1) )];
        end
    end   
     
for t = 1:n_T  
    %Chance cons of Nodal Voltage ↓↓↓↓↓LIM +-6%
    af_dev_Vp = [0 0 0 0 0 0 0 0 af_gt_1(t) 0 -af_eb_1(t) 0 0 (af_bsd_1(t)-af_bsc_1(t)) 0 0 af_chp_1(t) 0 -af_p2g(t) 0 af_chp_2(t) 0 0 af_chp_3(t) 0 (af_bsd_2(t)-af_bsc_2(t)) 0 -af_eb_2(t) 0 0 0 af_gt_2(t)];%决策向量堆叠
    af_dev_Vn = [0 0 0 0 0 0 0 0 af_gt_1(t) 0 -af_eb_1(t) 0 0 (af_bsd_1(t)-af_bsc_1(t)) 0 0 af_chp_1(t) 0 -af_p2g(t) 0 af_chp_2(t) 0 0 af_chp_3(t) 0 (af_bsd_2(t)-af_bsc_2(t)) 0 -af_eb_2(t) 0 0 0 af_gt_2(t)];%决策向量堆叠
    afQ_dev = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 afQ_chp_1(t) 0 0 0 afQ_chp_2(t) 0 0 afQ_chp_3(t) 0 0 0 0 0 0 0 0];%无功决策向量堆叠
    con_Test = con_Test + [af_nbus_VP(t,:) == af_dev_Vp] + [af_nbus_VN(t,:) == af_dev_Vn] + [afQ_nbus(t,:) == afQ_dev];

    con_Test = con_Test + [ones(11, 1)*0.94*V_0 <= V_bus(t,2:12)' - Bra_B(1:11,:)*((R*Bra_B*Conf_WT'+ X*Bra_B*Conf_WT'*0.75)*Wdf_max(t)*1e3)/V_0...
         - Bra_B(1:11,:)*Wdf_max(t)*1e3*(- R*Bra_B*af_nbus_VP(t,:)' - X*Bra_B*afQ_nbus(t,:)')/V_0 <= ones(11, 1)*1.06*V_0];
    con_Test = con_Test + [ones(11, 1)*0.94*V_0 <= V_bus(t,2:12)' - Bra_B(1:11,:)*((R*Bra_B*Conf_WT'+ X*Bra_B*Conf_WT'*0.75)*Wdf_min(t)*1e3)/V_0...
         - Bra_B(1:11,:)*Wdf_min(t)*1e3*(- R*Bra_B*af_nbus_VP(t,:)' - X*Bra_B*afQ_nbus(t,:)')/V_0 <= ones(11, 1)*1.06*V_0];
    con_Test = con_Test + [ones(15, 1)*0.94*V_0 <= V_bus(t,19:33)' - Bra_B(18:32,:)*((R*Bra_B*Conf_WT'+ X*Bra_B*Conf_WT'*0.75)*Wdf_max(t)*1e3)/V_0...
         - Bra_B(18:32,:)*Wdf_max(t)*1e3*(- R*Bra_B*af_nbus_VP(t,:)' - X*Bra_B*afQ_nbus(t,:)')/V_0 <= ones(15, 1)*1.06*V_0];
    con_Test = con_Test + [ones(15, 1)*0.94*V_0 <= V_bus(t,19:33)' - Bra_B(18:32,:)*((R*Bra_B*Conf_WT'+ X*Bra_B*Conf_WT'*0.75)*Wdf_min(t)*1e3)/V_0...
         - Bra_B(18:32,:)*Wdf_min(t)*1e3*(- R*Bra_B*af_nbus_VP(t,:)' - X*Bra_B*afQ_nbus(t,:)')/V_0 <= ones(15, 1)*1.06*V_0]; 
    tran_xi = repmat(risk_factor,n_T,1)/r_vol;   
     
    for n = 12:17
          con_Chance = con_Chance + [delta_vol(t,1)*tran_xi(t) + eposi_vol(t,1)*rno*1.5e3 + sum(gamma_vol(t,:,1))/N_Sub <= 0];
          con_Chance = con_Chance + [Zpos_vol(t,1) >= 0] + [Zneg_vol(t,1) >= 0] + [eposi_vol(t,1) >= 0];
          con_Chance = con_Chance + [eposi_vol(t,1) >= norm(- (Bra_B(n,:)/V_0)*((R*Bra_B*Conf_WT'+ X*Bra_B*Conf_WT'*0.75)...
         + Bra_B*(- R*Bra_B*af_nbus_VP(t,:)' - X*Bra_B*afQ_nbus(t,:)'))-Zpos_vol(t,1)+Zneg_vol(t,1),inf)];
          con_Chance = con_Chance + [eposi_vol(t,1) >= norm(-Zpos_vol(t,1)+Zneg_vol(t,1),inf)]; 
          con_Chance = con_Chance + [gamma_vol(t,:,1)' >= repmat((V_bus(t,n+1) - 1.0602*V_0) - delta_vol(t,1)+ Zpos_vol(t,1)*Wdf_max(t)*1e3 - Zneg_vol(t,1)*Wdf_min(t)*1e3,N_Sub,1)...
                 + diag(Wdf_S(:,t)*1e3)*( repmat(- (Bra_B(n,:)/V_0)*((R*Bra_B*Conf_WT'+ X*Bra_B*Conf_WT'*0.75)...
         + Bra_B*(- R*Bra_B*af_nbus_VP(t,:)' - X*Bra_B*afQ_nbus(t,:)')) - Zpos_vol(t,1) + Zneg_vol(t,1),N_Sub,1))];
          con_Chance = con_Chance + [gamma_vol(t,:,1)' >= repmat(Zpos_vol(t,1)*Wdf_max(t)*1e3 - Zneg_vol(t,1)*Wdf_min(t)*1e3,N_Sub,1)...
                 + diag(Wdf_S(:,t)*1e3)*(repmat(- Zpos_vol(t,1) + Zneg_vol(t,1),N_Sub,1) )];
             
          con_Chance = con_Chance + [delta_vol(t,2)*tran_xi(t) + eposi_vol(t,2)*rno*1.5e3 + sum(gamma_vol(t,:,2))/N_Sub <= 0];
          con_Chance = con_Chance + [Zpos_vol(t,2) >= 0] + [Zneg_vol(t,2) >= 0] + [eposi_vol(t,2) >= 0];
          con_Chance = con_Chance + [eposi_vol(t,2) >= norm( (Bra_B(n,:)/V_0)*((R*Bra_B*Conf_WT'+ X*Bra_B*Conf_WT'*0.75)...
         + Bra_B*(- R*Bra_B*af_nbus_VP(t,:)' - X*Bra_B*afQ_nbus(t,:)'))-Zpos_vol(t,2)+Zneg_vol(t,2),inf)];
          con_Chance = con_Chance + [eposi_vol(t,2) >= norm(-Zpos_vol(t,2)+Zneg_vol(t,2),inf)]; 
          con_Chance = con_Chance + [gamma_vol(t,:,2)' >= repmat((-V_bus(t,n+1) + 0.94*V_0) - delta_vol(t,2)+ Zpos_vol(t,2)*Wdf_max(t)*1e3 - Zneg_vol(t,2)*Wdf_min(t)*1e3,N_Sub,1)...
                 + diag(Wdf_S(:,t)*1e3)*( repmat( (Bra_B(n,:)/V_0)*((R*Bra_B*Conf_WT'+ X*Bra_B*Conf_WT'*0.75)...
         + Bra_B*(- R*Bra_B*af_nbus_VP(t,:)' - X*Bra_B*afQ_nbus(t,:)'))  - Zpos_vol(t,2) + Zneg_vol(t,2),N_Sub,1) )];
          con_Chance = con_Chance + [gamma_vol(t,:,2)' >= repmat(Zpos_vol(t,2)*Wdf_max(t)*1e3 - Zneg_vol(t,2)*Wdf_min(t)*1e3,N_Sub,1)...
                 + diag(Wdf_S(:,t)*1e3)*(repmat(- Zpos_vol(t,2) + Zneg_vol(t,2),N_Sub,1) )];
    end
     end
%% ===============================Define objective functions=============================== %%

%nominal sate
obj_MEMG = C_buy_E'*P_buy - C_sell_E(1:n_T)'*P_sell + C_chp_1'*P_chp_1 + C_chp_2'*P_chp_2 + C_chp_3'*P_chp_3...
    - C_BS_ch'*(P_bsc_1 + P_bsc_2) + C_BS_dch'*(P_bsd_1 + P_bsd_2) + C_P2G'*P_p2g + C_EB'*(P_eb_1 + P_eb_2) + C_Gas'*G_buy;%发电成本+购气成本
obj_MEMG = obj_MEMG + c_e_car'*(e_chp*(P_chp_1 + P_chp_2 + P_chp_3) + e_gt*(P_gt_1 + P_gt_2) + e_gas*G_L(1:n_T) + e_grid*P_buy - e_p2g*P_p2g); %碳成本

%under uncertainty
obj_MEMG =obj_MEMG + C_P2G'*(diag(Wdf_max(1:n_T))*af_p2g);
obj_MEMG = obj_MEMG + c_e_car(1:n_T)'*(-e_chp*(diag(Wdf_down(1:n_T))*af_chp_1 + diag(Wdf_down(1:n_T))*af_chp_2 + diag(Wdf_down(1:n_T))*af_chp_3)...
    - e_gt*diag(Wdf_down(1:n_T))*(af_gt_1+af_gt_2) -e_p2g*diag(Wdf_down(1:n_T))*af_p2g);

%Defining auxiliary variables of Objctives↓↓↓
con_Obj = [];
for n_chp = 1: N_chp
eval(['alpha_chp_',num2str(n_chp),'= sdpvar(n_T,1);']);
eval(['beta_chp_',num2str(n_chp),'= sdpvar(N_Sub,1);']);
end
for n_bs = 1: N_bs
eval(['alpha_bsc_',num2str(n_bs),'= sdpvar(n_T,1);']);
eval(['beta_bsc_',num2str(n_bs),'= sdpvar(N_Sub,1);']);
eval(['C_bsc_',num2str(n_bs),'= C_BS_ch;']);
eval(['alpha_bsd_',num2str(n_bs),'= sdpvar(n_T,1);']);
eval(['beta_bsd_',num2str(n_bs),'= sdpvar(N_Sub,1);']);
eval(['C_bsd_',num2str(n_bs),'= C_BS_dch;']);
end
for n_eb = 1: N_eb
eval(['alpha_eb_',num2str(n_eb),'= sdpvar(n_T,1);']);
eval(['beta_eb_',num2str(n_eb),'= sdpvar(N_Sub,1);']);
eval(['C_eb_',num2str(n_eb),'= C_EB;']);
end
for  t = 1:n_T
    for n_chp = 1: N_chp
    %Transformation objective function for CHP ↓↓↓
    eval(['obj_MEMG = obj_MEMG + alpha_chp_',num2str(n_chp),'(t)*rho+ 1e-9 + sum(beta_chp_',num2str(n_chp),')/N_Sub;']);
    eval(['con_Obj =  con_Obj +[beta_chp_',num2str(n_chp),'>= repmat(-C_chp_',num2str(n_chp),'(t)*af_chp_',num2str(n_chp),...
        '(t)*Wdf_min(t) + alpha_chp_',num2str(n_chp),'(t)*Wdf_min(t),N_Sub,1) - diag(Wdf_S(:,t))*repmat(alpha_chp_',num2str(n_chp),'(t),N_Sub,1)];']);
    eval(['con_Obj =  con_Obj +[beta_chp_',num2str(n_chp),'>= repmat(-C_chp_',num2str(n_chp),'(t)*af_chp_',num2str(n_chp),...
        '(t)*Wdf_max(t) + alpha_chp_',num2str(n_chp),'(t)*Wdf_max(t),N_Sub,1) - diag(Wdf_S(:,t))*repmat(alpha_chp_',num2str(n_chp),'(t),N_Sub,1) ];']);
    eval(['con_Obj =  con_Obj +[beta_chp_',num2str(n_chp),'>= diag(Wdf_S(:,t))*repmat(-C_chp_',num2str(n_chp),'(t)*af_chp_',num2str(n_chp),...
        '(t),N_Sub,1)];']);
    eval(['con_Obj =  con_Obj + [alpha_chp_',num2str(n_chp),'(t) >= 0];']);
    end
    %Transformation objective function for BS ↓↓↓
    for n_bs = 1: N_bs
    eval(['obj_MEMG = obj_MEMG + alpha_bsc_',num2str(n_bs),'(t)*rho+ 1e-9 + sum(beta_bsc_',num2str(n_bs),')/N_Sub;']);
    eval(['con_Obj =  con_Obj +[beta_bsc_',num2str(n_bs),'>= repmat(-C_bsc_',num2str(n_bs),'(t)*af_bsc_',num2str(n_bs),...
        '(t)*Wdf_min(t) + alpha_bsc_',num2str(n_bs),'(t)*Wdf_min(t),N_Sub,1) - diag(Wdf_S(:,t))*repmat(alpha_bsc_',num2str(n_bs),'(t),N_Sub,1)];']);
    eval(['con_Obj =  con_Obj +[beta_bsc_',num2str(n_bs),'>= repmat(-C_bsc_',num2str(n_bs),'(t)*af_bsc_',num2str(n_bs),...
        '(t)*Wdf_max(t) + alpha_bsc_',num2str(n_bs),'(t)*Wdf_max(t),N_Sub,1) - diag(Wdf_S(:,t))*repmat(alpha_bsc_',num2str(n_bs),'(t),N_Sub,1) ];']);
    eval(['con_Obj =  con_Obj +[beta_bsc_',num2str(n_bs),'>= diag(Wdf_S(:,t))*repmat(-C_bsc_',num2str(n_bs),'(t)*af_bsc_',num2str(n_bs),...
        '(t),N_Sub,1)];']);
    eval(['con_Obj =  con_Obj + [alpha_bsc_',num2str(n_bs),'(t) >= 0];']);    
    end
    for n_bs = 1: N_bs
    eval(['obj_MEMG = obj_MEMG + alpha_bsd_',num2str(n_bs),'(t)*rho+ 1e-9 + sum(beta_bsd_',num2str(n_bs),')/N_Sub;']);
    eval(['con_Obj =  con_Obj +[beta_bsd_',num2str(n_bs),'>= repmat(-C_bsd_',num2str(n_bs),'(t)*af_bsd_',num2str(n_bs),...
        '(t)*Wdf_min(t) + alpha_bsd_',num2str(n_bs),'(t)*Wdf_min(t),N_Sub,1) - diag(Wdf_S(:,t))*repmat(alpha_bsd_',num2str(n_bs),'(t),N_Sub,1)];']);
    eval(['con_Obj =  con_Obj +[beta_bsd_',num2str(n_bs),'>= repmat(-C_bsd_',num2str(n_bs),'(t)*af_bsd_',num2str(n_bs),...
        '(t)*Wdf_max(t) + alpha_bsd_',num2str(n_bs),'(t)*Wdf_max(t),N_Sub,1) - diag(Wdf_S(:,t))*repmat(alpha_bsd_',num2str(n_bs),'(t),N_Sub,1) ];']);
    eval(['con_Obj =  con_Obj +[beta_bsd_',num2str(n_bs),'>= diag(Wdf_S(:,t))*repmat(-C_bsd_',num2str(n_bs),'(t)*af_bsd_',num2str(n_bs),...
        '(t),N_Sub,1)];']);
    eval(['con_Obj =  con_Obj + [alpha_bsd_',num2str(n_bs),'(t) >= 0];']);  
    end
    %Transformation objective function for EB ↓↓↓
    for n_eb = 1: N_eb
    eval(['obj_MEMG = obj_MEMG + alpha_eb_',num2str(n_eb),'(t)*rho+ 1e-9 + sum(beta_eb_',num2str(n_eb),')/N_Sub;']);
    eval(['con_Obj =  con_Obj +[beta_eb_',num2str(n_eb),'>= repmat(C_eb_',num2str(n_eb),'(t)*af_eb_',num2str(n_eb),...
        '(t)*Wdf_min(t) + alpha_eb_',num2str(n_eb),'(t)*Wdf_min(t),N_Sub,1) - diag(Wdf_S(:,t))*repmat(alpha_eb_',num2str(n_eb),'(t),N_Sub,1)];']);
    eval(['con_Obj =  con_Obj +[beta_eb_',num2str(n_eb),'>= repmat(C_eb_',num2str(n_eb),'(t)*af_eb_',num2str(n_eb),...
        '(t)*Wdf_max(t) + alpha_eb_',num2str(n_eb),'(t)*Wdf_max(t),N_Sub,1) - diag(Wdf_S(:,t))*repmat(alpha_eb_',num2str(n_eb),'(t),N_Sub,1) ];']);
    eval(['con_Obj =  con_Obj +[beta_eb_',num2str(n_eb),'>= diag(Wdf_S(:,t))*repmat(-C_eb_',num2str(n_eb),'(t)*af_eb_',num2str(n_eb),...
        '(t),N_Sub,1)];']);
    eval(['con_Obj =  con_Obj + [alpha_eb_',num2str(n_eb),'(t) >= 0];']);
    end
end



constraints_all = con_Chance + con_Robust + con_General + con_Test + con_Obj;%integrate all MEMG constraints
disp("Modeling complete");
tic%start of timer
op=sdpsettings('solver','gurobi','savesolveroutput',1,'verbose',0);%set solver
disp("========================Calculation in process======================")

sol_all = optimize(constraints_all,obj_MEMG,op);
if sol_all.problem == 0
disp("Solved successfully!");
toc%end of timer
op_sol = sprintf('The minimized operation cost of MEMG is %d (W-DRO-method).',sol_all.solveroutput.result.objval);
disp(op_sol);
elseif sol_all.problem ~= 0
    [model,~] =  export(constraints_all, obj_MEMG,sdpsettings('solver','gurobi'));
    iis = gurobi_iis(model);
    gurobi_write(model , 'MyFile.lp');
    find(iis.Arows)
    msg = 'Model infeasible!';
    error(msg)
end


%Save output
af_MEMG = [value(af_chp_1)+value(af_chp_2)+value(af_chp_3) value(af_gt_1)+value(af_gt_2) value(af_eb_1)+value(af_eb_2) value(af_bsc_1)+value(af_bsc_2) value(af_bsd_1)+value(af_bsd_2)  value(af_p2g)];
xE_MEMG = [value(P_chp_1)+value(P_chp_2)+value(P_chp_3) value(P_gt_1)+value(P_gt_2) (-value(P_eb_1)-value(P_eb_2)) (value(P_bsd_1)-value(P_bsc_1)+value(P_bsd_2)-value(P_bsc_2)) -value(P_p2g) (value(P_buy)-value(P_sell))];
VE_MEMG = value(V_bus)/V_0;%to 1 p.u.
PL_MEMG = [value(P_flow_in) value(P_flow_out)];%to MW