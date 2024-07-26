% DO method: scripted by Zehao Cao
yalmip('clear');
%% Loading power network data
mpc = loadcase('case33mg');%modified IEEE 33-BUS
node = mpc.bus;
P_D = repmat(node(:,3),1,n_T);P_D = P_D*diag(P_L(1:n_T));%active load matrix
Q_D = repmat(node(:,4),1,n_T);Q_D = Q_D*diag(P_L(1:n_T));%reactive load matrix
%% Define dicision variables
for n_chp = 1 : N_chp
eval(['P_chp_',num2str(n_chp),'= sdpvar(n_T,1);']);% P&Q of CHP
eval(['Q_chp_',num2str(n_chp),'= sdpvar(n_T,1);']);
end
for n_gt = 1 : N_gt
eval(['P_gt_',num2str(n_gt),'= sdpvar(n_T,1);']);
eval(['af_gt_',num2str(n_gt),'= sdpvar(n_T,1);']);
end
for n_eb = 1 : N_eb
eval(['P_eb_',num2str(n_eb),'= sdpvar(n_T,1);']);
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
P_buy = sdpvar(n_T,1);%elec purchase
P_sell = sdpvar(n_T,1);%elec sell
af_p2g = sdpvar(n_T,1);
G_buy = sdpvar(n_T,1);%gas purchase
P_flow_out = sdpvar(n_T,n_bus);Q_flow_out = sdpvar(n_T,n_bus);%powerflow of lines "in"
P_flow_in = sdpvar(n_T,n_bus);Q_flow_in = sdpvar(n_T,n_bus);%powerflow of lines "out"
V_bus = sdpvar(n_T,n_bus);%nodal voltage
%% ===============================Define General Constraints=============================== %%
con_General = [];
disp("Modeling MEMG as DO...");
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
    %active power flow lim
    Q_dev = [Q_chp_1(t); Q_chp_2(t); Q_chp_3(t);];%decision vectors stacking
    Q_DER = [0.75*PW_1(t); 0.75*PW_2(t); 0.75*PW_3(t)];%DERs' vectors stacking
    for n = 1:n_bus
        if n == 1
        con_General = con_General + [ Q_flow_in(t,n) == 0 ]+[ Q_flow_out(t,n) == 0];%power flow of root
        elseif n == 3||n == 19
        con_General = con_General + [Q_flow_in(t,3) + Q_flow_in(t,19) == Q_flow_out(t,2)]+...
            [Q_flow_in(t,n)+ Conf_Q(1:3,n)'*Q_dev + Conf_Q(4:6,n)'*Q_DER == Q_flow_out(t,n) + Q_D(n,t)];%%power flow of forked lines
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
        con_General = con_General + [V_bus(t,n) == V_bus(t,n-1)-(R(n-1,n-1)*P_flow_in(t,n)*1e3+X(n-1,n-1)*Q_flow_in(t,n)*1e3)/V_0];%voltage of forked nodes
        end
     end    
    %ensuer variables non-negative
    for n_chp = 1 : N_chp
    eval([' con_General = con_General + [P_chp_',num2str(n_chp),'(t) >= 0];']);
    end
    for n_gt = 1 : N_gt
    eval([' con_General = con_General + [P_gt_',num2str(n_gt),'(t) >= 0];']);
    eval([' con_General = con_General + [0 <= af_gt_',num2str(n_gt),'(t) <= 1];']);
    end
    for n_eb = 1 : N_eb
    eval([' con_General = con_General + [P_eb_',num2str(n_eb),'(t) >= 0];']);
    end
    for n_bs = 1 : N_bs
    eval([' con_General = con_General + [P_bsc_',num2str(n_bs),'(t) >= 0];']);
    eval([' con_General = con_General + [0 <= af_bsc_',num2str(n_bs),'(t) <= 1];']);
    eval([' con_General = con_General + [P_bsd_',num2str(n_bs),'(t) >= 0];']);
    eval([' con_General = con_General + [0 <= af_bsd_',num2str(n_bs),'(t) <= 1];']);
    eval([' con_General = con_General + [0 <= af_bs_',num2str(n_bs),'(t) <= 1];']);
    eval([' con_General = con_General + [af_bsc_',num2str(n_bs),'(t) * af_bsd_',num2str(n_bs),'(t) == 0];']);
    eval([' con_General = con_General + [E_bs_',num2str(n_bs),'(t) >= 0];']);
    end
    con_General = con_General + [P_p2g(t) >= 0] + [0 <= af_p2g(t) <= 1];
    con_General = con_General + [700>= G_buy(t) >= 0];
    con_General = con_General + [P_buy(t) >= 0] + [P_sell(t) >= 0];
    %LDR lim, all t to 1
    con_General = con_General + [af_gt_1(t) + af_gt_2(t) + af_p2g(t) + af_bsc_1(t) + af_bsc_2(t) + af_bsd_1(t) + af_bsd_2(t) == 1];
    con_General = con_General + [af_gt_1(t) == af_gt_2(t)] + [af_gt_2(t) == af_p2g(t)] + [af_p2g(t) == af_bsc_1(t) + af_bsd_1(t)] +...
        [ af_bsc_1(t) + af_bsd_1(t) == af_bsc_2(t) + af_bsd_2(t)] + [af_bsc_2(t) + af_bsd_2(t) == af_gt_1(t)];
    
    %heat balance constraints
    con_General = con_General + [H_E*(P_chp_1(t) + P_chp_2(t) + P_chp_3(t)) + eta_EB*(P_eb_1(t) + P_eb_2(t)) == H_L(t,t)];
end
%% ===============================Define General Constraints=============================== %%
con_Robust = [];
for t = 1:n_T
    %ramping cons of CHP
    if t>=2%from t=2
            for n_chp = 1 : N_chp
            eval([' con_Robust = con_Robust + [ P_chp_',num2str(n_chp),'(t)  - P_chp_',num2str(n_chp),'(t-1) <= ramp_chp_',num2str(n_chp),'];']);
            eval([' con_Robust = con_Robust + [ P_chp_',num2str(n_chp),'(t)  - P_chp_',num2str(n_chp),'(t-1) >= -down_chp_',num2str(n_chp),'];']);
            end
    end
    %capacity lim of BS
    if t == 1
        for n_bs = 1 : N_bs
        eval(['con_Robust = con_Robust + [ E_bs_',num2str(n_bs),'(t)  == 0];']);
        end 
    end
    if t>=2
        for n_bs = 1 : N_bs
        eval(['con_Robust = con_Robust + [ E_bs_',num2str(n_bs),'(t-1) + eta_BS*P_bsc_',num2str(n_bs),'(t)*3600 - 3600*P_bsd_',num2str(n_bs),'(t)/eta_BS  <= E_BS_max(t)];']);
        eval(['con_Robust = con_Robust + [ E_bs_',num2str(n_bs),'(t-1) + eta_BS*P_bsc_',num2str(n_bs),'(t)*3600 - 3600*P_bsd_',num2str(n_bs),'(t)/eta_BS  >= 0];']);
        eval(['con_Robust = con_Robust + [ E_bs_',num2str(n_bs),'(t) ==  E_bs_',num2str(n_bs),'(t-1) + eta_BS*P_bsc_',num2str(n_bs),'(t)*3600 - 3600*P_bsd_',num2str(n_bs),'(t)/eta_BS];']);
        end 
    end
            
    %gas balance lim
    con_Robust = con_Robust + [G_L(t) <= eta_P2G*P2Gf*P_p2g(t) - ( a_GT_1*P_gt_1(t)+b_GT_1)- ( a_GT_2*P_gt_2(t)+b_GT_2) + G_buy(t) <= 1.1*G_L(t)]; 
    
    %power lim of P2G
    con_Robust = con_Robust + [P_p2g(t) <= Pmax_P2G(t)];
    con_Robust = con_Robust + [P_p2g(t) >= Pmin_P2G(t)];
    
    %power lim of BS
        for n_bs = 1 : N_bs
        eval(['con_Robust = con_Robust + [P_bsc_',num2str(n_bs),'(t) >= 0];']);
        eval(['con_Robust = con_Robust + [P_bsd_',num2str(n_bs),'(t) >= 0];']);
        end
    
    %power lim of EB
        for n_eb = 1 : N_eb
        eval(['con_Robust = con_Robust + [P_eb_',num2str(n_eb),'(t) >= Pmin_EB(t)];']);
        end
    
    %power lim of GT
        for n_gt = 1 : N_gt
        eval(['con_Robust = con_Robust + [P_gt_',num2str(n_gt),'(t) >= 0];']);
        end
end
%% ===============================Define General Constraints=============================== %%
con_Test = [];
for t = 1:n_T
    %power lim of CHP
     for n_chp = 1 : N_chp
     eval([' con_Test = con_Test + [P_chp_',num2str(n_chp),'(t) <= Pmax_chp_',num2str(n_chp),'(t)] +[P_chp_',num2str(n_chp),'(t) >= Pmin_chp_',num2str(n_chp),'(t)];']);
     eval([' con_Test = con_Test + [Q_chp_',num2str(n_chp),'(t) <= Qmax_chp_',num2str(n_chp),'(t)] +[Q_chp_',num2str(n_chp),'(t) >= Qmin_chp_',num2str(n_chp),'(t)];']);
     end
    %power lim of BS
     for n_bs = 1 : N_bs
     eval([' con_Test = con_Test + [P_bsc_',num2str(n_bs),'(t) <= Pmax_BS_C(t)];']);
     eval([' con_Test = con_Test + [P_bsd_',num2str(n_bs),'(t) <= Pmax_BS_D(t)];']);
     end
    %power lim of EB
    for n_eb = 1 : N_eb
    eval([' con_Test = con_Test + [P_eb_',num2str(n_eb),'(t) <= Pmax_EB(t)];']);
    end
    %power lim of GT
    con_Test = con_Test + [Pmin_GT_1 <= P_gt_1 <= Pmax_GT_1];
    con_Test = con_Test + [Pmin_GT_2<= P_gt_2 <= Pmax_GT_2];
    %general constraints of Nodal Voltage +-6%

    con_Test = con_Test + [0.94*repmat(V_0,n_bus,1) <= V_bus(t,:)' <= 1.06*repmat(V_0,n_bus,1) ];
end
    %general constraints of carbon emissions
    con_Test = con_Test + [repmat(e_chp,1,n_T)*(P_chp_1 + P_chp_2 + P_chp_3)...
        + repmat(e_gt,1,n_T)*(P_gt_1 + P_gt_2) + repmat(e_grid,1,n_T)*P_buy + repmat(e_gas,1,n_T)*G_L - repmat(e_p2g,1,n_T)*(P_p2g) <= (1-lambda_e)*e_max];
%% ===============================Define objective functions=============================== %%
constraints_all = con_Test + con_Robust + con_General;%integrate all MEMG constraints
%nominal sate
obj_MEMG = C_buy_E'*P_buy - C_sell_E(1:n_T)'*P_sell + C_chp_1'*P_chp_1 + C_chp_2'*P_chp_2 + C_chp_3'*P_chp_3...
    - C_BS_ch'*(P_bsc_1 + P_bsc_2) + C_BS_dch'*(P_bsd_1 + P_bsd_2) + C_P2G'*P_p2g + C_EB'*(P_eb_1 + P_eb_2) + C_Gas'*G_buy;
obj_MEMG = obj_MEMG + c_e_car'*(e_chp*(P_chp_1 + P_chp_2 + P_chp_3) + e_gt*(P_gt_1 + P_gt_2) + e_gas*G_L(1:n_T) + e_grid*P_buy - e_p2g*P_p2g);

disp("Modeling complete");
tic%start of timer
op=sdpsettings('solver','gurobi','savesolveroutput',1,'gurobi.ResultFile','MyFile.lp','verbose',0);%set solver
%op=sdpsettings('solver','cplex','savesolveroutput',1);
disp("========================Calculation in process======================")
sol_all = optimize(constraints_all,obj_MEMG,op);
if sol_all.problem == 0
disp("Solved successfully!");
op_sol = sprintf('The minimized operation cost of MEMG is %d (DO-method).',sol_all.solveroutput.result.objval);
disp(op_sol);
elseif sol_all.problem ~= 0
    [model,~] =  export(constraints_all, obj_MEMG,sdpsettings('solver','gurobi'));
    iis = gurobi_iis(model);
    gurobi_write(model , 'MyFile.lp');
    find(iis.Arows)
    msg = 'Model infeasible!';
    error(msg)
end
toc%end of timer

%save output
af_MEMG = [value(af_gt_1) value(af_gt_2) value(af_bsc_1) value(af_bsc_2) value(af_bsd_1) value(af_bsd_2) value(af_p2g)];
xE_MEMG = [value(P_chp_1) value(P_chp_2) value(P_chp_3) value(P_gt_1) value(P_gt_2) -value(P_eb_1) -value(P_eb_2)...
    -value(P_bsc_1) -value(P_bsc_2) value(P_bsd_1) value(P_bsd_2) -value(P_p2g) value(P_buy) -value(P_sell)];
VE_MEMG = value(V_bus)/V_0;%to 1 p.u.
PL_MEMG = [value(P_flow_in) value(P_flow_out)]/1e3;%to MW
af_dev_Vp = value([0 0 0 0 0 0 0 0 af_gt_1(t) 0 0 0 0 (af_bsd_1(t)-af_bsc_1(t)) 0 0 0 0 -af_p2g(t) 0 0 0 0 0 0 (af_bsd_2(t)-af_bsc_2(t)) 0 0 0 0 0 af_gt_2(t)]);%决策向量堆叠
af_dev_Vn = value([0 0 0 0 0 0 0 0 af_gt_1(t) 0 0 0 0 (af_bsd_1(t)-af_bsc_1(t)) 0 0 0 0 -af_p2g(t) 0 0 0 0 0 0 (af_bsd_2(t)-af_bsc_2(t)) 0 0 0 0 0 af_gt_2(t)]);%决策向量堆叠
afQ_dev = value([0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.333 0 0 0 0.333 0 0 0.333 0 0 0 0 0 0 0 0]);%无功决策向量堆叠
af_nbus_VP = repmat(af_dev_Vp,n_T,1);af_nbus_VN =  repmat(af_dev_Vn,n_T,1);afQ_nbus = repmat(afQ_dev,n_T,1);
car = value(repmat(e_chp,1,n_T)*(P_chp_1 + P_chp_2 + P_chp_3)...
        + repmat(e_gt,1,n_T)*(P_gt_1 + P_gt_2) + repmat(e_grid,1,n_T)*P_buy + repmat(e_gas,1,n_T)*G_L - repmat(e_p2g,1,n_T)*(P_p2g));