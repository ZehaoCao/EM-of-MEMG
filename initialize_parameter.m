%Initial parameters: scripted by Zehao Cao
%% Define the parameters of MEMG
cost_factor = 1;%cost coefficient
power_factor = 1;%power coefficient
%% Paremeters of CHP
N_chp = 3;
Pmax_chp_1 = ones(n_T, 1) * 15*power_factor;%power limit of CHP: MW
Pmin_chp_1 = ones(n_T, 1) * 3.5*power_factor;
Qmax_chp_1 = ones(n_T, 1) * 12*power_factor;%power limit of CHP: MVar
Qmin_chp_1 = ones(n_T, 1) * 0;
H_E = 1.2;%heat-elec ratio
C_chp_1 = ones(n_T, 1)*57.1*cost_factor*dt;%unit cost of CHP :$/MWh
ramp_chp_1 = 6*power_factor;%ramping limit of CHP: MW
down_chp_1 = 6*power_factor;

Pmax_chp_2 = ones(n_T, 1) * 10.5*power_factor;%power limit of CHP: MW
Pmin_chp_2 = ones(n_T, 1) * 4.5*power_factor;
Qmax_chp_2 = ones(n_T, 1) * 6*power_factor;%power limit of CHP: MVar
Qmin_chp_2 = ones(n_T, 1) * 0;
C_chp_2 = ones(n_T, 1)*60.8*cost_factor*dt;%unit cost of CHP :$/MWh
ramp_chp_2 = 4;%ramping limit of CHP: MW
down_chp_2 = 4;

Pmax_chp_3 = ones(n_T, 1) * 7.5*power_factor;%power limit of CHP: MW
Pmin_chp_3 = ones(n_T, 1) * 2*power_factor;
Qmax_chp_3 = ones(n_T, 1) * 7.5*power_factor;%power limit of CHP: MVar
Qmin_chp_3 = ones(n_T, 1) * 0;
C_chp_3 = ones(n_T, 1)*62.4*cost_factor*dt;%unit cost of CHP :$/MWh
ramp_chp_3 = 4;%ramping limit of CHP: MW
down_chp_3 = 4;
%% Paremeters of BS
N_bs = 2;
Pmax_BS_C = ones(n_T, 1) * 3*power_factor;%power limit of BSC: MW
Pmin_BS_C = zeros(n_T, 1);
Pmax_BS_D = ones(n_T, 1) * 3*power_factor;%power limit of BSD: MW
Pmin_BS_D = zeros(n_T, 1);
E_BS_max = ones(n_T, 1) * 6*power_factor*3600;%capacity limit of BS: MWh
E_BS_min = zeros(n_T, 1);
eta_BS = 0.8;%efficiency factor of BS
C_BS_ch = ones(n_T, 1)*56.2*cost_factor*dt;%unit cost of BS :$/MWh
C_BS_dch = ones(n_T, 1)*67.3*cost_factor*dt;
%% Paremeters of EB
N_eb = 2;
Pmax_EB = ones(n_T, 1) *6*power_factor;%power limit of EB: MW
Pmin_EB = ones(n_T, 1)*1*power_factor;
eta_EB = 0.9;%efficiency factor of EB
C_EB = ones(n_T, 1)*1.9*cost_factor*dt;%unit cost of EB :$/MWh

%% Paremeters of GT
N_gt = 2;
Pmax_GT_1 = ones(n_T, 1) * 6*power_factor;%power limit of GT: MW
Pmin_GT_1 = ones(n_T, 1) * 0;
a_GT_1 = 35; b_GT_1 = 8;%gas consumption coefficient of GT: kg/MW

Pmax_GT_2 = ones(n_T, 1) * 4.5*power_factor;
Pmin_GT_2 = ones(n_T, 1) * 0;
a_GT_2 = 30; b_GT_2 = 5;
%% Paremeters of P2G
N_p2g = 1;
Pmax_P2G = ones(n_T, 1) * 6*power_factor;%power limit of P2G: MW
Pmin_P2G = ones(n_T, 1) * 0;
eta_P2G = 0.80;%efficiency factor of P2G
P2Gf=(0.50*(3600/46));%elec-gas ratio
C_P2G = ones(n_T, 1)*48.9*cost_factor*dt;%unit cost of P2G :$/MWh
%% Paremeters of WT
N_wt = 3;
[PW_max,PW_min,PV_max,PV_min] = DERs_org_data;
PW_factor = [0.5 0.5 0.5];%set power ratio of WT: p.u.
PW_1 = PW_factor(1)*PW_max; PW_2 = PW_factor(2)*PW_max; PW_3 = PW_factor(3)*PW_max;%define WTs
PV = (1/32)*PV_max;%define PVs
%% Network configuration
[Bra_B,R,X,n_bus]=LDPF_parameter;%get power network data
n_dev = N_chp + N_bs*2 + N_eb + N_gt + N_p2g + N_wt;
Conf =zeros(n_dev,n_bus);%define network configration P
Conf(1,18) = 1;Conf(2,22) = 1;Conf(3,25) = 1;%set CHP at 18、22、25
Conf(4,15) = 1;Conf(5,15) = 1;Conf(6,27) = 1;Conf(7,27) = 1;%set BS at 15、27
Conf(8,12) = 1;Conf(9,29) = 1;%set EB at 12、29
Conf(10,10) = 1;Conf(11,33) = 1;%set GT at 10、33
Conf(12,20) = 1;%set P2G at 20
Conf(13,3) = 1;Conf(14,16) = 1;Conf(15,32) = 1;%set WT at 3、16、32
Conf(16,2:n_bus)=1;%set pv at all nodes 
Conf_Q = [Conf(1:3,:);Conf(13:15,:)];%define network configration Q
Conf_WT = sum(Conf(13:15,2:n_bus));
Conf_WT(Conf_WT~=0) = 1/3;%define uncertainty ratio of 3 WTs
Conf_DEV = sum(Conf(1:12,2:n_bus));
Conf_DEV(Conf_DEV~=0) = 1;
V_0 = 12.66/1e3;% base voltage
%% Parameters of purchase&sell
C_buy_E = ones(n_T, 1)*118.1*cost_factor;%price of elec purchase
C_sell_E = [ones(7, 1)*115.2*0.9;ones(13, 1)*95.2;ones(4, 1)*115.2*0.8]*cost_factor;%price of elec sell
C_Gas = ones(n_T, 1)*32*0.0434*cost_factor*dt;%price of gas purchase

%% Parameters of carbon emissions
e_chp = 0.536;%carbon emissions of CHP: ton/MWh
e_gt = 0.461;%carbon emissions of GT: ton/MWh
e_grid = 0.571;%carbon emissions of GRID: ton/MWh
e_gas = (18.5/3600);%carbon emissions of GAS: 单位ton/MWh
e_p2g = 0.415;%carbon emissions of P2G: ton/MWh
c_e_car = ones(n_T, 1)*20*cost_factor*dt;%%unit cost of carbon emissions: %/ton
lambda_e = 0.10;%carbon reduction ratio
e_max = 440*power_factor;%max carbon emissions

%% Parameters of Loads
%elec load
P_L = [222.35
228.28
227.74
229.38
231.25
235.58
248.71
256.32
274.85
282.96
287.48
291.48
268.94
265.88
269.49
272.23
276.32
281.39
294.78
283.88
268.65
246.51
243.44
249
];
P_L = P_L/(max(P_L));%Normalized to take the shape of its load curve
%heat load
H_L = 5*ones(1, 24)+0.3*[58.3
55.70
54.40
64.87
70.64
70.70
65.29
56.88
47.30
45.95
38.54
32.37
42.43
48.49
54.20
61.14
63.07
65.25
73.66
75.53
70.41
63.59
62.29
53.76];
H_L = H_L*1.5*power_factor;
%gas load
G_L = 20*[5.02
5.04
5.6
8
11.8
12
13.8
13.2
11
9.6
11.1
13.2
12.6
11.2
10.6
11.2
14.2
15.2
13.4
13.08
10.8
10
8
5.2];
G_L = G_L*0.6*power_factor;


