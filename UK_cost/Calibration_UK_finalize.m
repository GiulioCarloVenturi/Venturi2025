clear;
clc;
load Calibration_UK_f.mat

%% Parameters and Data between 2000 and 2012
left = find(year(Data_TT.Dates) == 2000);
right = find(year(Data_TT.Dates) == 2012);
pi_m = mean(Data_TT.Inflation(left:right));
i_l = mean(Data_TT.LoanRate(left:right))/100;
i_D = mean(Data_TT.DepositRate(left:right))/100;
alpha1 = 0.75*0.02 + 0.25*0.00; %%%
alpha2 = 0.75*0.12 + 0.25*1.00; %%%
alpha3 = 0.75*0.86 + 0.25*0.00; %%%
epsilon = 0.001; 
beta = 0.99; 
i_reserves = 0.00; %mean(Data_TT.BankRate(left:right))/100; 
r_reserves = (1+i_reserves)/(1+pi_m)-1;
i_m = 1/beta*(1+pi_m)-1;
theta = parameters_UK_f(2);
sigma = parameters_UK_f(3);
B = parameters_UK_f(4);
merchant_fee = 0.0033; %%%
mu_d = 0.0090 + merchant_fee; %%%
mu_c = 0.0180; %%% 
mu_db = merchant_fee;
%% Cost of Handling Deposits
Deposit_Cost = rmmissing(readtimetable("Data_Calibration_UK.xlsx", "UseExcel", false,Sheet="Deposit Cost"));
Total_Assets = convert2annual(rmmissing(readtimetable("Data_Calibration_UK.xlsx", "UseExcel", false,Sheet="Assets")),"Aggregation","mean");
Cost_TT = synchronize(Deposit_Cost,Total_Assets,'intersection');
cost = mean((Cost_TT.Costs)./Cost_TT.Assets);
%% Functions
u = @(x) ((x+epsilon).^(1-sigma)-epsilon.^(1-sigma))./(1-sigma);
du = @(x)(x+epsilon).^(-sigma); 
d2u = @(x)(x+epsilon).^(-sigma-1)*(-sigma);
c = @(x)x;                                                         
dc = @(x)1;
d2c = @(x)0;
y_star = 1- epsilon;                    
g_func = @(x) (1-theta).*u(x)+theta.*c(x);                      
dg_func = @(x) (1-theta).*du(x)+theta.*dc(x);
d2g_func = @(x)(1-theta).*d2u(x)+theta.*d2c(x);
y_grid = linspace(0,y_star,1e4);
p_grid = g_func(y_grid);
p_grid_d = g_func(y_grid)/(1-mu_d);
p_grid_c = g_func(y_grid)/(1-mu_c);
p_star_d = p_grid_d(end);
p_star_c = p_grid_c(end);
y_func = @(p) interp1(p_grid,y_grid,min(p,max(p_grid)));          
p_func = @(x) min(x,max(p_grid));
lambda_func = @(z) max(du(y_func(z))./dg_func(y_func(z))-1,0);          
dlambda_func = @(z) ((d2u(y_func(z)).*dg_func(y_func(z))-du(y_func(z)).*d2g_func(y_func(z)))./dg_func(y_func(z)).^2).*(z<=p_grid(end));
D_grid = linspace(0,max(p_grid_d),5e3);
xi_func = @(rho) (1-rq_ratio)*max([1+rho,1+r_reserves])+rq_ratio*(1+r_reserves) - cost;
rho_bar = fzero(@(x) xi_func(x)-1/beta,[1e-10,2]);
rho_grid = linspace(r_reserves,rho_bar,100);
df = @(k) A*eta*k.^(eta-1); 
df_inv = @(rho) (A*eta./(1+rho)).^(1/(1-eta));
%% Pass Variables
otherpar.beta = beta;
otherpar.epsilon = epsilon;
otherpar.cost = cost; 
otherpar.pi_m = pi_m;
otherpar.spread = (i_l-i_D)/(1+pi_m);
otherpar.i_D = i_D;
otherpar.i_l = i_l;
otherpar.alpha1 = alpha1;
otherpar.alpha2 = alpha2;
otherpar.alpha3 = alpha3;
otherpar.rq_ratio = rq_ratio;
otherpar.eta = eta; 
otherpar.lambda_func = lambda_func;
otherpar.dlambda_func = dlambda_func;
otherpar.p_star_d = p_star_d;
otherpar.p_star_c = p_star_c;
otherpar.rho_bar = rho_bar;
otherpar.r_reserves= r_reserves;
otherpar.i_m = i_m;
otherpar.mu_d = mu_d;
otherpar.mu_c = mu_c;
otherpar.df_inv = df_inv;
otherpar.df = df;
otherpar.mu_db = mu_db;
%% Calibrate A and Nb to match loan-deposit spread
psiz_grid = D_inv_demand(i_m,D_grid,otherpar);
A = fzero(@(x)calibration_func_v3(x,otherpar,psiz_grid),[0.5,2]);
[drate,nb,eq_final,d2c]= calibration_func_v3(A,otherpar,psiz_grid);
%% Save Parameters
parameters_UK_final_f=[sigma,B,nb,alpha1,alpha2,alpha3,beta,cost,eta,A,theta];
save parameters_UK_final_f.mat parameters_UK_final_f Data_TT rq_ratio