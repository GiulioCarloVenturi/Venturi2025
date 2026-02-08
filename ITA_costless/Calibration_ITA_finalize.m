clear;
clc;
load Calibration_ITA_nf.mat

%% Parameters and Data between 2000 and 2012
left = find(year(Data_TT.Dates) == 2000);
right = find(year(Data_TT.Dates) == 2012);
pi_m = prod(1+(Data_TT.Inflation(left:right))/100)^(1/(right-left+1)) - 1;
i_l = mean(Data_TT.LoanRate(left:right))/100;
i_D = mean(Data_TT.DepositRate(left:right))/100;
alpha1 = 0.7083*0.31 + 0.2917*0.00; 
alpha2 = 0.7083*0.03 + 0.2917*1.00; 
alpha3 = 0.7083*0.66 + 0.2917*0.00; 
epsilon = 0.001; 
beta = 0.99;
rq_ratio = 0.02; 
i_reserves = mean(Data_TT.MRRRate(left:right))/100; 
r_reserves = (1+i_reserves)/(1+pi_m)-1;
i_m = 1/beta*(1+pi_m)-1;
pi_eff = (1+pi_m)/(1+i_reserves)-1;
trade_prob = parameters_ITA_nf(end);
alpha1 = alpha1*trade_prob;
alpha2 = alpha2*trade_prob;
alpha3 = alpha3*trade_prob;
theta = parameters_ITA_nf(2);
sigma = parameters_ITA_nf(3);
B = parameters_ITA_nf(4);
mu_d = 0.00;
%% Cost of Handling Deposits
Deposit_Cost = rmmissing(readtimetable("Data_Calibration_ITA.xlsx", "UseExcel", false,Sheet="Deposit Cost"));
Total_Assets = convert2annual(rmmissing(readtimetable("Data_Calibration_ITA.xlsx", "UseExcel", false,Sheet="Total Assets")),"Aggregation","mean");
Cost_TT = synchronize(Deposit_Cost,Total_Assets,'intersection');
cost = mean((Cost_TT.Operativi+Cost_TT.Personale)./Cost_TT.TotalAssets);
%% Calibrate eta
Commercial_Loan = flipud(rmmissing(readtimetable("Data_Calibration_ITA.xlsx",Sheet="Commercial Loan Rate")));
left = find(year(Commercial_Loan.Dates) == 1999 & month(Commercial_Loan.Dates) == 1);
right = find(year(Commercial_Loan.Dates) == 2020 & month(Commercial_Loan.Dates) == 2);
Real_Loans = Commercial_Loan.LogRealLoans(left:right);
Real_Rate = Commercial_Loan.LogRealRate(left:right);
lm = fitlm(Real_Rate,Real_Loans);
elasticity = lm.Coefficients.Estimate(2);
eta = 1 + 1/elasticity;
%% Functions
u = @(x) ((x+epsilon).^(1-sigma)-epsilon.^(1-sigma))./(1-sigma);
du = @(x)(x+epsilon).^(-sigma); 
d2u = @(x)(x+epsilon).^(-sigma-1)*(-sigma);
c = @(x)x;                                                         
dc = @(x)1;
d2c = @(x)0;
y_star = fminsearch(@(x)(du(x)-dc(x)).^2,1);                    
g_func = @(x) (1-theta).*u(x)+theta.*c(x);                      
dg_func = @(x) (1-theta).*du(x)+theta.*dc(x);
d2g_func = @(x)(1-theta).*d2u(x)+theta.*d2c(x);
y_grid = linspace(0,y_star,1e4);
p_grid = g_func(y_grid);
p_grid_d = g_func(y_grid)/(1-mu_d);
p_star_d = p_grid_d(end);
y_func = @(p) interp1(p_grid,y_grid,min(p,max(p_grid)));          
p_func = @(x) min(x,max(p_grid));
y_func_d = @(p) interp1(p_grid_d,y_grid,min(p,max(p_grid_d)));          
p_func_d = @(x) min(x,max(p_grid_d));
lambda_func = @(z) max(du(y_func(z))./dg_func(y_func(z))-1,0);          
dlambda_func = @(z) ((d2u(y_func(z)).*dg_func(y_func(z))-du(y_func(z)).*d2g_func(y_func(z)))./dg_func(y_func(z)).^2).*(z<=p_grid(end));
D_grid = linspace(0,max(p_grid_d),5e3);
xi_func = @(rho) (1-rq_ratio)*max([1+rho,1+r_reserves])+rq_ratio*(1+r_reserves) - cost;
rho_bar = fzero(@(x) xi_func(x)-1/beta,[1e-10,2]);
rho_grid = linspace(r_reserves,rho_bar,100);
df = @(k)A*eta*k.^(eta-1); 
df_inv = @(rho)(A*eta./(1+rho)).^(1/(1-eta));
%% Pass Variables
otherpar.beta = beta;
otherpar.epsilon = epsilon;
otherpar.cost = cost; 
otherpar.pi_eff = pi_eff;
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
otherpar.rho_bar = rho_bar;
otherpar.r_reserves= r_reserves;
otherpar.i_m = i_m;
otherpar.mu_d = mu_d;
otherpar.df_inv = df_inv;
otherpar.df = df;
%% Calibrate A and Nb to match loan-deposit spread
psiz_grid = D_inv_demand(i_m,D_grid,otherpar);
A = fzero(@(x)calibration_func_v3(x,otherpar,psiz_grid),[1,2]);
[drate,nb,eq_final,d2c]= calibration_func_v3(A,otherpar,psiz_grid);
%% Save Parameters
parameters_ITA_final_nf=[sigma,B,nb,alpha1,alpha2,alpha3,beta,cost,eta,A,theta];
save parameters_ITA_final_nf.mat parameters_ITA_final_nf Data_TT