clear;
clc;
options = optimset('maxiter',5e5,'tolX',1e-10,'tolfun',1e-10,'maxfunevals',5e5);
load parameters_UK_final_f.mat;

sigma=parameters_UK_final_f(1);
B=parameters_UK_final_f(2);
Nb = parameters_UK_final_f(3);
alpha1=parameters_UK_final_f(4);
alpha2=parameters_UK_final_f(5);
alpha3=parameters_UK_final_f(6);
beta=parameters_UK_final_f(7);
cost=parameters_UK_final_f(8);
eta=parameters_UK_final_f(9);                                   
A=parameters_UK_final_f(10); 
theta=parameters_UK_final_f(11);
epsilon=0.001;
opt = 1;
fees = [0.0180,0.0095,0.0033,0.00434]; % Cash, card, IF, debit
left = find(year(Data_TT.Dates) == 2000,1,"first");
right = find(year(Data_TT.Dates) == 2012,1,"last");
i_reserve = mean(Data_TT.BankRate(left:right))/100; 
pi_m = mean(Data_TT.Inflation(left:right));
pi_cbdc = pi_m;
%% Cost of Payments
mu_c = fees(1);
mu_d = fees(2)+fees(3);
mu_db = fees(3);
mu_cb = fees(4);
%% Real Interest Rate
r = 1/beta - 1;
%% Nominal Interest Rate
i_m = (1+r)*(1+pi_m)-1;
%% Real Interest Rate on Reserves
r_reserves = (1+i_reserve)./(1+pi_m)-1;
%% DM Utility Function and Derivatives
u = @(x) ((x+epsilon).^(1-sigma)-epsilon.^(1-sigma))./(1-sigma);
du = @(x) (x+epsilon).^(-sigma);
d2u = @(x) (x+epsilon).^(-sigma-1)*(-sigma);
%% DM Production Cost and Derivatives
c = @(x) x;                         
dc = @(x) 1;
d2c= @(x) 0;
%% Y Star
y_star = fminsearch(@(x)(du(x)-dc(x)).^2,1);                     
%% Trading Protocol with Kalai Bargaining
g_func = @(x) (1-theta).*u(x)+theta.*c(x);                        
dg_func = @(x) (1-theta).*du(x)+theta.*dc(x);
d2g_func = @(x) (1-theta).*d2u(x)+theta.*d2c(x);
%% Trading Quantity and Payment
y_grid = 0:0.001:y_star;
p_grid = g_func(y_grid);
p_grid_d = g_func(y_grid)/(1-mu_d);
p_grid_c = g_func(y_grid)/(1-mu_c);
p_star_d = p_grid_d(end);
p_star_c = p_grid_c(end);
p_star = max(p_star_c,p_star_d);
y_func = @(p) interp1(p_grid,y_grid,min(p,max(p_grid)));          
p_func = @(x) min(x,max(p_grid));
%% Liquidity Premium
lambda_func = @(z) max(du(y_func(z))./dg_func(y_func(z))-1,0); 
dlambda_func = @(z) (d2u(y_func(z)).*dg_func(y_func(z))-du(y_func(z)).*d2g_func(y_func(z)))./...
    (dg_func(y_func(z)).^2).*(z<=p_grid(end));
%% Production Function
f = @(k) A*k.^(eta);
df = @(k) A*eta*k.^(eta-1);
%% Loand Demand, page 17 
df_inv = @(rho) ((1+rho)./(A.*eta)).^(1/(eta-1)); 
%% Xi Equation (11)
xi_func = @(rho) (1-rq_ratio)*max(1+rho,1+r_reserves)+rq_ratio.*(1+r_reserves)-cost;
%% Maximum Return on Loans: R_l <= 1/beta
rho_bar = fzero(@(x) xi_func(x)-1/beta,[1e-10,2]);
rho_bar = min(rho_bar,1/beta-1+cost);
rho_grid = linspace(-.05,rho_bar,1000);
%% Pass Variables   
func_var.beta=beta;
func_var.alpha1=alpha1;
func_var.alpha2=alpha2;
func_var.alpha3=alpha3;
func_var.lambda_func=lambda_func;
func_var.dlambda_func=dlambda_func;
func_var.rq_ratio=rq_ratio;
func_var.cost=cost;
func_var.df_inv = df_inv;
func_var.df=df;
func_var.p_star=p_star;
func_var.p_star_d=p_star_d;
func_var.p_star_c=p_star_c;
func_var.rho_bar = rho_bar;
func_var.Nb=Nb;
func_var.i_m=i_m;
func_var.r_reserves = r_reserves;
func_var.mu_d = mu_d;
func_var.mu_c = mu_c;
func_var.mu_db = mu_db;
func_var.mu_cb = mu_cb;
%% Compute inverse demand for Deposits without CBDC
D_grid = linspace(0,max(p_grid_d),5e3);
ygrid = D_inv_demand(i_m,D_grid,func_var);
%% Compute Loand Demand and Supply
[LS,LD,rho_grid] = LSLD(rho_grid,func_var,Nb,ygrid);
%% Equilibrium when i_cbdc = 0
icbdc_grid=linspace(0,i_m,1000);
z_grid_c = linspace(0,p_grid_c(end),100000);
z_grid_d = linspace(0,p_grid_d(end),100000);
lambda_grid = lambda_func((1-mu_c)*z_grid_c + (1-mu_d)*z_grid_d);
lambda_grid_d = lambda_func((1-mu_d)*z_grid_d);
lambda_grid_c = lambda_func((1-mu_c)*z_grid_c);
y=SS_eq_noCBDC(func_var,Nb,ygrid);
psi = beta*(1-mu_d)*(1 + alpha2*lambda_func((1-mu_d)*y(2)) + alpha3*lambda_func((1-mu_c)*y(1)+(1-mu_d)*y(2)));
drate = 1/psi-1;
loan = fzero(@(x) df(x)-(1+y(end)),[1e-10,100]);
eq_final_nocbdc(1,:) = [y,drate,loan];
%% Equilibrium when i_cbdc >0
eq_cbdc = zeros(size(icbdc_grid,2),5);
drate_cbdc = zeros(size(icbdc_grid,2),1);
eliquidity = zeros(size(icbdc_grid,2),1);
for j=1:length(icbdc_grid)
    icbdc_real  = (1/beta)*(1+pi_m)/(1+icbdc_grid(j))/(1-mu_cb) - 1;
	im_real  = (1/beta)*(1+pi_m)/(1-mu_c) - 1;
    r_cbdc = (1+icbdc_grid(j))/(1+pi_m) - 1;
    % \hat R_d > R_cbdc
    if (1-mu_d)*(1+drate) > (1-mu_cb)*(1+r_cbdc)
        eq_cbdc(j,:) = [y,drate,loan];
        drate_cbdc(j,1) = drate;
        eliquidity(j,1) = y(2);
    % \hat R_d <= R_cbdc
    else
        drate_cbdc(j,1) = (1-mu_cb)*(1+r_cbdc)/(1-mu_d)-1;
        z_eq = SS_eq( alpha1,alpha2,alpha3,z_grid_c,z_grid_d,lambda_grid,lambda_grid_c,lambda_grid_d,im_real,icbdc_real,func_var);    
        lsupply_max = (1-rq_ratio)*(z_eq(1,2)./(1+drate_cbdc(j,1)));
        rho_hat = max(1+r_reserves, (1-mu_db)*(1+drate_cbdc(j,1))+cost) - 1;
        ld_low = fzero(@(x) df(x)-(1+rho_hat),[1e-10,100]); 
        if ld_low<lsupply_max
            % R_d*d
            deposit =  ld_low/(1-rq_ratio)*(1+drate_cbdc(j,1));
            eq_cbdc(j,:) = [z_eq(1),deposit,rho_hat,drate_cbdc(j,1),ld_low];
        else
            rhocbdc = df(lsupply_max)-1;
            eq_cbdc(j,:) = [z_eq,rhocbdc,drate_cbdc(j,1),lsupply_max];  
        end
        eliquidity(j,1) = z_eq(2);
    end
    if j/200 - floor(j/200) == 0
    display(['Equilibrium with CBDC: Iteration over CBDC grid number ' num2str(j)])
    else
    end
end
cash = eq_cbdc(:,1);
deposit = eq_cbdc(:,2);
cbdc = eliquidity - deposit;
loan = eq_cbdc(:,5);
deposit_rate = eq_cbdc(:,4);
loan_rate = eq_cbdc(:,3);
reserves = deposit./(1+deposit_rate) - loan;
Output_CBDC(:,1) = alpha1*y_func((1-mu_c)*cash)+...
    alpha2*y_func(((1-mu_cb)*cbdc + (1-mu_d)*deposit))+...
    alpha3*y_func(min((1-mu_c)*cash(:,1)+(1-mu_d)*deposit+(1-mu_cb)*cbdc,p_star))+...
    B+loan+f(loan)-(1-mu_db)*deposit...
    +(deposit./(1+deposit_rate)-loan)*(1+r_reserves)';
output_index(:,1)=(Output_CBDC(:,1)/Output_CBDC(1,1)-1)*100;
psi_D(:,1) = deposit./(1+deposit_rate);
liquidity = [cbdc./(1+deposit_rate),deposit./(1+deposit_rate),cash./(1+pi_m)];
%% Evaluate Social Costs of Payments without CBDC
GDP = B+alpha1*(1-mu_c)*eq_final_nocbdc(1)+alpha2*(1-mu_d)*eq_final_nocbdc(2)+...
        alpha3*min([(1-mu_c)*eq_final_nocbdc(1)+(1-mu_d)*eq_final_nocbdc(2),p_star])+...
        (1+r_reserves)*rq_ratio*eq_final_nocbdc(2)/(1+drate)-(1-mu_db)*eq_final_nocbdc(2)+...
        (1+eq_final_nocbdc(3)+eta)*(1-rq_ratio)*eq_final_nocbdc(2)/(1+drate)/eta;
socialcosts_d = ((mu_d-mu_db)*eq_final_nocbdc(2)/GDP)*100;
socialcosts_c = (mu_c*eq_final_nocbdc(1)/GDP)*100;
disp(['Without CBDC, the model generates social costs of digital bank deposits worth ' num2str(socialcosts_d) ' % of GDP']);
disp(['Without CBDC, the model generates social costs of cash worth ' num2str(socialcosts_c) ' % of GDP'])
cash2deposit = eq_final_nocbdc(1)/(eq_final_nocbdc(1)+eq_final_nocbdc(2));
%% Welfare
% Entrepreneur Welfare: f(L) - R_lL
welfare_e = f(loan)-(1+loan_rate).*loan;
% Banker Welfare: R_lL + R_r - ((1-mu_db)R_d + cost)D
welfare_banker = (1+loan_rate).*loan +...
    (1+r_reserves)*(deposit./(1+deposit_rate) - loan) -...
    ((1-mu_db)*(1+deposit_rate)+cost).*(deposit./(1+deposit_rate));
% Seller Welfare
liq_sell = alpha1*(p_func((1-mu_c)*cash)) + alpha2*p_func(((1-mu_d)*deposit+(1-mu_cb)*cbdc))+...
            alpha3*p_func(((1-mu_c)*cash+(1-mu_d)*deposit+(1-mu_cb)*cbdc));
welfare_seller_func = @(x) ...
    alpha1*(-c(y_func((1-mu_c)*cash(1)))+p_func(((1-mu_c)*cash(1)))) +...
    alpha2*(-c(y_func((1-mu_d)*deposit(1)))+p_func(((1-mu_d)*deposit(1)))) +...
    alpha3*(-c(y_func((1-mu_c)*cash(1)+(1-mu_d)*deposit(1))))+...
    B/2*log(x*B/2)-B/2 + liq_sell(1); 
welfare_seller = ...
    alpha1*(-c(y_func((1-mu_c)*cash))+p_func(((1-mu_c)*cash))) +...
    alpha2*(-c(y_func((1-mu_d)*deposit+(1-mu_cb)*cbdc))+p_func(((1-mu_d)*deposit+(1-mu_cb)*cbdc))) +...
    alpha3*(-c(y_func((1-mu_c)*cash+(1-mu_d)*deposit+(1-mu_cb)*cbdc)))+...
    B/2*log(B/2)-B/2 + liq_sell; 
welfare_seller_comp = zeros(size(icbdc_grid,2),1);
for i=1:length(icbdc_grid)
    welfare_seller_comp(i) = fzero(@(x) welfare_seller_func(x) - welfare_seller(i),[0.1,2]);
end
% Buyer Welfare
tax = (cash.*(1-(1/(1+pi_m))) + cbdc.*(1-(1+icbdc_grid')./(1+pi_m)) + reserves.*(1-(1+r_reserves)));
welfare_buyer_func = @(x)...
    alpha1*(u(x*y_func((1-mu_c)*cash(1)))-p_func(((1-mu_c)*cash(1)))) +...
    alpha2*(u(x*y_func((1-mu_d)*deposit(1)))-p_func(((1-mu_d)*deposit(1)))) +...
    alpha3*(u(x*y_func((1-mu_c)*cash(1)+(1-mu_d)*deposit(1)))-p_func(((1-mu_c)*cash(1)+(1-mu_d)*deposit(1)))) +...
    B/2*log(x*B/2) - B/2 + tax(1) +...
    cash(1).*((1-mu_c)*(1/(1+pi_m))-1) + cbdc(1).*((1-mu_cb)*(1+icbdc_grid(1))./(1+pi_m)-1) + deposit(1).*(1-mu_d-1./(1+deposit_rate(1))); 
 welfare_buyer = alpha1*(u(y_func((1-mu_c)*cash))-p_func(((1-mu_c)*cash))) +...
    alpha2*(u(y_func((1-mu_d)*deposit+(1-mu_cb)*cbdc))-p_func(((1-mu_d)*deposit+(1-mu_cb)*cbdc))) +...
    alpha3*(u(y_func((1-mu_c)*cash+(1-mu_d)*deposit+(1-mu_cb)*cbdc))-p_func(((1-mu_c)*cash+(1-mu_d)*deposit+(1-mu_cb)*cbdc))) +...
    B/2*log(B/2) - B/2 + tax +...
    cash.*((1-mu_c)*(1/(1+pi_m))-1) + cbdc.*((1-mu_cb)*(1+icbdc_grid')./(1+pi_m)-1) + deposit.*((1-mu_d)-1./(1+deposit_rate)); 
  
welfare_buyer_comp = zeros(size(icbdc_grid,2),1);
for i=1:length(icbdc_grid)
    welfare_buyer_comp(i) = fzero(@(x)welfare_buyer_func(x)-welfare_buyer(i),[0.001,1.5]);
end
welfare_hh = welfare_buyer_comp + welfare_seller_comp - 1;
%% Summarize Results
[loan_max,nmax] = max(loan./loan(1));
i_max=icbdc_grid(nmax);
id = find(loan./loan(1)>=1,1,'last');
id1 = find(loan./loan(1)>1,1,'first');
i_max1=icbdc_grid(id);
i_min1=icbdc_grid(id1);
[output_max,nmax] = max(output_index);
ioutput_max=icbdc_grid(nmax);
id = find(output_index>=0,1,'last');
id1 = find(output_index>0,1,'first');
iy_max1=icbdc_grid(id);
iy_min1=icbdc_grid(id1);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('1.28% (0.33% merchant fee) deposit cost, 1.80% cash cost, 0.--% CBDC cost');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(['CBDC increases Bank Lending if its rate is between ' num2str(i_min1*100) '% and ' num2str(i_max1*100) '%']);
disp(['It increases output if its rate is between ' num2str(iy_min1*100) '% and ' num2str(iy_max1*100) '%']);
disp(['The CBDC rate that maximizes lending is ' num2str(i_max*100) '% and that maximizes output is ' num2str(ioutput_max*100) '%' ]);
disp(['The maximum increase in lending is ' num2str((loan_max-1)*100) '% and the maximum increase in output is ' num2str(output_max) '%' ]);

%% Save Results
save simulation_fees.mat eq_cbdc icbdc_grid output_index i_m pi_m psi_D...
    welfare_buyer_comp welfare_seller_comp welfare_banker welfare_e welfare_hh eliquidity...
    D_grid rho_grid LS rho_bar rq_ratio r_reserves ygrid Nb mu_d mu_cb mu_c mu_db cost...
    LD beta cost liquidity