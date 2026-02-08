clear;
clc;
options = optimset('maxiter',5e5,'tolX',1e-10,'tolfun',1e-10,'maxfunevals',5e5);
load parameters_US_final_nf.mat;

sigma=parameters_US_final_nf(1);
B=parameters_US_final_nf(2);
Nb = parameters_US_final_nf(3);
alpha1=parameters_US_final_nf(4);
alpha2=parameters_US_final_nf(5);
alpha3=parameters_US_final_nf(6);
beta=parameters_US_final_nf(7);
eta=parameters_US_final_nf(9);                                   
A=parameters_US_final_nf(10); 
theta=parameters_US_final_nf(11);
epsilon=0.001;
rq_ratio = 0.056; 
opt = 1;
i_reserve = 0.0102;

%% Fees
mu_d = 0.00;
mu_c = 0.00;
mu_db = 0.00;
mu_cb = 0.00;
%% Real Interest Rate
r = 1/beta - 1;
%% Inflation Rate
left = find(year(Data_TT.Dates) == 2000);
right = find(year(Data_TT.Dates) == 2012);
pi_m = (prod(1+Data_TT.Inflation(left:right)/100))^(1/right) - 1; %% Updated, original 0.01515
pi_cbdc = pi_m;
cost=parameters_US_final_nf(8);
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
p_star = p_grid(end);
p_star_d = p_grid_d(end);
y_func = @(p) interp1(p_grid,y_grid,min(p,max(p_grid)));          
y_func_d = @(p) interp1(p_grid_d,y_grid,min(p,max(p_grid_d)));          
p_func = @(x) min(x,max(p_grid));
p_func_d = @(x) min(x,max(p_grid_d));
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
func_var.rho_bar = rho_bar;
func_var.Nb=Nb;
func_var.i_m=i_m;
func_var.r_reserves = r_reserves;
func_var.mu_d = mu_d;
%% Compute inverse demand for Deposits without CBDC
D_grid = linspace(0,max(p_grid_d),5e3);
ygrid = D_inv_demand(i_m,D_grid,func_var);
%% Compute Loand Demand and Supply
[LS,LD,rho_grid] = LSLD(rho_grid,func_var,pi_m,Nb,ygrid);
%% Equilibrium when i_cbdc = 0
icbdc_grid=linspace(0,i_m,1000);
z_grid = linspace(0,p_grid(end),100000);
z_grid_d = linspace(0,p_grid_d(end),100000);
lambda_grid = lambda_func(z_grid);
lambda_grid_d = lambda_func((1-mu_d)*z_grid_d);
y=SS_eq_noCBDC(func_var,Nb,ygrid);
psi = beta*(1-mu_d)*(1 + alpha2*lambda_func((1-mu_d)*y(2)) + alpha3*lambda_func(y(1)+(1-mu_d)*y(2)));
drate = 1/psi-1;
loan = fzero(@(x) df(x)-(1+y(end)),[1e-10,100]);
eq_final_nocbdc(1,:) = [y,drate,loan];
eq_cbdc = zeros(size(icbdc_grid,2),5);
drate_cbdc = zeros(size(icbdc_grid,2),1);
eliquidity = zeros(size(icbdc_grid,2),1);
%% Equilibrium when i_cbdc >0
for j=1:length(icbdc_grid)
    icbdc_real  = (1/beta)*(1+pi_cbdc)/(1+icbdc_grid(j)) - 1; %LHS Eq. 8, p 10 (1-mu)R_d=R_cbdc
    r_cbdc = (1+icbdc_grid(j))/(1+pi_cbdc) - 1;
    % \hat R_d > R_cbdc
    if (1-mu_d)*(1+drate) > (1+r_cbdc)
        eq_cbdc(j,:) = [y,drate,loan];
        drate_cbdc(j,1) = drate;
        eliquidity(j,1) = y(2);
    % \hat R_d <= R_cbdc
    else
        drate_cbdc(j,1) = (1+r_cbdc)/(1-mu_d)-1;
        z_eq = SS_eq( alpha1,alpha2,alpha3,z_grid,z_grid_d,lambda_grid,lambda_grid_d,i_m,icbdc_real);    
        lsupply_max = (1-rq_ratio)*(z_eq(1,2)./(1+drate_cbdc(j,1)));
        rho_hat = max(1+r_reserves, (1-mu_d)*(1+drate_cbdc(j,1))+cost) - 1;
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
Output_CBDC(:,1) = alpha1*(1-mu_c)*cash+...
    alpha2*((1-mu_cb)*cbdc + (1-mu_d)*deposit)+...
    alpha3*min((1-mu_c)*cash(:,1)+(1-mu_d)*deposit+(1-mu_cb)*cbdc,p_star)+...
    B+loan+f(loan)-(1-mu_db)*deposit...
    +(deposit./(1+deposit_rate)-loan)*(1+r_reserves)';
output_index(:,1)=(Output_CBDC(:,1)/Output_CBDC(1,1)-1)*100;
psi_D(:,1) = deposit./(1+deposit_rate);
liquidity = [cbdc./(1+deposit_rate),deposit./(1+deposit_rate),cash./(1+pi_m)];
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
    alpha3*(-c(y_func((1-mu_c)*cash(1)+(1-mu_d)*deposit(1)))+p_func(((1-mu_c)*cash(1)+(1-mu_d)*deposit(1)))) +...
    B/2*log(x*B/2)-B/2 + liq_sell(1);
welfare_seller = ...
    alpha1*(-c(y_func((1-mu_c)*cash))+p_func(((1-mu_c)*cash))) +...
    alpha2*(-c(y_func((1-mu_d)*deposit+(1-mu_cb)*cbdc))+p_func(((1-mu_d)*deposit+(1-mu_cb)*cbdc))) +...
    alpha3*(-c(y_func((1-mu_c)*cash+(1-mu_d)*deposit+(1-mu_cb)*cbdc))+p_func(((1-mu_c)*cash+(1-mu_d)*deposit+(1-mu_cb)*cbdc))) +...
    B/2*log(B/2)-B/2 + liq_sell;
welfare_seller_comp = zeros(size(icbdc_grid,2),1);
for i=1:length(icbdc_grid)
    welfare_seller_comp(i) = fzero(@(x) welfare_seller_func(x) - welfare_seller(i),[0.5,1.5]);
end
% Buyer Welfare
tax = (cash.*(1-(1/(1+pi_m))) + cbdc.*(1-(1+icbdc_grid')./(1+pi_m)) + reserves.*(1-(1+r_reserves)));
welfare_buyer_func = @(x)...
    alpha1*(u(x*y_func((1-mu_c)*cash(1)))-p_func(((1-mu_c)*cash(1)))) +...
    alpha2*(u(x*y_func((1-mu_d)*deposit(1)))-p_func(((1-mu_d)*deposit(1)))) +...
    alpha3*(u(x*y_func((1-mu_c)*cash(1)+(1-mu_d)*deposit(1)))-p_func(((1-mu_c)*cash(1)+(1-mu_d)*deposit(1)))) +...
    B/2*log(x*B/2) - B/2 + tax(1) +...
    cash(1).*((1/(1+pi_m))-1) + cbdc(1).*((1+icbdc_grid(1))./(1+pi_m)-1) + deposit(1).*(1-1./(1+deposit_rate(1))); 
 welfare_buyer = alpha1*(u(y_func((1-mu_c)*cash))-p_func(((1-mu_c)*cash))) +...
    alpha2*(u(y_func((1-mu_d)*deposit+(1-mu_cb)*cbdc))-p_func(((1-mu_d)*deposit+(1-mu_cb)*cbdc))) +...
    alpha3*(u(y_func((1-mu_c)*cash+(1-mu_d)*deposit+(1-mu_cb)*cbdc))-p_func(((1-mu_c)*cash+(1-mu_d)*deposit+(1-mu_cb)*cbdc))) +...
    B/2*log(B/2) - B/2 + tax +...
    cash.*((1/(1+pi_m))-1) + cbdc.*((1+icbdc_grid')./(1+pi_m)-1) + deposit.*(1-1./(1+deposit_rate)); 
welfare_buyer_comp = zeros(size(icbdc_grid,2),1);
for i=1:length(icbdc_grid)
    welfare_buyer_comp(i) = fzero(@(x)welfare_buyer_func(x)-welfare_buyer(i),[0.001,1.5]);
end
welfare_hh = welfare_buyer_comp + welfare_seller_comp - 1;
%% Save Results
save simulation_nofees.mat eq_cbdc icbdc_grid output_index i_m pi_m psi_D...
    welfare_buyer_comp welfare_seller_comp welfare_banker welfare_e welfare_hh eliquidity liquidity
%% Summarize Results
[loan_max,nmax] = max(eq_cbdc(:,end)./eq_cbdc(1,end));
i_max=icbdc_grid(nmax);
id = find(eq_cbdc(:,end)./eq_cbdc(1,end)>=1,1,'last');
id1 = find(eq_cbdc(:,end)./eq_cbdc(1,end)>1,1,'first');
i_max1=icbdc_grid(id);
i_min1=icbdc_grid(id1);
[output_max,nmax] = max(output_index(:,1));
ioutput_max=icbdc_grid(nmax);
id = find(output_index(:,1)>=0,1,'last');
id1 = find(output_index(:,1)>0,1,'first');
iy_max1=icbdc_grid(id);
iy_min1=icbdc_grid(id1);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('Costless Payments');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(['CBDC increases Bank Lending if its rate is between ' num2str(i_min1*100) '% and ' num2str(i_max1*100) '%']);
disp(['It increases output if its rate is between ' num2str(iy_min1*100) '% and ' num2str(iy_max1*100) '%']);
disp(['The CBDC rate that maximizes lending is ' num2str(i_max*100) '% and that maximizes output is ' num2str(ioutput_max*100) '%' ]);
disp(['The maximum increase in lending is ' num2str((loan_max-1)*100) '% and the maximum increase in output is ' num2str(output_max) '%' ]);