clc;
clear;
load simulation_fees.mat liquidity eq_cbdc eliquidity
cash_fee = eq_cbdc(:,1);
deposit_fee = eq_cbdc(:,2);
cbdc_fee = eliquidity - deposit_fee;
loan_fee = eq_cbdc(:,5);
deposit_rate_fee = eq_cbdc(:,4);
loan_rate_fee = eq_cbdc(:,3);
clear liquidity eq_cbdc eliquidity
load simulation_nofees.mat liquidity icbdc_grid i_m eq_cbdc eliquidity
cash_nofee = eq_cbdc(:,1);
deposit_nofee = eq_cbdc(:,2);
cbdc_nofee = eliquidity - deposit_nofee;
loan_nofee = eq_cbdc(:,5);
deposit_rate_nofee = eq_cbdc(:,4);
loan_rate_nofee = eq_cbdc(:,3);
clear liquidity eq_cbdc eliquidity
liquidity_fee = [cbdc_fee./(1+deposit_rate_fee),deposit_fee./(1+deposit_rate_fee),cash_fee];
liquidity_nofee = [cbdc_nofee./(1+deposit_rate_nofee),deposit_nofee./(1+deposit_rate_nofee),cash_nofee];
ID_cbdc_nofee = find(liquidity_nofee(:,1) > 0,1,'first');
ID_cbdc_fee = find(liquidity_fee(:,1) > 0,1,'first');
ID_deposit = find(liquidity_nofee(:,2) > liquidity_nofee(1,2),1,'first');
[~,ID_max_nofee] = max(liquidity_nofee(:,2));
[~,ID_max_fee] = max(liquidity_fee(:,2));
ID_min_nofee = find(liquidity_nofee(:,2)<liquidity_nofee(1,2),1,'first');
ID_min_fee = find(liquidity_fee(:,2)<liquidity_fee(1,2),1,'first');
xt = [icbdc_grid(ID_cbdc_nofee) icbdc_grid(ID_cbdc_fee)...
    icbdc_grid(ID_max_nofee) icbdc_grid(ID_max_fee) icbdc_grid(ID_deposit)...
    icbdc_grid(ID_min_nofee) icbdc_grid(ID_min_fee)]*100;
yt = [0 0 liquidity_nofee(ID_max_nofee,2) liquidity_fee(ID_max_fee,2)...
    liquidity_nofee(ID_deposit,2) liquidity_nofee(ID_min_nofee,2) liquidity_fee(ID_min_fee,2)];
str = {'A','B','C','D','E','F','G'};
f=figure;
hold on
plot(icbdc_grid*100,liquidity_fee(:,2),'k',LineWidth=3)
plot(icbdc_grid*100,liquidity_fee(:,1),'b',LineWidth=3)
plot(icbdc_grid*100,liquidity_nofee(:,2),'k--',LineWidth=3)
plot(icbdc_grid*100,liquidity_nofee(:,1),'b--',LineWidth=3)
plot(icbdc_grid(ID_cbdc_nofee)*100,0,'k.','MarkerSize',30)
plot(icbdc_grid(ID_cbdc_fee)*100,0,'k.','MarkerSize',30)
plot(icbdc_grid(ID_max_nofee)*100,liquidity_nofee(ID_max_nofee,2),'k.','MarkerSize',30)
plot(icbdc_grid(ID_max_fee)*100,liquidity_fee(ID_max_fee,2),'k.','MarkerSize',30)
plot(icbdc_grid(ID_deposit)*100,liquidity_nofee(ID_deposit,2),'k.','MarkerSize',30)
plot(icbdc_grid(ID_min_nofee)*100,liquidity_nofee(ID_min_nofee,2),'k.','MarkerSize',30)
plot(icbdc_grid(ID_min_fee)*100,liquidity_fee(ID_min_fee,2),'k.','MarkerSize',30)
text(xt,yt+0.05,str,'Interpreter','latex','FontSize',15)
xlim([0,i_m*100])
xlabel('$i_{CBDC}$', 'interpreter','latex')
xtickformat("percentage");
ax=gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
set(f,'PaperType','A4');
set(f,'PaperOrientation','landscape');
print('GER_eliquidity','-dpdf','-r200','-fillpage')


load parameters_GER_final_f.mat;
sigma=parameters_GER_final_f(1);
B=parameters_GER_final_f(2);
Nb = parameters_GER_final_f(3);
alpha1=parameters_GER_final_f(4);
alpha2=parameters_GER_final_f(5);
alpha3=parameters_GER_final_f(6);
beta=parameters_GER_final_f(7);
cost=parameters_GER_final_f(8);
eta=parameters_GER_final_f(9);                                   
A=parameters_GER_final_f(10); 
theta=parameters_GER_final_f(11);
epsilon=0.001;
rq_ratio = 0.01; 
opt = 1;
fees = [0.0178,0.0092,0.0037,0.0043]; % Cash, card, merchant fee, debit
mu_c = fees(1);
mu_d = fees(2)+fees(3);
mu_db = fees(3);
mu_cb = fees(4);
%% Real Interest Rate
r = 1/beta - 1;
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

ZZZ = mean((1+deposit_rate_fee)./(1+deposit_rate_nofee)-1)*100;
ZZZ1 = mean(lambda_func((1-mu_d)*deposit_fee)./lambda_func(deposit_nofee));
ZZZ2 = mean(lambda_func((1-mu_d)*deposit_fee+(1-mu_cb)*cbdc_fee)./lambda_func(deposit_nofee+cbdc_nofee));