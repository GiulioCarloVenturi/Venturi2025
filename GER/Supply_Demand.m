clear;
clc;
options = optimset('maxiter',5e5,'tolX',1e-10,'tolfun',1e-10,'maxfunevals',5e5);
load parameters_ITA_final_nf.mat;

sigma=parameters_ITA_final_nf(1);
B=parameters_ITA_final_nf(2);
Nb = parameters_ITA_final_nf(3);
alpha1=parameters_ITA_final_nf(4);
alpha2=parameters_ITA_final_nf(5);
alpha3=parameters_ITA_final_nf(6);
beta=parameters_ITA_final_nf(7);
eta=parameters_ITA_final_nf(9);                                   
A=parameters_ITA_final_nf(10); 
theta=parameters_ITA_final_nf(11);
epsilon=0.001;
rq_ratio = 0.01; %% Updated, reserve requirement after 2012 
opt = 1;
% ECB average deposit rate between 2000 and 2008
ECB_rate = readtimetable("Data_Calibration_ITA.xlsx",Sheet="ECB Rate");
left = find(year(ECB_rate.Dates) == 2000,1,"first");
right = find(year(ECB_rate.Dates) < 2009,1,"last");
time_step = daysact(ECB_rate.Dates(1:end-1),ECB_rate.Dates(2:end))/365;
rate = exp(ECB_rate.DepositRate(1:end-1)/100.*time_step);
time_window = daysact(ECB_rate.Dates(left),ECB_rate.Dates(right))/365;
i_reserve = log(prod(rate(left:right)))/time_window; %% Updated, original 0.0102

%%
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
cost=parameters_ITA_final_nf(8);
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
func_var.p_star_c=p_star_d;
func_var.rho_bar = rho_bar;
func_var.Nb=Nb;
func_var.i_m=i_m;
func_var.r_reserves = r_reserves;
func_var.mu_d = mu_d;
func_var.mu_c = mu_c;
%% Compute inverse demand for Deposits without CBDC
D_grid = linspace(0,max(p_grid_d),5e3);
ygrid = D_inv_demand(i_m,D_grid,func_var);
%% Compute Loand Demand and Supply
[LS] = LSLD(rho_grid,func_var,Nb,ygrid);
[LS_pc,LD,rho_grid] = LSLD_pc(rho_grid,func_var,ygrid);
LS = [LS(1:end-1),LS(end)*2];
LS_pc = [LS_pc(1:end-1),LS(end)*2];
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
%% Loan Supply with CBDC
i_cbdc = 0.01;
[~,ID] = min(abs(icbdc_grid - i_cbdc));
rd_cbdc = eq_cbdc(ID,4);
D_cbdc = eq_cbdc(ID,2);
L_cbdc = eq_cbdc(ID,5);
D_bar = max(D_grid);
xi_func = @(rho) (1-rq_ratio)*max(1+rho,1+r_reserves)+rq_ratio.*(1+r_reserves)-cost;
psi_func = @(x) interp1(ygrid(:,1),ygrid(:,3),min(x,D_bar)).*(x<=D_bar)+beta.*(x>D_bar);
dpsi_func = @(x) interp1(ygrid(:,1),ygrid(:,4),min(x,D_bar)).*(x<=D_bar);
eq_D_func = @(x,y) xi_func(y).*(psi_func(x)./((1-mu_db).*(1-dpsi_func(x).*x/Nb)) )-1;
R_under = fzero(@(x) (1-rq_ratio)*(1+x)+rq_ratio*(1+r_reserves)- ((1-mu_db)*(1+rd_cbdc)+cost),[-1,rho_bar]);
R_upper = fzero(@(x) eq_D_func(D_cbdc,x),[-1,rho_bar]);
[~,ID_upper] = min(abs(rho_grid-R_upper));
[~,ID_under] = min(abs(rho_grid-R_under));
LS_cbdc = [zeros(ID_under,1);L_cbdc*ones(ID_upper-ID_under,1);LS(1,ID_upper+1:end)'];
%% Generate Plot
ID = find(rho_grid == r_reserves);
DS_mon_res = LS(ID)/(1-rq_ratio);
DS_pc_res = LS_pc(ID)/(1-rq_ratio);
DS_mon = max(0,LS/(1-rq_ratio) + DS_mon_res*(LS/(1-rq_ratio) < DS_mon_res));
DS_pc = max(0,LS_pc/(1-rq_ratio) + DS_pc_res*(LS_pc/(1-rq_ratio) < DS_pc_res));
DS_cbdc = max(0,LS_cbdc/(1-rq_ratio));
% Loan Supply
f=figure;
hold on
plot(rho_grid*100,LS_pc ,'r','LineWidth',2)
plot(rho_grid*100,LS,'k','LineWidth',2)
plot(rho_grid*100,LS_cbdc,'k--','LineWidth',2)
plot(rho_grid*100,LD,'b','LineWidth',2)
ylim([0 2])
xlabel('$R_l$','Interpreter','latex')
xtickformat("percentage")
ylabel('$L$','Interpreter','latex')
legend({'Oligopoly','Perfect Competition','CBDC','Loan Demand'},'Location','northwest','Interpreter','latex')
ax=gca;
ax.FontSize = 20;
set(f,'PaperType','A4');
set(f,'PaperOrientation','landscape');
print('ITA_loan_supply','-dpdf','-r200','-fillpage')
% Deposit Supply
f=figure;
hold on
plot(rho_grid(1:end-1)*100,DS_pc(1:end-1),'r','LineWidth',2)
plot(rho_grid(1:end-1)*100,DS_mon(1:end-1),'k','LineWidth',2)
plot(rho_grid(1:end-1)*100,DS_cbdc(1:end-1),'k--','LineWidth',2)
plot([rho_grid(end)*100 rho_grid(end)*100],[DS_mon(end)/2 DS_mon(end)],'k','LineWidth',2')
ylim([0 2])
xlabel('$R_l$','Interpreter','latex')
xtickformat("percentage")
ylabel('$D$','Interpreter','latex')
legend({'Oligopoly', 'Perfect Competition','CBDC'},'Location','northwest','Interpreter','latex')
ax=gca;
ax.FontSize = 20;
set(f,'PaperType','A4');
set(f,'PaperOrientation','landscape');
print('ITA_deposit_supply','-dpdf','-r200','-fillpage')
%% Generate Zoomed Graph
[~,ID_mon] = min(abs(LS - LD));
[~,ID_pc] = min(abs(LS_pc - LD));
[~,ID_cbdc] = min(abs(LS_cbdc' - LD));
xt = [rho_grid(ID_mon)*100 rho_grid(ID_pc)*100 rho_grid(ID_cbdc)*100 ...
    R_under*100 R_upper*100 R_under*100 R_upper*100 r_reserves*100];
yt = [LD(ID_mon)-0.02 LD(ID_pc) LD(ID_cbdc) L_cbdc L_cbdc 1 1 1];
ytd = [DS_mon(ID_mon)-0.02 DS_pc(ID_pc) DS_cbdc(ID_cbdc) L_cbdc/(1-rq_ratio) L_cbdc/(1-rq_ratio) 1 1 1];
str = {'A', 'B', 'C' ,'D', 'E','$\underline R$','$\overline R$','$R_r$'};
% Loan Zoom
f=figure;
hold on
plot(rho_grid*100,LS,'k','LineWidth',2)
plot(rho_grid*100,LS_pc,'r','LineWidth',2)
plot(rho_grid*100,LS_cbdc,'k--','LineWidth',2)
plot(rho_grid*100,LD,'b','LineWidth',2)
plot(rho_grid(ID_mon)*100,LD(ID_mon),'k.','MarkerSize',30)
plot(rho_grid(ID_pc)*100,LD(ID_pc),'k.','MarkerSize',30)
plot(rho_grid(ID_cbdc)*100,LD(ID_cbdc),'k.','MarkerSize',30)
plot(R_under*100,L_cbdc,'k.','MarkerSize',30)
plot(R_upper*100,L_cbdc,'k.','MarkerSize',30)
plot(R_under*100,1,'k.','MarkerSize',30)
plot(R_upper*100,1,'k.','MarkerSize',30)
plot(r_reserves*100,1,'k.','MarkerSize',30)
text(xt,yt+0.01,str,'Interpreter','latex','FontSize',15)
ylim([1 1.15])
xlim([-0.5 +2.5])
xlabel('$R_l$','Interpreter','latex')
xtickformat("percentage")
ylabel('$L$','Interpreter','latex')
legend({'Oligopoly','Perfect Competition','CBDC','Loan Demand'},'Location','northwest','Interpreter','latex')
ax=gca;
ax.FontSize = 20;
set(f,'PaperType','A4');
set(f,'PaperOrientation','landscape');
print('ITA_loan_zoom','-dpdf','-r200','-fillpage')
% Deposit Zoom
f = figure;
hold on
plot(rho_grid(1:end-1)*100,DS_mon(1:end-1),'k','LineWidth',2)
plot(rho_grid(1:end-1)*100,DS_pc(1:end-1),'r','LineWidth',2)
plot(rho_grid(1:end-1)*100,DS_cbdc(1:end-1),'k--','LineWidth',2)
plot([rho_grid(end)*100 rho_grid(end)*100],[DS_mon(end)/2 DS_mon(end)],'k','LineWidth',2')
plot(rho_grid(ID_mon)*100,DS_mon(ID_mon),'k.','MarkerSize',30)
plot(rho_grid(ID_pc)*100,DS_pc(ID_pc),'k.','MarkerSize',30)
plot(rho_grid(ID_cbdc)*100,DS_cbdc(ID_cbdc),'k.','MarkerSize',30)
plot(R_under*100,L_cbdc/(1-rq_ratio),'k.','MarkerSize',30)
plot(R_upper*100,L_cbdc/(1-rq_ratio),'k.','MarkerSize',30)
plot(R_under*100,1,'k.','MarkerSize',30)
plot(R_upper*100,1,'k.','MarkerSize',30)
plot(r_reserves*100,1,'k.','MarkerSize',30)
text(xt,ytd+0.01,str,'Interpreter','latex','FontSize',15)
ylim([1 1.15])
xlim([-0.5 +2.5])
xlabel('$R_l$','Interpreter','latex')
xtickformat("percentage")
ylabel('$D$','Interpreter','latex')
%title('\textbf{Deposit Supply}','Interpreter','latex')
legend({'Oligopoly', 'Perfect Competition','CBDC'},'Location','northwest','Interpreter','latex')
ax=gca;
ax.FontSize = 20;
set(f,'PaperType','A4');
set(f,'PaperOrientation','landscape');
print('ITA_deposit_zoom','-dpdf','-r200','-fillpage')