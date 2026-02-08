clear;
clc;
options = optimset('maxiter',5e5,'tolX',1e-10,'tolfun',1e-10,'maxfunevals',5e5);
load parameters_GER_final_f.mat;

sigma=parameters_GER_final_f(1);
B=parameters_GER_final_f(2);
Nb = parameters_GER_final_f(3);
beta=parameters_GER_final_f(7);
cost=parameters_GER_final_f(8);
eta=parameters_GER_final_f(9);                                   
A=parameters_GER_final_f(10); 
theta=parameters_GER_final_f(11);
epsilon=0.001;
rq_ratio = 0.01; 
opt = 1;
fees = [0.0178,0.0092,0.0037,0.0043]; % Cash, card, merchant fee
mu_c = fees(1);
mu_d = fees(2)+fees(3);
mu_db = fees(3);
mu_cb = fees(4);
ECB_rate = readtimetable("Data_Calibration_GER.xlsx",Sheet="ECB Rate");
left = find(year(ECB_rate.Dates) == 2000,1,"first");
right = find(year(ECB_rate.Dates) < 2009,1,"last");
time_step = daysact(ECB_rate.Dates(1:end-1),ECB_rate.Dates(2:end))/365;
rate = exp(ECB_rate.DepositRate(1:end-1)/100.*time_step);
time_window = daysact(ECB_rate.Dates(left),ECB_rate.Dates(right))/365;
i_reserve = log(prod(rate(left:right)))/time_window; 
r = 1/beta - 1;
left = find(year(Data_TT.Dates) == 2003);
right = find(year(Data_TT.Dates) == 2012);
pi_m = (prod(1+Data_TT.Inflation(left:right)/100))^(1/right) - 1; 
i_m = (1+r)*(1+pi_m)-1;
r_reserves = (1+i_reserve)./(1+pi_m)-1;
u = @(x) ((x+epsilon).^(1-sigma)-epsilon.^(1-sigma))./(1-sigma);
du = @(x) (x+epsilon).^(-sigma);
d2u = @(x) (x+epsilon).^(-sigma-1)*(-sigma);
c = @(x) x;                         
dc = @(x) 1;
d2c= @(x) 0;
y_star = fminsearch(@(x)(du(x)-dc(x)).^2,1);                     
g_func = @(x) (1-theta).*u(x)+theta.*c(x);                        
dg_func = @(x) (1-theta).*du(x)+theta.*dc(x);
d2g_func = @(x) (1-theta).*d2u(x)+theta.*d2c(x);
y_grid = 0:0.001:y_star;
p_grid = g_func(y_grid);
p_grid_d = g_func(y_grid)/(1-mu_d);
p_grid_c = g_func(y_grid)/(1-mu_c);
p_star_d = p_grid_d(end);
p_star_c = p_grid_c(end);
p_star = max(p_star_c,p_star_d);
y_func = @(p) interp1(p_grid,y_grid,min(p,max(p_grid)));          
p_func = @(x) min(x,max(p_grid));
lambda_func = @(z) max(du(y_func(z))./dg_func(y_func(z))-1,0); 
dlambda_func = @(z) (d2u(y_func(z)).*dg_func(y_func(z))-du(y_func(z)).*d2g_func(y_func(z)))./...
    (dg_func(y_func(z)).^2).*(z<=p_grid(end));
f = @(k) A*k.^(eta);
df = @(k) A*eta*k.^(eta-1);
df_inv = @(rho) ((1+rho)./(A.*eta)).^(1/(eta-1)); 
xi_func = @(rho) (1-rq_ratio)*max(1+rho,1+r_reserves)+rq_ratio.*(1+r_reserves)-cost;
rho_bar = fzero(@(x) xi_func(x)-1/beta,[1e-10,2]);
rho_bar = min(rho_bar,1/beta-1+cost);
rho_grid = linspace(-.05,rho_bar,1000);
%% Pass Variables   
func_var.beta=beta;
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

alpha1_grid = linspace(parameters_GER_final_f(4)/2,parameters_GER_final_f(4),5);
for n = 3:size(alpha1_grid,2)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(['Cash Reduction ' num2str((alpha1_grid(end)-alpha1_grid(n))*100) 'pp']);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
alpha1=alpha1_grid(n);
alpha2=parameters_GER_final_f(5)+(alpha1_grid(end)-alpha1_grid(n));
alpha3=parameters_GER_final_f(6);
func_var.alpha1=alpha1;
func_var.alpha2=alpha2;
func_var.alpha3=alpha3;
%% Compute inverse demand for Deposits without CBDC
D_grid = linspace(0,max(p_grid_d),5e3);
ygrid = D_inv_demand_fee(i_m,D_grid,func_var);
[LS,LD,rho_grid] = LSLD_fee(rho_grid,func_var,Nb,ygrid);
icbdc_grid=linspace(0,i_m,1000);
z_grid_c = linspace(0,p_grid_c(end),100000);
z_grid_d = linspace(0,p_grid_d(end),100000);
lambda_grid = lambda_func((1-mu_c)*z_grid_c + (1-mu_d)*z_grid_d);
lambda_grid_d = lambda_func((1-mu_d)*z_grid_d);
lambda_grid_c = lambda_func((1-mu_c)*z_grid_c);
y=SS_eq_noCBDC_fee(func_var,Nb,ygrid);
psi = beta*(1-mu_d)*(1 + alpha2*lambda_func((1-mu_d)*y(2)) + alpha3*lambda_func((1-mu_c)*y(1)+(1-mu_d)*y(2)));
drate = 1/psi-1;
loan_no = fzero(@(x) df(x)-(1+y(end)),[1e-10,100]);
eq_final_nocbdc(1,:) = [y,drate,loan_no];
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
        eq_cbdc(j,:) = [y,drate,loan_no];
        drate_cbdc(j,1) = drate;
        eliquidity(j,1) = y(2);
    % \hat R_d <= R_cbdc
    else
        drate_cbdc(j,1) = (1-mu_cb)*(1+r_cbdc)/(1-mu_d)-1;
        z_eq = SS_eq_fee( alpha1,alpha2,alpha3,z_grid_c,z_grid_d,lambda_grid,lambda_grid_c,lambda_grid_d,im_real,icbdc_real,func_var);    
        lsupply_max = (1-rq_ratio)*(z_eq(1,2)./(1+drate_cbdc(j,1)));
        rho_hat = max(1+r_reserves, (1-mu_db)*(1+drate_cbdc(j,1))+cost) - 1;
        ld_low = fzero(@(x) df(x)-(1+rho_hat),[1e-10,100]); 
        if ld_low<lsupply_max
            % R_d*d
            deposit_ =  ld_low/(1-rq_ratio)*(1+drate_cbdc(j,1));
            eq_cbdc(j,:) = [z_eq(1),deposit_,rho_hat,drate_cbdc(j,1),ld_low];
        else
            rhocbdc = df(lsupply_max)-1;
            eq_cbdc(j,:) = [z_eq,rhocbdc,drate_cbdc(j,1),lsupply_max];  
        end
        eliquidity(j,1) = z_eq(2);
    end
end
cash(:,n) = eq_cbdc(:,1);
deposit(:,n) = eq_cbdc(:,2);
cbdc(:,n) = eliquidity - deposit(:,n);
loan(:,n) = eq_cbdc(:,5);
deposit_rate(:,n) = eq_cbdc(:,4);
loan_rate(:,n) = eq_cbdc(:,3);
reserves(:,n) = deposit(:,n)./(1+deposit_rate(:,n)) - loan(:,n);
Output_CBDC(:,n) = alpha1*y_func(min((1-mu_c)*cash(:,n),p_star_c))+...
    alpha2*y_func(min(((1-mu_cb)*cbdc(:,n) + (1-mu_d)*deposit(:,n)),p_star_d))+...
    alpha3*y_func((1-mu_c)*cash(:,n)+(1-mu_d)*deposit(:,n)+(1-mu_cb)*cbdc(:,n))+...
    B+loan(:,n)+f(loan(:,n))-(1-mu_db)*deposit(:,n)...
    +(deposit(:,n)./(1+deposit_rate(:,n))-loan(:,n))*(1+r_reserves)';
output_index(:,n)=(Output_CBDC(:,n)/Output_CBDC(1,n)-1)*100;
%% Summarize Results
[loan_max(n),nmax] = max(loan(:,n)./loan(1,n));
i_max(n)=icbdc_grid(nmax);
id = find(loan(:,n)./loan(1,n)>=1,1,'last');
id1 = find(loan(:,n)./loan(1,n)>1,1,'first');
i_max1(n)=icbdc_grid(id);
i_min1(n)=icbdc_grid(id1);
[output_max(n),nmax] = max(output_index(:,n));
ioutput_max(n)=icbdc_grid(nmax);
id = find(output_index(:,n)>=0,1,'last');
id1 = find(output_index(:,n)>0,1,'first');
iy_max1(n)=icbdc_grid(id);
iy_min1(n)=icbdc_grid(id1);
disp(['CBDC increases Bank Lending if its rate is between ' num2str(i_min1(n)*100) '% and ' num2str(i_max1(n)*100) '%']);
disp(['It increases output if its rate is between ' num2str(iy_min1(n)*100) '% and ' num2str(iy_max1(n)*100) '%']);
disp(['The CBDC rate that maximizes lending is ' num2str(i_max(n)*100) '% and that maximizes output is ' num2str(ioutput_max(n)*100) '%' ]);
disp(['The maximum increase in lending is ' num2str((loan_max(n)-1)*100) '% and the maximum increase in output is ' num2str(output_max(n)) '%' ]);
end

cash_delta = (alpha1_grid(end)-alpha1_grid)*100;
f=figure;
hold on
plot(cash_delta,i_min1*100,'--k','LineWidth',2)
plot(cash_delta,i_max*100,'k','LineWidth',2)
plot(cash_delta,i_max1*100,'--k','LineWidth',2)
ytickformat("percentage")
xlabel('PP Reduction of $\gamma_1$ ','Interpreter','latex','FontSize',15)
ylabel('$i_{cbdc}^*$','Interpreter','latex','FontSize',30)
set(f,'PaperType','A4');
set(f,"PaperOrientation","landscape")
print('GER_cashcomp_rate','-dpdf','-r200','-fillpage')
f=figure;
hold on
plot(cash_delta,(loan_max-1)*100,'k','LineWidth',2)
ytickformat("percentage")
xlabel('PP Reduction of $\gamma_1$ ','Interpreter','latex','FontSize',15)
ylabel('$\Delta$ Loan','Interpreter','latex','FontSize',15)
set(f,'PaperType','A4');
set(f,"PaperOrientation","landscape")
print('GER_cashcomp_loan','-dpdf','-r200','-fillpage')
f=figure;
hold on
plot(cash_delta,output_max,'k','LineWidth',2)
ytickformat("percentage")
xlabel('PP Reduction of $\gamma_1$ ','Interpreter','latex','FontSize',15)
ylabel('$\Delta$ Output','Interpreter','latex','FontSize',15)
set(f,'PaperType','A4');
set(f,"PaperOrientation","landscape")
print('GER_cashcomp_output','-dpdf','-r200','-fillpage')
f=figure;
hold on
plot(icbdc_grid*100,deposit,'LineWidth',2)
xtickformat("percentage")
xlabel('$i_{CBDC}$ ','Interpreter','latex','FontSize',15)
ylabel('Deposits','Interpreter','latex','FontSize',15)
legend({ num2str(round(alpha1_grid(1)/alpha1_grid(end),2)) ...
    num2str(round(alpha1_grid(2)/alpha1_grid(end),2))...
    num2str(round(alpha1_grid(3)/alpha1_grid(end),2))...
    num2str(round(alpha1_grid(4)/alpha1_grid(end),2))...
    num2str(round(alpha1_grid(5)/alpha1_grid(end),2))})
set(f,'PaperType','A4');
set(f,"PaperOrientation","landscape")
print('GER_cashcomp_deposits','-dpdf','-r200','-fillpage')
f=figure;
hold on
plot(icbdc_grid*100,((1+pi_m)*(1+deposit_rate)-1)*100,'LineWidth',2)
xtickformat("percentage")
ytickformat("percentage")
xlabel('$i_{CBDC}$ ','Interpreter','latex','FontSize',15)
ylabel('Deposit Rate','Interpreter','latex','FontSize',15)
legend({ num2str(round(alpha1_grid(1)/alpha1_grid(end),2)) ...
    num2str(round(alpha1_grid(2)/alpha1_grid(end),2))...
    num2str(round(alpha1_grid(3)/alpha1_grid(end),2))...
    num2str(round(alpha1_grid(4)/alpha1_grid(end),2))...
    num2str(round(alpha1_grid(5)/alpha1_grid(end),2))})
set(f,'PaperType','A4');
set(f,"PaperOrientation","landscape")
print('GER_cashcomp_deprate','-dpdf','-r200','-fillpage')
f=figure;
hold on
plot(icbdc_grid*100,loan,'LineWidth',2)
xtickformat("percentage")
xlabel('$i_{CBDC}$ ','Interpreter','latex','FontSize',15)
ylabel('Loans','Interpreter','latex','FontSize',15)
legend({ num2str(round(alpha1_grid(1)/alpha1_grid(end),2)) ...
    num2str(round(alpha1_grid(2)/alpha1_grid(end),2))...
    num2str(round(alpha1_grid(3)/alpha1_grid(end),2))...
    num2str(round(alpha1_grid(4)/alpha1_grid(end),2))...
    num2str(round(alpha1_grid(5)/alpha1_grid(end),2))})
set(f,'PaperType','A4');
set(f,"PaperOrientation","landscape")
print('GER_cashcomp_loans','-dpdf','-r200','-fillpage')
f=figure;
hold on
plot(icbdc_grid*100,((1+pi_m)*(1+loan_rate)-1)*100,'LineWidth',2)
xtickformat("percentage")
ytickformat("percentage")
xlabel('$i_{CBDC}$ ','Interpreter','latex','FontSize',15)
ylabel('Loan Rate','Interpreter','latex','FontSize',15)
legend({ num2str(round(alpha1_grid(1)/alpha1_grid(end),2)) ...
    num2str(round(alpha1_grid(2)/alpha1_grid(end),2))...
    num2str(round(alpha1_grid(3)/alpha1_grid(end),2))...
    num2str(round(alpha1_grid(4)/alpha1_grid(end),2))...
    num2str(round(alpha1_grid(5)/alpha1_grid(end),2))})
set(f,'PaperType','A4');
set(f,"PaperOrientation","landscape")
print('GER_cashcomp_loanrate','-dpdf','-r200','-fillpage')
f=figure;
hold on
plot(icbdc_grid*100,(1+pi_m)*(loan_rate-deposit_rate)*100,'LineWidth',2)
xtickformat("percentage")
ytickformat("percentage")
xlabel('$i_{CBDC}$ ','Interpreter','latex','FontSize',15)
ylabel('Spread','Interpreter','latex','FontSize',15)
legend({ num2str(round(alpha1_grid(1)/alpha1_grid(end),2)) ...
    num2str(round(alpha1_grid(2)/alpha1_grid(end),2))...
    num2str(round(alpha1_grid(3)/alpha1_grid(end),2))...
    num2str(round(alpha1_grid(4)/alpha1_grid(end),2))...
    num2str(round(alpha1_grid(5)/alpha1_grid(end),2))})
set(f,'PaperType','A4');
set(f,"PaperOrientation","landscape")
print('GER_cashcomp_spread','-dpdf','-r200','-fillpage')
f=figure;
hold on
plot(icbdc_grid*100,Output_CBDC,'LineWidth',2)
xtickformat("percentage")
xlabel('$i_{CBDC}$ ','Interpreter','latex','FontSize',15)
ylabel('GDP','Interpreter','latex','FontSize',15)
legend({ num2str(round(alpha1_grid(1)/alpha1_grid(end),2)) ...
    num2str(round(alpha1_grid(2)/alpha1_grid(end),2))...
    num2str(round(alpha1_grid(3)/alpha1_grid(end),2))...
    num2str(round(alpha1_grid(4)/alpha1_grid(end),2))...
    num2str(round(alpha1_grid(5)/alpha1_grid(end),2))})
set(f,'PaperType','A4');
set(f,"PaperOrientation","landscape")
print('GER_cashcomp_gdp','-dpdf','-r200','-fillpage')