function eq=SS_eq_noCBDC(func_var,Nb,y)
%% Parameters
rq_ratio = func_var.rq_ratio;
cost = func_var.cost;
rho_bar= func_var.rho_bar;
beta = func_var.beta;
df_inv = func_var.df_inv;
df = func_var.df;
r_reserves = func_var.r_reserves;
mu_d = func_var.mu_d;
%% Functions
xi_func = @(rho) (1-rq_ratio).*max(1+rho,1+r_reserves) + rq_ratio*(1+r_reserves) - cost;
D_bar = max(y(:,1));
psi_func = @(x) interp1(y(:,1),y(:,3),min(x,D_bar)).*(x<=D_bar)+beta.*(x>D_bar);
dpsi_func = @(x) interp1(y(:,1),y(:,4),min(x,D_bar)).*(x<=D_bar);
z_func = @(x) interp1(y(:,1),y(:,2),min(x,(1-mu_d)*D_bar));
%% Loan Supply and Demand at Lowest R_L
Dlow = fzero(@(x) xi_func(r_reserves).*(psi_func(x)./((1-mu_d).*(1-dpsi_func(x).*x/Nb)) )-1,[0,10]);
LS_low = (1-rq_ratio)*Dlow*psi_func(Dlow);
LD_high = df_inv(r_reserves);
%% Loan Supply and Demand at Highest R_L
Dbar = fzero(@(x) xi_func(rho_bar).*(psi_func(x)./((1-mu_d).*(1-dpsi_func(x).*x/Nb)) )-1,[0,D_bar-1e-3]);
LS_bar = (1-rq_ratio)*Dbar*psi_func(Dbar);
LD_low = df_inv(rho_bar);
%% Loan Supply and Demand at R_L
rho_func = @(x) df(x)-1; 
eq_D_func = @(x) xi_func(rho_func(psi_func(x).*x*(1-rq_ratio))).*(psi_func(x)./((1-mu_d).*(1-dpsi_func(x).*x/Nb)) )-1;

if LS_low>LD_high
    eq = [z_func(Dlow),Dlow,r_reserves];
elseif LS_bar<LD_low
    eq = [z_func(Dbar),Dbar,rho_bar];
else
    D = fzero(@(x) eq_D_func(x),[Dlow,Dbar]);
    eq = [z_func(D),D,rho_func(D*psi_func(D)*(1-rq_ratio))];
end