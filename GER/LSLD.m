%% Compute Loan Demand and Supply
function [LS, LD, rho_grid]=LSLD(rho_grid,func_var,Nb,y)
mu_d = func_var.mu_d;
rq_ratio = func_var.rq_ratio;
cost = func_var.cost;
beta = func_var.beta;
df_inv = func_var.df_inv;
r_reserves = func_var.r_reserves;
xi_func = @(rho) (1-rq_ratio)*max(1+rho,1+r_reserves)+rq_ratio.*(1+r_reserves)-cost;
D_bar = max(y(:,1));
psi_func = @(x) interp1(y(:,1),y(:,3),min(x,D_bar)).*(x<=D_bar)+beta.*(x>D_bar);
dpsi_func = @(x) interp1(y(:,1),y(:,4),min(x,D_bar)).*(x<=D_bar);
eq_D_func = @(x,y) xi_func(y).*(psi_func(x)./((1-mu_d).*(1-dpsi_func(x).*x/Nb)) )-1;

%% Compute Loan Supply for equation 12 case 3: R_r < R_l < 1/beta
if rho_grid(1) < r_reserves
    ID = find(rho_grid > r_reserves,1,'first');
    rho_grid = [rho_grid(1:ID-1),r_reserves,rho_grid(ID:end)];
end
for i=1:length(rho_grid)
    [D_grid(i),fval(i)] = fzero(@(x) eq_D_func(x,rho_grid(i)), [1e-4,D_bar]);
    LS(i) = (1-rq_ratio)*D_grid(i)*psi_func(D_grid(i));
end
%% Compute Loan Supply for equation 12 case 3: R_l < R_r
if rho_grid(1) < r_reserves
    ID = find(rho_grid < r_reserves,1,'last');
    LS(1:ID)=0;
end
%% Compute Loan Demand equation, page 17
LD = df_inv(rho_grid);