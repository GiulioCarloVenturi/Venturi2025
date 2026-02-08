%% Compute inverse demand function for deposits without CBDC 

function y = D_inv_demand(i_m,D_grid,func_var)
%% Load Parameters
mu_d = func_var.mu_d;
D_grid = D_grid(:); 
alpha1 = func_var.alpha1; 
alpha2 = func_var.alpha2; 
alpha3 = func_var.alpha3;
lambda_func = func_var.lambda_func; 
dlambda_func = func_var.dlambda_func;
Dgrid_size = length(D_grid);
p_star_d = func_var.p_star_d;
beta = func_var.beta;
z = zeros(Dgrid_size,1);
%% Compute Cash balances for each amount of Deposits in Dgrid, equation 7
for i =1:Dgrid_size
    if alpha1*lambda_func(0)+alpha3*lambda_func(0+(1-mu_d)*D_grid(i))<i_m
       z(i)=0;
    else
    if i<=5
         z(i) = fzero(@(z) alpha1*lambda_func(z) + alpha3*lambda_func(z+(1-mu_d)*D_grid(i)) - i_m,[0,p_star_d]);
    else
        z(i) = fzero(@(z) alpha1*lambda_func(z) + alpha3*lambda_func(z+(1-mu_d)*D_grid(i)) - i_m,[0,z(i-1)+0.01]);
    end
    end
end
%% 1/R_d
psi = beta*(1-mu_d)*(1+ alpha2*lambda_func((1-mu_d)*D_grid) +...
    alpha3*lambda_func(z+(1-mu_d)*D_grid));
%% Derivative of 1/R_d
% Differential of line 25
dzdD=-alpha3*dlambda_func(z+(1-mu_d)*D_grid)*(1-mu_d)./...
    (alpha1*dlambda_func(z) + alpha3*dlambda_func(z+(1-mu_d)*D_grid));
dpsi = beta*(1-mu_d)*(alpha2*dlambda_func((1-mu_d)*D_grid)*(1-mu_d) + alpha3*dlambda_func(z+(1-mu_d)*D_grid).*(1-mu_d+dzdD));
%% Output
y = [D_grid,z,psi,dpsi];