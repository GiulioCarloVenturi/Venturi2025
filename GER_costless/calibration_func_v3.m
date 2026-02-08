function [drate_nominal_diff,yy,eq_final,d2c]=calibration_func_v3(A,otherpar,psiz_grid)

beta = otherpar.beta;
pi_m = otherpar.pi_m; 
alpha2 = otherpar.alpha2;
alpha3 = otherpar.alpha3;
spread = otherpar.spread;
mu_d = otherpar.mu_d;
eta = otherpar.eta;
df = @(k) A*eta*k.^(eta-1); 
df_inv = @(rho) (A*eta./(1+rho)).^(1/(1-eta));
otherpar.df_inv = df_inv;
otherpar.df = df;
Nb_grid = 1:1000;
lambda_func = otherpar.lambda_func;
for i=1:length(Nb_grid)
    Nb = Nb_grid(i);
    y = SS_eq_noCBDC(otherpar,Nb,psiz_grid );
    eq_final(i,:) = y;
    psi = (1-mu_d)*beta*(1 + alpha2*lambda_func((1-mu_d)*eq_final(i,2)) + alpha3*lambda_func(eq_final(i,1)+(1-mu_d)*eq_final(i,2)));
    drate(i) = 1./psi - 1;
    if eq_final(i,3) - drate(i) < spread
            yy = Nb;
        break;
    end
end
d2c = eq_final(yy,2)*psi/eq_final(yy,1); % Deposit-to-Cash Ratio
drate_nominal_diff = (1+drate(yy))*(1+pi_m)- (1+otherpar.i_D);