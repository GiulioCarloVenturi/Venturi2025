function [c,y] = eqZD_func_con(x,i_m,i_D,lambda_func,alpha1,alpha2,alpha3,mu_d,mu_c)
y(1)=alpha1*lambda_func((1-mu_c)*x(1))+alpha3*lambda_func((1-mu_c)*x(1)+(1-mu_d)*x(2))-i_m;
y(2)=alpha2*lambda_func((1-mu_d)*x(2))+alpha3*lambda_func((1-mu_c)*x(1)+(1-mu_d)*x(2))-i_D;
c=[];
end