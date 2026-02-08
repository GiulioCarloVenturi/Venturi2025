function [yyy, output,pars, model]=calibrate_money_demand_fuc(ppp,tradep,otherpar,data)
alpha1 = otherpar.alpha1;
alpha2 = otherpar.alpha2;
alpha3 = otherpar.alpha3;
otherpar.theta = ppp; 
trading_prob = tradep; 
alpha1 = alpha1*trading_prob;
alpha2 = alpha2*trading_prob;
alpha3= alpha3*trading_prob;
x0=[0.8,4]; % x=[sigma,B]
xlow=[0,0]; 
xhigh=[5,10]; 
tic
x = fmincon(@(x)cal_MD_v2(x,data,otherpar,alpha1,alpha2,alpha3),x0,[],[],[],[],xlow,xhigh);
toc
[yy,model,output]=cal_MD_v2(x,data,otherpar,alpha1,alpha2,alpha3); 
yyy= output(1) - otherpar.markup;
pars = [yy,ppp,x];