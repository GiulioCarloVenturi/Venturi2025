function [yy,M1GDP,output]=cal_MD_v2(parameters,data,otherpar,alpha1,alpha2,alpha3)
options=optimset('tolX',1e-10,'tolFun',1e-10,'display','off');
sigma=parameters(1);
B=parameters(2);
beta = otherpar.beta;
epsilon = otherpar.epsilon;
theta= otherpar.theta;
eta = otherpar.eta;
rq_ratio = otherpar.rq_ratio;
r_reserves = otherpar.r_reserves;
u = @(x) ((x+epsilon).^(1-sigma)-epsilon.^(1-sigma))./(1-sigma); 
du = @(x)(x+epsilon).^(-sigma);
c = @(x) x;                                                  
dc = @(x) 1;
y_star = fminsearch(@(x)(du(x)-dc(x)).^2,1);                    
g_func = @(x) (1-theta).*u(x)+theta.*c(x);                     
dg_func = @(x) (1-theta).*du(x)+theta.*dc(x);
y_grid = linspace(0,y_star,10001);
p_grid = g_func(y_grid);
y_func = @(p)interp1(p_grid,y_grid,min(p,max(p_grid)));         
lambda_func = @(z) max(du(y_func(z))./dg_func(y_func(z))-1,0);     
xlow = [0,0];
xhigh = [p_grid(end),p_grid(end)];
x0 = [0.5,0.5];
Ndata = length(data);
M1GDP = zeros(Ndata,1);
x = zeros(Ndata,2);
M1 = zeros(Ndata,1);
GDP = zeros(Ndata,1);
markup = zeros(Ndata,1);
DMoutput_ratio = zeros(Ndata,1);
nominal_rho = data(:,5);
inflation = data(:,end);
for i=1:Ndata
    i_m = data(i,1);
    i_d = data(i,4);
    rho = (1+nominal_rho(i))/(1+inflation(i))-1;
    x(i,:)=fmincon(@(x) 1,x0,[],[],[],[],xlow,xhigh,@(x) ...
        eqZD_func_con(x,i_m,(1+i_m)/(1+i_d)-1,lambda_func,alpha1,alpha2,alpha3),options);
    eq = x(i,:);
    psi = (1+i_m)/(1+i_d)*beta;
    M1(i) = x(i,1)+x(i,2)*psi;
    GDP(i) = B+alpha1*x(i,1)+alpha2*x(i,2)+alpha3*min([sum(x(i,:)),p_grid(end)])...
        +(1+r_reserves)*rq_ratio*x(i,2)*psi-x(i,2)+(1+rho+eta)*(1-rq_ratio)*...
        x(i,2)*psi/eta;
    M1GDP(i) = M1(i)/GDP(i);
    x0 = x(i,:);
    markup(i) = (alpha1*min(eq(1),max(p_grid))/y_func(eq(1)) +...
        alpha2*min(eq(2),max(p_grid))/y_func(eq(2)) ...
        +alpha3*min(eq(2)+eq(1),max(p_grid))/y_func(eq(2)+eq(1)))/(alpha1+alpha2+alpha3);
    DMoutput_ratio(i) = (alpha1*eq(1)+alpha2*eq(2)+...
        alpha3*min([sum(eq),p_grid(end)]))/GDP(i);
end
yy = sum((M1GDP(:)-data(:,2)).^2);
M1GDP=M1GDP(:);
output = [mean(markup), mean(DMoutput_ratio)];