clear;
clc;

load parameters_ITA_final_f.mat parameters_ITA_final_f
load parameters_ITA_final_nf.mat parameters_ITA_final_nf

param = [parameters_ITA_final_nf;parameters_ITA_final_f];

%% Production Function
figure;
A_nofee = param(1,10);
A_fee = param(2,10);
eta = param(1,9);
loan = linspace(0,5,1000);
fun_nofee = A_nofee*loan.^eta;
fun_fee = A_fee*loan.^eta;
subplot(1,3,1)
hold on
plot(loan, fun_nofee,LineWidth=2)
plot(loan,fun_fee,LineWidth=2)
xlabel('$L$','Interpreter','latex',FontSize=20)
ylabel('$AL^\eta$','Interpreter','latex','FontSize',20)
legend('Costless Payments', 'Costly Payments',Location='southeast')
title('Production Function','Interpreter','latex','FontSize',20)
%% Utility Function
% Cm Utility Function
B_nofee = param(1,2);
B_fee = param(2,2);
CM_cons = linspace(0,B_fee,1000);
CM_nofee = B_nofee*log(CM_cons);
CM_fee = B_fee*log(CM_cons);
subplot(1,3,2)
hold on
plot(CM_cons, CM_nofee,LineWidth=2)
plot(CM_cons, CM_fee,LineWidth=2)
xlabel('$x$','Interpreter','latex',FontSize=20)
ylabel('$B\log(x)$','Interpreter','latex','FontSize',20)
legend('Costless Payments', 'Costly Payments',Location='southeast')
title('CM Utility Function','Interpreter','latex','FontSize',20)
% DM Utility Function
sigma_nofee = param(1,1);
sigma_fee = param(2,1);
epsilon = 0.001;
y_star = 1 - epsilon;
DM_cons = linspace(0,y_star,100);
DM_nofee = ((DM_cons+epsilon).^(1-sigma_nofee) - epsilon.^(1-sigma_nofee))./(1-sigma_nofee);
DM_fee = ((DM_cons+epsilon).^(1-sigma_fee) - epsilon.^(1-sigma_fee))./(1-sigma_fee);
subplot(1,3,3)
hold on
plot(DM_cons(2:end), DM_nofee(2:end)./DM_nofee(2),LineWidth=2)
plot(DM_cons(2:end), DM_fee(2:end)./DM_fee(2),LineWidth=2)
xlabel('$y$','Interpreter','latex',FontSize=20)
ylabel('$\frac{(y+\varepsilon)^{1-\sigma} - \varepsilon^{1-\sigma}}{1-\sigma}$','Interpreter','latex','FontSize',20)
legend('Costless Payments', 'Costly Payments',Location='southeast')
title('DM Utility Function','Interpreter','latex','FontSize',20)