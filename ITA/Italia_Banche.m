clear;
clc;

Assets = (readtimetable("Bilancio_Italia_Banche.xlsx","Sheet","Attivo_Sintetico"));
Assets_extended = readtimetable("Bilancio_Italia_Banche.xlsx",Sheet="Attivo");
Liabilities = (readtimetable("Bilancio_Italia_Banche.xlsx","Sheet","Passivo_Sintetico"));
Inflation = readtimetable("Bilancio_Italia_Banche.xlsx","Sheet","Inflation");
Deposits = readtimetable("Bilancio_Italia_Banche.xlsx","Sheet","Depositi_Sintetico");
Loans = readtimetable("Bilancio_Italia_Banche.xlsx","Sheet","Prestiti Sintetico");
ContoEconomico = readtimetable("Bilancio_Italia_Banche.xlsx",Sheet="Conto Economico");
deposit_rate = readtimetable("Bilancio_Italia_Banche.xlsx","Sheet","Depositi_Tassi");
loan_rate = readtimetable("Bilancio_Italia_Banche.xlsx","Sheet","Prestiti_Tassi");
time = Assets.Data;
asset_label = [{'Loans'} {'Securities'} {'Shares'} {'Other'}];
liability_label = [{'Foreign Deposits'} {'National Deposits'} {'Bonds'} {'Capital and Reserves'} {'Other'}];
deposit_label = Deposits.Properties.VariableNames;
deprate_label = deposit_rate.Properties.VariableNames;
Assets_real = array2table(Assets.Variables.*Inflation.RevaluationCoeff);
Liabilities_real = array2table(Liabilities.Variables.*Inflation.RevaluationCoeff);
Deposits_real = Deposits.Variables.*Inflation.RevaluationCoeff;
Loans_real = Loans.Variables.*Inflation.RevaluationCoeff;

f = figure;
subplot(1,2,1)
area(time,Assets_real.Variables*1000000)
legend(asset_label,FontSize=15,Location="northwest",Interpreter="latex")
title('Assets',FontSize=20,Interpreter='latex')
subplot(1,2,2)
area(time,Liabilities_real.Variables*1000000)
legend(liability_label,FontSize=15,Location="northwest",Interpreter="latex")
title('Liabilities',FontSize=20,Interpreter='latex')
set(f,'PaperType','A2')
set(f,"PaperOrientation","landscape")
print('ITA_bs','-dpdf','-r200','-fillpage')

f = figure;
subplot(1,3,1)
area(time,Deposits_real(:,2:3)*1000000)
legend(deposit_label(2:3),FontSize=15,Location="northwest",Interpreter="latex")
title('Foreign EU Deposits',FontSize=20)
subplot(1,3,2)
area(time,Deposits_real(:,4:end)*1000000)
legend(deposit_label(2:3),FontSize=15,Location="northwest",Interpreter="latex")
title('National Deposits',FontSize=20)
subplot(1,3,3)
area(time,[Deposits_real(:,2)+Deposits_real(:,4),Deposits_real(:,3)+Deposits_real(:,5)]*1000000)
legend(deposit_label(2:3),FontSize=15,Location="northwest",Interpreter="latex")
title('Total Deposits',FontSize=20)
set(f,'PaperType','A2')
set(f,"PaperOrientation","landscape")
print('ITA_bs_deposits','-dpdf','-r200','-fillpage')

f = figure;
subplot(1,2,1)
plot(time,[deposit_rate.TassoCC,deposit_rate.TassoTempo],LineWidth=2)
yline(1.77,'k--')
title('Deposit Rates','FontSize',20,'Interpreter','latex')
legend({'Checkable','Savings'},FontSize=15,Interpreter="latex")
subplot(1,2,2)
plot(loan_rate.Data,mean(loan_rate.Variables,2),LineWidth=2)
title('Loan Rates','FontSize',20,'Interpreter','latex')
yline(3.66,'k--')
set(f,'PaperType','A4')
set(f,"PaperOrientation","landscape")
print('ITA_bs_rate','-dpdf','-r200','-fillpage')
%% Cost Evaluation
dep_increase = [1+1.25/100,1+1.69/100];
rate_l = [3.78/100,3.66/100];
rate_d = [0.76/100,0.90/100];
rate_l_new = [3.34/100,3.07/100];
rate_d_new = [1.43/100,1.77/100];
res_req = 0.01;
res_rate = 2.56/100;
t = find( month(Assets.Data) == 12 & year(Assets.Data) == 2007);
tt = find( year(ContoEconomico.Data) == 2007);
% New Loans
deposits_t = mean(Deposits_real(t:t+11,2)+Deposits_real(t:t+11,4))*1000000;
deposits_t_new = deposits_t.*dep_increase;
loans_t = (mean(Loans_real(t:t+11,7)) - mean(Loans_real(t+12:t+12+11,7)))*1000000; 
loans_t_new = loans_t+deposits_t.*(dep_increase-1)*(1-res_req);
reserves_new = deposits_t.*(dep_increase-1)*res_req;
int_act_delta = rate_l_new.*loans_t_new - rate_l.*loans_t + reserves_new*res_rate;
int_pass_delta = rate_d_new.*deposits_t_new-rate_d.*deposits_t;
delta_profit = int_act_delta - int_pass_delta;
disp(['Riduzione Interesse Attivi ' num2str(int_act_delta/(ContoEconomico.InteressiAttivi(tt)*1000000)*100) '%'])
disp(['Aumento Interesse Passivi ' num2str(int_pass_delta/(ContoEconomico.InteressiPassivi(tt)*1000000)*100) '%'])
disp(['Riduzione Margine Intermediazione ' num2str(delta_profit/(ContoEconomico.MargineDiIntermediazione(tt)*1000000)*100) '%'])
disp(['Riduzione Utile Netto ' num2str(delta_profit/(ContoEconomico.UtileO___PerditaNetti(tt)*1000000)*100) '%'])
disp(['Riduzione Utile Netto in mld' num2str(delta_profit/1000000000)])
% 2007 Loans
deposits_t = mean(Deposits_real(t:t+11,2)+Deposits_real(t:t+11,4))*1000000;
deposits_t_new = deposits_t.*dep_increase;
loans_t = mean(Loans_real(t:t+11,7))*1000000; 
loans_t_new = loans_t+deposits_t.*(dep_increase-1)*(1-res_req);
reserves_new = deposits_t.*(dep_increase-1)*res_req;
int_act_delta = rate_l_new.*loans_t_new - rate_l.*loans_t + reserves_new*res_rate;
int_pass_delta = rate_d_new.*deposits_t_new-rate_d.*deposits_t;
delta_profit = int_act_delta - int_pass_delta;
disp(['Riduzione Interesse Attivi ' num2str(int_act_delta/(ContoEconomico.InteressiAttivi(tt)*1000000)*100) '%'])
disp(['Aumento Interesse Passivi ' num2str(int_pass_delta/(ContoEconomico.InteressiPassivi(tt)*1000000)*100) '%'])
disp(['Riduzione Margine Intermediazione ' num2str(delta_profit/(ContoEconomico.MargineDiIntermediazione(tt)*1000000)*100) '%'])
disp(['Riduzione Utile Netto ' num2str(delta_profit/(ContoEconomico.UtileO___PerditaNetti(tt)*1000000)*100) '%'])
disp(['Riduzione Utile Netto in mld' num2str(delta_profit/1000000000)])


Z = retime(timetable(time,Assets_real.Variables),"yearly","mean");
ZZ = Z.Variables;
figure;
subplot(1,4,1)
plot(Z.time(2:end),100*(diff(ZZ(:,1))./ZZ(1:end-1,1)))
title('Prestiti')
subplot(1,4,2)
plot(Z.time(2:end),100*(diff(ZZ(:,2))./ZZ(1:end-1,2)))
title('Titoli')
subplot(1,4,3)
plot(Z.time(2:end),100*(diff(ZZ(:,3))./ZZ(1:end-1,3)))
title('Azioni')
subplot(1,4,4)
plot(Z.time(2:end),100*(diff(ZZ(:,4))./ZZ(1:end-1,4)))
title('Altro')