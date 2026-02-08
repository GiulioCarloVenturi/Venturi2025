clear;
clc;

load simulation_fees.mat
eq_fee = eq_cbdc;
output_fee = output_index;
deposit_fee = psi_D;
clearvars -except eq_fee output_fee deposit_fee
load simulation_nofees.mat
eq_nofee = eq_cbdc;
output_nofee = output_index;
deposit_nofee = psi_D;
clearvars -except eq_fee output_fee deposit_fee eq_nofee output_nofee deposit_nofee icbdc_grid pi_m

%% Plot Macro Variables
f = figure;
% Deposit Rate
hold on
plot(icbdc_grid*100,((1+eq_fee(:,4))*(1+pi_m)-1)*100,'r','linewidth',4)
plot(icbdc_grid*100,((1+eq_nofee(:,4))*(1+pi_m)-1)*100,'k','linewidth',4)
xlim([0,3])
xlabel('$i_{CBDC}$', 'interpreter','latex')
ytickformat("percentage");
xtickformat("percentage");
ax=gca;
ax.FontSize = 20;
set(f,'PaperType','A4');
set(f,"PaperOrientation","landscape")
print('ITA_deprate','-dpdf','-r200','-fillpage')
% Deposits
f = figure;
hold on
plot(icbdc_grid*100,(deposit_fee(:,1)./deposit_fee(1,1)-1)*100,'r','linewidth',4);
plot(icbdc_grid*100,(deposit_nofee(:,1)./deposit_nofee(1,1)-1)*100,'k','linewidth',4);
plot(icbdc_grid*100,zeros(size(icbdc_grid,2)),'k--')
xlim([0,3])
xlabel('$i_{CBDC}$', 'interpreter','latex')
xtickformat("percentage");
ytickformat("percentage");
ax=gca;
ax.FontSize = 20;
set(f,'PaperType','A4');
set(f,"PaperOrientation","landscape")
print('ITA_deposit','-dpdf','-r200','-fillpage')
% Loan Rate
f=figure;
hold on
plot(icbdc_grid*100,((1+eq_fee(:,3))*(1+pi_m)-1)*100,'r','linewidth',4);
plot(icbdc_grid*100,((1+eq_nofee(:,3))*(1+pi_m)-1)*100,'k','linewidth',4);
xlim([0,3])
xlabel('$i_{CBDC}$', 'interpreter','latex')
ytickformat("percentage");
xtickformat("percentage");
ax=gca;
ax.FontSize = 20;
set(f,'PaperType','A4');
set(f,"PaperOrientation","landscape")
print('ITA_loanrate','-dpdf','-r200','-fillpage')
% Loans
f=figure;
hold on
plot(icbdc_grid*100,(eq_fee(:,end)./eq_fee(1,end)-1)*100,'r','linewidth',4)
plot(icbdc_grid*100,(eq_nofee(:,end)./eq_nofee(1,end)-1)*100,'k','linewidth',4)
plot(icbdc_grid*100,zeros(size(icbdc_grid,2)),'k--');
xlim([0,3])
xlabel('$i_{CBDC}$', 'interpreter','latex')
xtickformat("percentage");
ytickformat("percentage");
ax=gca;
ax.FontSize = 20;
set(f,'PaperType','A4');
set(f,"PaperOrientation","landscape")
print('ITA_loan','-dpdf','-r200','-fillpage')
% Loan-Deposit Spread
ndrate_fee=((1+eq_fee(:,4))*(1+pi_m)-1)*100;
nlrate_fee=((1+eq_fee(:,3))*(1+pi_m)-1)*100;
ndrate_nofee=((1+eq_nofee(:,4))*(1+pi_m)-1)*100;
nlrate_nofee=((1+eq_nofee(:,3))*(1+pi_m)-1)*100;
f=figure;
hold on
plot(icbdc_grid*100,nlrate_fee-ndrate_fee,'r','linewidth',4)
plot(icbdc_grid*100,nlrate_nofee-ndrate_nofee,'k','linewidth',4)
xlim([0,3])
ylim([1,3.5])
xlabel('$i_{CBDC}$', 'interpreter','latex')
ytickformat("percentage");
xtickformat("percentage");
ax=gca;
ax.FontSize = 20;
set(f,'PaperType','A4');
set(f,"PaperOrientation","landscape")
print('ITA_spread','-dpdf','-r200','-fillpage')
% Output
f=figure;
hold on
plot(icbdc_grid*100,output_fee(:,1),'r','linewidth',4)
plot(icbdc_grid*100,output_nofee(:,1),'k','linewidth',4)
plot(icbdc_grid*100,zeros(size(icbdc_grid,2)),'k--')
xlim([0,3])
xlabel('$i_{CBDC}$', 'interpreter','latex')
ytickformat("percentage");
xtickformat("percentage");
ax=gca;
ax.FontSize = 20;
set(f,'PaperType','A4');
set(f,"PaperOrientation","landscape")
print('ITA_GDP','-dpdf','-r200','-fillpage')
%%
% No fees
[loan_max_nofee,nmax_nofee] = max(eq_nofee(:,end)./eq_nofee(1,end));
i_max_nofee=icbdc_grid(nmax_nofee);
id_nofee = find(eq_nofee(:,end)./eq_nofee(1,end)>=1,1,'last');
id1_nofee = find(eq_nofee(:,end)./eq_nofee(1,end)>1,1,'first');
i_max1_nofee=icbdc_grid(id_nofee);
i_min1_nofee=icbdc_grid(id1_nofee);
[output_max_nofee,nmax_nofee] = max(output_nofee(:,1));
ioutput_max_nofee=icbdc_grid(nmax_nofee);
id_nofee = find(output_nofee(:,1)>=0,1,'last');
id1_nofee = find(output_nofee(:,1)>0,1,'first');
iy_max1_nofee=icbdc_grid(id_nofee);
iy_min1_nofee=icbdc_grid(id1_nofee);
% Fees
[loan_max_fee,nmax_fee] = max(eq_fee(:,end)./eq_fee(1,end));
i_max_fee=icbdc_grid(nmax_fee);
id_fee = find(eq_fee(:,end)./eq_fee(1,end)>=1,1,'last');
id1_fee = find(eq_fee(:,end)./eq_fee(1,end)>1,1,'first');
i_max1_fee=icbdc_grid(id_fee);
i_min1_fee=icbdc_grid(id1_fee);
[output_max_fee,nmax_fee] = max(output_fee(:,1));
ioutput_max_fee=icbdc_grid(nmax_fee);
id_fee = find(output_fee(:,1)>=0,1,'last');
id1_fee = find(output_fee(:,1)>0,1,'first');
iy_max1_fee=icbdc_grid(id_fee);
iy_min1_fee=icbdc_grid(id1_fee);
%
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('0% fee on transaction settled with deposits');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(['CBDC increases Bank Lending if its rate is between ' num2str(i_min1_nofee*100) '% and ' num2str(i_max1_nofee*100) '%']);
disp(['It increases output if its rate is between ' num2str(iy_min1_nofee*100) '% and ' num2str(iy_max1_nofee*100) '%']);
disp(['The CBDC rate that maximizes lending is ' num2str(i_max_nofee*100) '% and that maximizes output is ' num2str(ioutput_max_nofee*100) '%' ]);
disp(['The maximum increase in lending is ' num2str((loan_max_nofee-1)*100) '% and the maximum increase in output is ' num2str(output_max_nofee) '%' ]);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('2.53% (0.73% merchant fee) fee on deposits, 4.70% fee on cash, 1.30% fee on CBDC');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(['CBDC increases Bank Lending if its rate is between ' num2str(i_min1_fee*100) '% and ' num2str(i_max1_fee*100) '%']);
disp(['It increases output if its rate is between ' num2str(iy_min1_fee*100) '% and ' num2str(iy_max1_fee*100) '%']);
disp(['The CBDC rate that maximizes lending is ' num2str(i_max_fee*100) '% and that maximizes output is ' num2str(ioutput_max_fee*100) '%' ]);
disp(['The maximum increase in lending is ' num2str((loan_max_fee-1)*100) '% and the maximum increase in output is ' num2str(output_max_fee) '%' ]);