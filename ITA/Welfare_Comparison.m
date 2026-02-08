clear;
clc;

load simulation_nofees.mat welfare_banker welfare_e welfare_buyer_comp welfare_seller_comp
entrepreneur_nofee = welfare_e;
banker_nofee = welfare_banker;
seller_nofee = welfare_seller_comp;
buyer_nofee = welfare_buyer_comp;
clearvars welfare_buyer_comp welfare_seller_comp welfare_e welfare_banker
load simulation_fees.mat welfare_e welfare_banker welfare_buyer_comp welfare_seller_comp icbdc_grid i_m
entrepreneur_fee = welfare_e;
banker_fee = welfare_banker;
seller_fee = welfare_seller_comp;
buyer_fee = welfare_buyer_comp;
clearvars welfare_buyer_comp welfare_seller_comp welfare_e welfare_banker
% Entrepreneurs
f = figure;
hold on
plot(icbdc_grid*100,entrepreneur_nofee./entrepreneur_nofee(1),'k','linewidth',4)
plot(icbdc_grid*100,entrepreneur_fee./entrepreneur_fee(1),'r','linewidth',4)
xlabel('$i_{CBDC}$','interpreter','latex',FontSize=20)
xlim([0,3])
xtickformat("percentage")
ax=gca;
ax.FontSize = 20;
set(f,'PaperType','A4');
set(f,"PaperOrientation","landscape")
print('ITA_welfare_entre','-dpdf','-r200','-fillpage')
% Bankers
f=figure;
hold on
plot(icbdc_grid*100,banker_nofee./banker_nofee(1),'k','linewidth',4)
plot(icbdc_grid*100,banker_fee./banker_fee(1),'r','linewidth',4)
xlabel('$i_{CBDC}$','interpreter','latex',FontSize=20)
xlim([0,3])
xtickformat("percentage")
ax=gca;
ax.FontSize = 20;
set(f,'PaperType','A4');
set(f,"PaperOrientation","landscape")
print('ITA_welfare_banker','-dpdf','-r200','-fillpage')
% Households
f=figure;
hold on
plot(icbdc_grid*100,(seller_nofee+buyer_nofee)./(seller_nofee(1)+buyer_nofee(1)),'k','linewidth',4)
plot(icbdc_grid*100,(seller_fee+buyer_fee)./(seller_fee(1)+buyer_fee(1)),'r','linewidth',4)
xlabel('$i_{CBDC}$','interpreter','latex',FontSize=20)
xlim([0,3])
xtickformat("percentage")
ax=gca;
ax.FontSize = 20;
set(f,'PaperType','A4');
set(f,"PaperOrientation","landscape")
print('ITA_welfare_hh','-dpdf','-r200','-fillpage')