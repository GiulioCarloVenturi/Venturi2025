clear;
clc;

%% Import Data
Deposit = readtimetable("Data_Calibration_US.xlsx", "UseExcel", false,Sheet="Deposit Rate");
Loan = readtimetable("Data_Calibration_US.xlsx", "UseExcel", false,Sheet="Loan Rate");
Inflation = readtimetable("Data_Calibration_US.xlsx", "UseExcel", false,Sheet="Inflation");
M1 = readtimetable("Data_Calibration_US.xlsx", "UseExcel", false,Sheet="M1");
GDP = readtimetable("Data_Calibration_US.xlsx", "UseExcel", false,Sheet="GDP");
Yield = readtimetable("Data_Calibration_US.xlsx", "UseExcel", false,Sheet="Yield3M");
%% Aggregate Data on Annual Basis
Deposit_Annual = convert2annual(Deposit,"Aggregation","mean");
Loan_Annual = convert2annual(Loan,"Aggregation","mean");
Inflation_Annual = convert2annual(Inflation,"Aggregation","mean");
M1_Annual = convert2annual(M1,"Aggregation","mean");
GDP_Annual = convert2annual(GDP,"Aggregation","mean");
Yield_Annual = convert2annual(Yield,"Aggregation","mean");
M1GDP = readtimetable("Data_Calibration_US.xlsx","UseExcel",false,Sheet="M1GDP");
%% Build Dataset
Data_TT = synchronize(Deposit_Annual,Loan_Annual,Inflation_Annual,M1_Annual,GDP_Annual,Yield_Annual,M1GDP,'intersection');
save Calibration_US_nf.mat Data_TT
first_year = 2000;
last_year = 2012;
ID_firstyear = find(year(Data_TT.Dates) == first_year);
ID_lastyear = find(year(Data_TT.Dates) == last_year);
Data_TT = Data_TT(ID_firstyear:ID_lastyear,:);
Time = year(Data_TT.Dates);
data = [Data_TT.Yield3M/100,Data_TT.M1GDP,Time,Data_TT.DepositRate,Data_TT.LoanRate,Data_TT.Inflation/100];
Variables = Data_TT.Properties.VariableNames;
%% Calibrate eta
Commercial_Loan = rmmissing(readtimetable("Data_Calibration_US.xlsx",Sheet="Commercial Loan Rate"));
left = find(year(Commercial_Loan.Dates) == 1995 & month(Commercial_Loan.Dates) == 1);
right = find(year(Commercial_Loan.Dates) == 2020 & month(Commercial_Loan.Dates) == 2);
Real_Loans = Commercial_Loan.LogRealLoans(left:right);
Real_Rate = Commercial_Loan.LogRealRate(left:right);
lm = fitlm(Real_Rate,Real_Loans);
elasticity = lm.Coefficients.Estimate(2);
eta = 1 + 1/elasticity;
%% Define Parameters
theta_low = 0.50;
rq_ratio = 0.02; 
pi_m = mean(Data_TT.Inflation)/100; 
i_l = mean(Data_TT.LoanRate);
i_D = mean(Data_TT.DepositRate);
i_reserves = 0.00; %mean(Data_TT.MRRRate)/100; 
alpha1 = 0.0458; 
alpha2 = 0.2644; 
alpha3 = 0.6898; 
epsilon = 0.001; 
beta = 0.99;
markup = 1.20; 
r_reserves = (1+i_reserves)/(1+pi_m)-1;
otherpar.markup = markup; 
otherpar.r_reserves = r_reserves; 
otherpar.beta = beta;
otherpar.epsilon = epsilon;
otherpar.pi_m = pi_m;
otherpar.spread = i_l-i_D;
otherpar.i_D = i_D;
otherpar.i_l = i_l;
otherpar.alpha1 = alpha1;
otherpar.alpha2 = alpha2;
otherpar.alpha3 = alpha3;
otherpar.rq_ratio = rq_ratio;
otherpar.eta = eta; 
%% Calibration
tradep = 1;
theta_cand = fzero(@(x) calibrate_money_demand_fuc(x,tradep,otherpar,data),[theta_low,1]);
[yyy, output,pars, model] = calibrate_money_demand_fuc(theta_cand,tradep,otherpar,data);
WW = [pars,output,tradep];
fileID = fopen('output.txt','w');
fprintf(fileID,'%12.8f ',WW);
fprintf(fileID,'\n');
fclose(fileID);

parameters_US_nf = importdata('output.txt');
save Calibration_US_nf.mat parameters_US_nf Data_TT