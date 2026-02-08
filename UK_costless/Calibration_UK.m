clear;
clc;

%% Import Data
Deposit = readtimetable("Data_Calibration_UK.xlsx", "UseExcel", false,Sheet="Deposit Rate");
Loan = readtimetable("Data_Calibration_UK.xlsx", "UseExcel", false,Sheet="Loan Rate");
Inflation = readtimetable("Data_Calibration_UK.xlsx", "UseExcel", false,Sheet="Inflation");
M1 = readtimetable("Data_Calibration_UK.xlsx", "UseExcel", false,Sheet="M1");
GDP = readtimetable("Data_Calibration_UK.xlsx", "UseExcel", false,Sheet="GDP");
Yield = readtimetable("Data_Calibration_UK.xlsx", "UseExcel", false,Sheet="Yield5M");
MRR_rate = fillmissing(readtimetable("Data_Calibration_UK.xlsx", "UseExcel", false,Sheet="Bank Rate"),'previous');
%% Aggregate Data on Annual Basis
Deposit_Annual = convert2annual(Deposit,"Aggregation","mean");
Loan_Annual = convert2annual(Loan,"Aggregation","mean");
Inflation_Annual = convert2annual(Inflation,"Aggregation","mean");
M1_Annual = convert2annual(M1,"Aggregation","mean");
GDP_Annual = convert2annual(GDP,"Aggregation","sum");
Yield_Annual = convert2annual(Yield,"Aggregation","mean");
MRR_rate_annual = convert2annual(MRR_rate,"Aggregation","mean");
%% Build Dataset
Data_TT = synchronize(Deposit_Annual,Loan_Annual,Inflation_Annual,M1_Annual,GDP_Annual,Yield_Annual,MRR_rate_annual,'intersection');
M1GDP = timetable(Data_TT.Dates,Data_TT.M1./Data_TT.GDP);
M1GDP.Properties.VariableNames = {'M1GDP'};
Data_TT = [Data_TT,M1GDP];
save Calibration_UK_nf.mat Data_TT
first_year = 2000;
last_year = 2012;
ID_firstyear = find(year(Data_TT.Dates) == first_year);
ID_lastyear = find(year(Data_TT.Dates) == last_year);
Data_TT = Data_TT(1:ID_lastyear,:);
Time = year(Data_TT.Dates);
data = [Data_TT.Yield5M/100,Data_TT.M1GDP,Time,Data_TT.DepositRate/100,Data_TT.LoanRate/100,Data_TT.Inflation];
Variables = Data_TT.Properties.VariableNames;
%% Calibrate eta 1999-2020
Commercial_Loan = rmmissing(readtimetable("Data_Calibration_UK.xlsx",Sheet="Commercial Loan Rate"));
left = 1;
right = find(year(Commercial_Loan.Dates) == 2020 & month(Commercial_Loan.Dates) == 2); %%
Real_Loans = Commercial_Loan.LogRealLoans(left:right);
Real_Rate = Commercial_Loan.LogRealRate(left:right);
lm = fitlm(Real_Rate,Real_Loans);
elasticity = lm.Coefficients.Estimate(2);
eta = 1 + 1/elasticity;
%% Reserve Ratio
CB_Cost = rmmissing(readtimetable("Data_Calibration_UK.xlsx", "UseExcel", false,Sheet="CB Cost"));
Total_Deposits = convert2annual(rmmissing(readtimetable("Data_Calibration_UK.xlsx", "UseExcel", false,Sheet="Deposits")),"Aggregation","mean");
Total_Assets = convert2annual(rmmissing(readtimetable("Data_Calibration_UK.xlsx", "UseExcel", false,Sheet="Assets")),"Aggregation","mean");
Cost_TT = synchronize(CB_Cost,Total_Deposits,Total_Assets,'intersection');
rq_ratio = mean((Cost_TT.Costs)./Cost_TT.Deposits); %%%
rq_ratio = 0.00;
%% Define Parameters
theta_low = 0.90;
pi_m = mean(Data_TT.Inflation); 
i_l = mean(Data_TT.LoanRate)/100;
i_D = mean(Data_TT.DepositRate)/100;
i_reserves = 0.00;
alpha1 = 0.75*0.02 + 0.25*0.00; %%%
alpha2 = 0.75*0.12 + 0.25*1.00; %%%
alpha3 = 0.75*0.86 + 0.25*0.00; %%%
epsilon = 0.001; 
beta = 0.99;
markup = 1.20; %%%
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
parameters_UK_nf = [pars,output,tradep];

save Calibration_UK_nf.mat parameters_UK_nf Data_TT eta rq_ratio