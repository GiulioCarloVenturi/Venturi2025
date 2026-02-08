clear;
clc;

%% Import Data
Deposit = readtimetable("Data_Calibration_GER.xlsx", "UseExcel", false,Sheet="Deposit Rate");
Loan = readtimetable("Data_Calibration_GER.xlsx", "UseExcel", false,Sheet="Loan Rate Original");
Inflation = readtimetable("Data_Calibration_GER.xlsx", "UseExcel", false,Sheet="Inflation");
M1 = readtimetable("Data_Calibration_GER.xlsx", "UseExcel", false,Sheet="M1");
GDP = readtimetable("Data_Calibration_GER.xlsx", "UseExcel", false,Sheet="GDP");
Yield = readtimetable("Data_Calibration_GER.xlsx", "UseExcel", false,Sheet="Yield5M");
MRR_rate = fillmissing(readtimetable("Data_Calibration_GER.xlsx", "UseExcel", false,Sheet="MRR"),'previous');
%% Aggregate Data on Annual Basis
Deposit_Annual = convert2annual(Deposit,"Aggregation","mean");
Loan_Annual = convert2annual(Loan,"Aggregation","mean");
Inflation_Annual = convert2annual(Inflation,"Aggregation","mean");
M1_Annual = convert2annual(M1,"Aggregation","mean");
GDP_Annual = convert2annual(GDP,"Aggregation","mean");
Yield_Annual = convert2annual(Yield,"Aggregation","mean");
MRR_rate_annual = convert2annual(MRR_rate,"Aggregation","mean");
%% Build Dataset
Data_TT = synchronize(Deposit_Annual,Loan_Annual,Inflation_Annual,M1_Annual,GDP_Annual,Yield_Annual,MRR_rate_annual,'intersection');
M1GDP = timetable(Data_TT.Dates,Data_TT.M1./Data_TT.GDP);
M1GDP.Properties.VariableNames = {'M1GDP'};
Data_TT = [Data_TT,M1GDP];
save Calibration_GER_nf.mat Data_TT
first_year = max(2000,min(year(Data_TT.Dates)));
last_year = 2012;
ID_firstyear = find(year(Data_TT.Dates) == first_year);
ID_lastyear = find(year(Data_TT.Dates) == last_year);
Data_TT = Data_TT(ID_firstyear:ID_lastyear,:);
Time = year(Data_TT.Dates);
data = [Data_TT.Yield5M/100,Data_TT.M1GDP,Time,Data_TT.DepositRate/100,Data_TT.LoanRate/100,Data_TT.Inflation/100];
Variables = Data_TT.Properties.VariableNames;
Reserves_Requirement = timetable(Data_TT.Dates,ones(size(Time,1),1) + ones(size(Time,1),1).*(Time < 2012));
%% Calibrate eta
Commercial_Loan = readtimetable("Data_Calibration_GER.xlsx",Sheet="Commercial Loan Rate");
left = 1;
right = find(year(Commercial_Loan.Dates) == 2020 & month(Commercial_Loan.Dates) == 4);
Real_Loans = Commercial_Loan.LogRealLoans(left:right);
Real_Rate = Commercial_Loan.LogRealRate(left:right);
lm = fitlm(Real_Rate,Real_Loans);
elasticity = lm.Coefficients.Estimate(2);
eta = 1 + 1/elasticity;
%% Define Parameters
rq_ratio = 0.02; 
pi_m = mean(Data_TT.Inflation)/100; 
i_l = mean(Data_TT.LoanRate)/100;
i_D = mean(Data_TT.DepositRate)/100;
i_reserves = mean(Data_TT.MRRRate)/100;
alpha1 = 0.7083*0.08 + 0.2917*0.00; 
alpha2 = 0.7083*0.0427 + 0.2917*1.00; 
alpha3 = 0.7083*0.8773 + 0.2917*0.00; 
epsilon = 0.001; 
beta = 0.99;
markup = round(1/(1-0.19),2); 
r_reserves = (1+i_reserves)/(1+pi_m)-1;
IF = 0.0037;
mu_d = 0.0092+IF;
mu_c = 0.0178; 
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
otherpar.mu_d = mu_d;
otherpar.mu_c = mu_c;
otherpar.mu_db = IF;
%% Calibration
tradep = 1;
theta_low = 0.7;
theta_cand = fzero(@(x) calibrate_money_demand_fuc(x,tradep,otherpar,data),[theta_low,1]);
[yyy, output,pars, model] = calibrate_money_demand_fuc(theta_cand,tradep,otherpar,data);
WW = [pars,output,tradep];
fileID = fopen('output.txt','w');
fprintf(fileID,'%12.8f ',WW);
fprintf(fileID,'\n');
fclose(fileID);

parameters_GER_f = importdata('output.txt');
save Calibration_GER_f.mat parameters_GER_f Data_TT eta