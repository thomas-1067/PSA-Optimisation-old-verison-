clear,clc
%% MAIN
% Define bounds
addpath('NSGA-II')

% Setup parallel pool if needed
if isempty(gcp('nocreate'))
    parpool('local', 12); % Adjust number of workers based on your system
end


% Define your material properties and type here, if they're not variable
c = 10; % number of cycles
N = 10; % number of nodes


% x(1) = adsorptionTime;
% x(3) = adsorptionFeedVelocity;
% x(4) = purgeFeedVelocity;
% x(5) = adsorptionPressure;

% Define the NSGA-II objective function
Function = @(x)PSA_Process(x,N,c);


% NSGA-II options (adjust these as needed)
options         = nsgaopt(); % Create default options structure
options.popsize = 40; % Population size
options.maxGen  = 3; % Max generation

options.vartype    = [1, 1, 1, 1]         ;
options.outputfile = 'UTSA-16_Process.txt' ;

options.numObj  = 2; % Number of objectives
options.numVar  = 4; % Number of design variables (adjust if different)
options.numCons = 0; % Number of constraints (adjust based on your function)
%%              AdT/PurT  Adv  PurV     AdP
options.lb      = [ 10, 0.05, 0.05,    5.5e5]; % Lower bound of variables
options.ub      = [ 500,  1,    1,      10e5]; % Upper bound of variables

options.nameObj = {'-Purity_H2', 'Recovery_H2'} ;% the objective names are showed in GUI window.
options.objfun = Function; % Objective function handle

options.useParallel = 'yes'; % Use parallel computation
options.poolsize     = 12   ;                            % number of worker processes
%% Run NSGA-II
result = nsga2(options);
%%

% Accessing variable inputs for all individuals
variableInputs = [result.pops.var]; % Assuming 'result' is your struct array

% Accessing objective function outputs for all individuals
objectiveOutputs = [result.pops.obj]; % This will concatenate the 'obj' fields of all structs in an array

rankOneIndices = find([result.pops.rank] == 1);

% Initialize an array to hold the objective pairs for rank 1 individuals
rankOneObjectivePairs = zeros(length(rankOneIndices), 2);
rankOneVariables = zeros(length(rankOneIndices), 4);
% Loop over each rank 1 individual to extract their objective pairs
for i = 1:length(rankOneIndices)
    % Current rank 1 individual index
    idx = rankOneIndices(i);
    
    % Extract objectives for the current individual, assuming 'obj' field exists
    objectives = result.pops(idx).obj;
    variablesout = result.pops(idx).var;
    % Since objectives are ordered in columns (odd indices for purity, even for recovery),
    % we can directly assign them to our array
    rankOneObjectivePairs(i, 1) = -objectives(1); % Purity
    rankOneObjectivePairs(i, 2) = -objectives(2); % Recovery
    
    rankOneVariables(i, 1) = variablesout(1);
    rankOneVariables(i, 2) = variablesout(2);
    rankOneVariables(i, 3) = variablesout(3);
    rankOneVariables(i, 4) = variablesout(4)/1e5;
end

disp(objectives)
input_var = result.pops(idx).var;
disp(input_var)

% Plotting the Pareto front
figure(1);
set(gcf, 'Position', [90, 90, 500, 300], 'Color', 'w');
scatter(rankOneObjectivePairs(:,1), rankOneObjectivePairs(:,2), 'bo','MarkerFaceColor', 'b');
xlabel('Productivity (mol_{H2out}/mol_{in}s)');
ylabel('Energy (kWh/mol_{H2out})');
ax = gca; % Current axis


ax.LineWidth = 1.2; % Making axes lines thicker
ax.FontSize = 12; % Increasing font size for better readability
axis tight;
% Add a box around the plot for better framing
box on;

function  [objectives, con] = PSA_Process(x,N,c)
% Physical Parameters
con = 0;
% [Purity_CH4, Purity_H2, Recovery_CH4, Recovery_H2, Work, MB]
L   = 3;            % Column length m 
t = 0;

adsorptionTime =x(1);
adsorptionFeedVelocity = x(2);
purgeFeedVelocity = x(3) ;
adsorptionPressure = x(4) ;
% make sure time for adsortion is the same as time for desorption
steptimes   = [20 adsorptionTime 20 adsorptionTime];     % Time vector
stepend = 4;
cycleend = c; % max cycle limit
n   = N + 1;
dz  = L/N;          % Mesh size
z   = 0:dz:L;       % Mesh generation



%%%%%%% blended gas parameters %%%%%%%%%%%
TPF = [adsorptionFeedVelocity 300 adsorptionPressure 1e5 purgeFeedVelocity];   % Feed: Velocity (m/s), tempurature (K) and pressure (Bar) [Vfeed Tfeed PH PI PL] 
yF  = 0.80;              % Mole fraction for methane 

% y is an array size n*5 of y1 = 1:n, q1 = n+1:2*n,  
% q2 = 2*n+1:3*n, T = 3*n+1:4*n, P = 4*n+1:5*n


[solPr, solAd, solDeP, solPur, tim, yadoot, phigh, plow] = adsorptionSolver(t,z,yF,TPF,steptimes,stepend,cycleend);

% for i = 1:size(sol.y,2)
%   [~,velocity(:,i),~,~,~] = adsorptionModel(sol.x(i),sol.y(:,i),z(2)-z(1),numel(z),yF,TPF,stepNo); 
% end

for i = 1:size(solPr.y,2)
    stepNo = 1;
    
  [~,Prvhalf(:,i),~,~] = adsorptionModel(solPr.x(i),solPr.y(:,i),z(2)-z(1),numel(z),yF,TPF,stepNo,yadoot,tim ,phigh, plow);  
end

for i = 1:size(solAd.y,2)
    stepNo = 2;
  [~,Advhalf(:,i),~,~] = adsorptionModel(solAd.x(i),solAd.y(:,i),z(2)-z(1),numel(z),yF,TPF,stepNo,yadoot,tim ,phigh, plow);  
end

for i = 1:size(solDeP.y,2)
    stepNo = 3;
  [~,DePvhalf(:,i),~,~] = adsorptionModel(solDeP.x(i),solDeP.y(:,i),z(2)-z(1),numel(z),yF,TPF,stepNo,yadoot,tim ,phigh, plow);  
end

for i = 1:size(solPur.y,2)
    stepNo=4;
  [~,Purvhalf(:,i),~,~] = adsorptionModel(solPur.x(i),solPur.y(:,i),z(2)-z(1),numel(z),yF,TPF,stepNo,yadoot,tim,phigh, plow);  
end



% sol.x is time steps
% sol.y is 

for assign = 1
% sol.x is time steps
% sol.y is 
MoleyPr2 = 1-solPr.y(1:n,1:end);
MoleyAd2 = 1-solAd.y(1:n,1:end);
MoleyDeP2 = 1-solDeP.y(1:n,1:end);
MoleyPur2 = 1-solPur.y(1:n,1:end);


y1_pr = solPr.y(1:n,1:end);
q1_pr = solPr.y(n+1:2*n,1:end);
q2_pr = solPr.y(2*n+1:3*n,1:end);
T_pr = solPr.y(3*n+1:4*n,1:end);
P_pr = solPr.y(4*n+1:5*n,1:end);
time_pr = solPr.x;

y1_ad = solAd.y(1:n,1:end);
q1_ad = solAd.y(n+1:2*n,1:end);
q2_ad = solAd.y(2*n+1:3*n,1:end);
T_ad = solAd.y(3*n+1:4*n,1:end);
P_ad = solAd.y(4*n+1:5*n,1:end);
time_ad = solAd.x;

y1_DeP = solDeP.y(1:n,1:end);
q1_DeP = solDeP.y(n+1:2*n,1:end);
q2_DeP = solDeP.y(2*n+1:3*n,1:end);
T_DeP = solDeP.y(3*n+1:4*n,1:end);
P_DeP = solDeP.y(4*n+1:5*n,1:end);
time_DeP = solDeP.x;

y1_pur = solPur.y(1:n,1:end); 
q1_pur = solPur.y(n+1:2*n,1:end);
q2_pur = solPur.y(2*n+1:3*n,1:end);
T_pur = solPur.y(3*n+1:4*n,1:end);
P_pur = solPur.y(4*n+1:5*n,1:end);
time_pur = solPur.x;

end
format long


%%% purity and recovery calcs
R = 8.3145;                 % ideal gas constant J/mol/K

%% Flux calculations
% y is an array size n*5 of y1 = 1:n, q1 = n+1:2*n,  
% q2 = 2*n+1:3*n, T = 3*n+1:4*n, P = 4*n+1:5*n


% Pressurisation
    PrVstart = (Prvhalf(1,1:end)+Prvhalf(2,1:end))/2;
    PrScoeff = PrVstart .* P_pr(1,1:end) ./ T_pr(1,1:end)/R;
% Adsorption    
    AdVstart = (Advhalf(1,1:end)+Advhalf(2,1:end))/2;
    AdScoeff = AdVstart .* P_ad(1,1:end) ./ T_ad(1,1:end)/R;
    
    AdVend = (Advhalf(n,1:end)+Advhalf(n+1,1:end))/2;
    AdEcoeff = AdVend .* P_ad(n,1:end) ./ T_ad(n,1:end)/R;
% Depressurisation
    DePVend = (DePvhalf(n,1:end)+DePvhalf(n+1,1:end))/2;
    DePEcoeff = DePVend .* P_DeP(n,1:end) ./ T_DeP(n,1:end)/R;
% Purge 
    PurVstart = (Purvhalf(1,1:end)+Purvhalf(2,1:end))/2;
    PurScoeff = PurVstart .* P_pur(1,1:end) ./ T_pur(1,1:end)/R;
    
    PurgVend = (Purvhalf(n,1:end)+Purvhalf(n+1,1:end))/2;
    PurgEcoeff = PurgVend .* P_pur(n,1:end) ./ T_pur(n,1:end)/R;

    
tf = (cycleend-1)*sum(steptimes) + sum(steptimes(1:2));   

%% Moles out
% should be 16
epl2 = 0.37; % do we want this???
Rbi = 0.5641895835;
A = pi*Rbi^2; % per m^2 atm
rhop = 850;               % Particle density kg/m3
rhob = (1-epl2)*rhop;

% Pressurisation
Pr_CH4_i = epl2 * A .* trapz(time_pr, PrScoeff .* y1_pr(1,1:end));
Pr_H2_i = epl2 * A .* trapz(time_pr, PrScoeff .* MoleyPr2(1,1:end));

Pr_IN_molvol = Pr_H2_i+Pr_CH4_i;

Pr_IN_H2 = Pr_H2_i./(Pr_H2_i+Pr_CH4_i);
Pr_IN_CH4 =Pr_CH4_i./(Pr_H2_i+Pr_CH4_i);
% Adsorption  
Ad_CH4_i = epl2 * A .* trapz(time_ad, AdScoeff .* y1_ad(1,1:end));
Ad_H2_i = epl2 * A .* trapz(time_ad, AdScoeff .* MoleyAd2(1,1:end));
Ad_CH4_o = epl2 * A .* trapz(time_ad, AdEcoeff .* y1_ad(n,1:end));
Ad_H2_o = epl2 * A .* trapz(time_ad, AdEcoeff .* MoleyAd2(n,1:end));

Ad_IN_molvol = Ad_CH4_i+Ad_H2_i;
Ad_OUT_molvol = Ad_CH4_o + Ad_H2_o;

Ad_IN_H2 = Ad_H2_i./(Ad_CH4_i+Ad_H2_i);
Ad_OUT_H2 = Ad_H2_o ./ (Ad_CH4_o + Ad_H2_o);
Ad_IN_CH4 = Ad_CH4_i./(Ad_CH4_i+Ad_H2_i);
Ad_OUT_CH4 = Ad_CH4_o ./ (Ad_CH4_o + Ad_H2_o);
% Depressurisation
DeP_CH4_o = epl2 * A .* trapz(time_DeP, DePEcoeff .* y1_DeP(n,1:end));
DeP_H2_o = epl2 * A .* trapz(time_DeP, DePEcoeff .* MoleyDeP2(n,1:end));

DeP_OUT_molvol = DeP_H2_o + DeP_CH4_o;

DeP_OUT_H2 = DeP_H2_o./(DeP_H2_o+ DeP_CH4_o);
DeP_OUT_CH4 = DeP_CH4_o./(DeP_H2_o+ DeP_CH4_o);
% Purge 
Pur_CH4_i = epl2 * A .* trapz(time_pur, PurScoeff .* y1_pur(1,1:end));
Pur_H2_i = epl2 * A .* trapz(time_pur, PurScoeff .* MoleyPur2(1,1:end));
Pur_CH4_o = epl2 * A .* trapz(time_pur, PurgEcoeff .* y1_pur(n,1:end));
Pur_H2_o = epl2 * A .* trapz(time_pur, PurgEcoeff .* MoleyPur2(n,1:end));

PUR_IN_molvol = Pur_H2_i+Pur_CH4_i;
PUR_OUT_molvol = Pur_H2_o+Pur_CH4_o;

PUR_IN_H2 = Pur_H2_i./(Pur_H2_i+Pur_CH4_i);
PUR_IN_CH4 = Pur_CH4_i./(Pur_H2_i+Pur_CH4_i);
PUR_OUT_H2 = Pur_H2_o./(Pur_H2_o+Pur_CH4_o);
PUR_OUT_CH4 = Pur_CH4_o./(Pur_H2_o+Pur_CH4_o);

%% mass balance
% accum at start of cycle
Acc_g_CH4 = epl2 * A .* (y1_pr(1:end,1) .* P_pr(1:end,1)./T_pr(1:end,1)./R);
gas_start_CH4 = trapz(z,Acc_g_CH4);
Acc_g_H2 = epl2 * A .* (MoleyPr2(1:end,1) .* P_pr(1:end,1)./T_pr(1:end,1)./R);
gas_start_H2 = trapz(z,Acc_g_H2);

% Adsorbate
Acc_a_CH4 = rhob * q1_pr(1:end,1);
absorbed_start_CH4 = trapz(z,Acc_a_CH4);
Acc_a_H2 = rhob * q2_pr(1:end,1);
absorbed_start_H2 = trapz(z,Acc_a_H2);

Mole_acc_start_CH4 = gas_start_CH4 +  absorbed_start_CH4;
Mole_acc_start_H2 = gas_start_H2 +  absorbed_start_H2;

% accum at end of cycle
Acc_ge_CH4 = epl2 * A .* (y1_pur(1:end,end) .* P_pur(1:end,end)./T_pur(1:end,end)./R);
gas_end_CH4 = trapz(z,Acc_ge_CH4);
Acc_ge_H2 = epl2 * A .* (MoleyPur2(1:end,end) .* P_pur(1:end,end)./T_pur(1:end,end)./R);
gas_end_H2 = trapz(z,Acc_ge_H2);
% Adsorbate
Acc_ae_CH4 = rhob * q1_pur(1:end,end);
Acc_ae_H2 = rhob * q2_pur(1:end,end);
absorbed_end_CH4 = trapz(z,Acc_ae_CH4);
absorbed_end_H2 = trapz(z,Acc_ae_H2);

Mole_acc_end_CH4 = gas_end_CH4 +  absorbed_end_CH4;
Mole_acc_end_H2 = gas_end_H2 +  absorbed_end_H2;


% MB = 100*(Ad_OUT_molvol + DeP_OUT_molvol + PUR_OUT_molvol + Mole_acc_end_CH4 + Mole_acc_end_H2) ./ (Pr_IN_molvol + Ad_IN_molvol + PUR_IN_molvol + Mole_acc_start_CH4 + Mole_acc_start_H2);
MB = 100*(Ad_OUT_molvol + DeP_OUT_molvol + PUR_OUT_molvol) ./ (Pr_IN_molvol + Ad_IN_molvol + PUR_IN_molvol);
%% Purity and recovery
Purity_CH4 = (DeP_OUT_CH4 * DeP_OUT_molvol + PUR_OUT_CH4 * PUR_OUT_molvol)/(DeP_OUT_molvol + PUR_OUT_molvol);
Purity_H2 = Ad_OUT_H2;

cons = Purity_H2 - 0.9985            ;
if cons < 0
    con = abs(cons) ;
end


% Recovery_CH4 = (DeP_OUT_CH4 * DeP_OUT_molvol + PUR_OUT_CH4 * PUR_OUT_molvol) / (Ad_IN_CH4*Ad_IN_molvol + Pr_IN_CH4*Pr_IN_molvol + PUR_IN_CH4*PUR_IN_molvol);
Recovery_H2 = (Ad_OUT_H2*Ad_OUT_molvol) / (Ad_IN_H2*Ad_IN_molvol + Pr_IN_H2*Pr_IN_molvol + PUR_IN_H2*PUR_IN_molvol);
%% Work equations
gamma = 1.4;
eff = 0.8;
% Rbi above

W_int_Pr = PrVstart.*P_pr(1,1:end).*((P_pr(1,1:end)./TPF(4)).^(gamma/(gamma-1))-1);
W_int_Ad = AdVstart.*P_ad(1,1:end).*((P_ad(1,1:end)./TPF(4)).^(gamma/(gamma-1))-1);

Work_Pr = 1/eff*pi*Rbi^2* (gamma/(gamma-1)) * trapz(time_pr,W_int_Pr) / 3.6e6; % kWh required
Work_Ad = 1/eff*pi*Rbi^2* (gamma/(gamma-1)) * trapz(time_ad,W_int_Ad) / 3.6e6; % kWh required
Work = Work_Pr + Work_Ad;
Work_per_H2mol = (Work_Pr + Work_Ad)/(Ad_OUT_H2*Ad_OUT_molvol);
%% Productivity
Productivity = (Ad_OUT_H2*Ad_OUT_molvol) / ((Ad_IN_molvol + Pr_IN_H2 + PUR_IN_H2) * sum(steptimes))

objectives = [-Purity_H2, -Recovery_H2 ];
return
for plotingstuff = 1
%%%%% Press plot
    figure(1)
    title('Pressurization')
    subplot(2,3,1)
    mesh(solPr.x,z,solPr.y(1:n,1:end))
    xlabel('t')
    ylabel('bed length')
    zlabel('mole fraction y1')
    title('mole fraction of methane')
    
    subplot(2,3,2)
    mesh(solPr.x,z,MoleyPr2)
    xlabel('t')
    ylabel('bed length')
    zlabel('mole fraction y2')
    title('mole fraction of hydrogen')
    
    subplot(2,3,4)
    mesh(solPr.x,z,solPr.y(n+1:2*n,1:end))
    xlabel('t')
    ylabel('bed length')
    zlabel('loading q1 mol/kg')
    title('loading of methane')
    
    subplot(2,3,5)
    mesh(solPr.x,z,solPr.y(2*n+1:3*n,1:end))
    xlabel('t')
    ylabel('bed length')
    zlabel('loading q2 mol/kg')
    title('loading of hydrogen')
    
    subplot(2,3,3)
    mesh(solPr.x,z,solPr.y(3*n+1:4*n,1:end))
    xlabel('t')
    ylabel('bed length')
    zlabel('Temp')
    title('temp of system')
    
    subplot(2,3,6)
    mesh(solPr.x,z,solPr.y(4*n+1:5*n,1:end))
    xlabel('time')
    ylabel('bed length')
    zlabel('Pressure')
    title('Pressure of system')
 %%%%% Adsorp plot
    figure(2)
    title('Adsorption')
    subplot(2,3,1)
    mesh(solAd.x,z,solAd.y(1:n,1:end))
    xlabel('t')
    ylabel('bed length')
    zlabel('mole fraction y1')
    title('mole fraction of methane')
    
    subplot(2,3,2)
    mesh(solAd.x,z,Moley2)
    xlabel('t')
    ylabel('bed length')
    zlabel('mole fraction y2')
    title('mole fraction of hydrogen')
    
    subplot(2,3,4)
    mesh(solAd.x,z,solAd.y(n+1:2*n,1:end))
    xlabel('t')
    ylabel('bed length')
    zlabel('loading q1 mol/kg')
    title('loading of methane')
    
    subplot(2,3,5)
    mesh(solAd.x,z,solAd.y(2*n+1:3*n,1:end))
    xlabel('t')
    ylabel('bed length')
    zlabel('loading q2 mol/kg')
    title('loading of hydrogen')
    
    subplot(2,3,3)
    mesh(solAd.x,z,solAd.y(3*n+1:4*n,1:end))
    xlabel('t')
    ylabel('bed length')
    zlabel('Temp')
    title('temp of system')
    
    subplot(2,3,6)
    mesh(solAd.x,z,solAd.y(4*n+1:5*n,1:end))
    xlabel('time')
    ylabel('bed length')
    zlabel('Pressure')
    title('Pressure of system')
%%%%% Depress plot   
    figure(3)
    title('DePressurization')
    subplot(2,3,1)
    mesh(solDeP.x,z,flip(solDeP.y(1:n,1:end)))
    xlabel('t')
    ylabel('bed length')
    zlabel('mole fraction y1')
    title('mole fraction of methane')
    
    subplot(2,3,2)
    mesh(solDeP.x,z,flip(MoleyDeP2))
    xlabel('t')
    ylabel('bed length')
    zlabel('mole fraction y2')
    title('mole fraction of hydrogen')
    
    subplot(2,3,4)
    mesh(solDeP.x,z,flip(solDeP.y(n+1:2*n,1:end)))
    xlabel('t')
    ylabel('bed length')
    zlabel('loading q1 mol/kg')
    title('loading of methane')
    
    subplot(2,3,5)
    mesh(solDeP.x,z,flip(solDeP.y(2*n+1:3*n,1:end)))
    xlabel('t')
    ylabel('bed length')
    zlabel('loading q2 mol/kg')
    title('loading of hydrogen')
    
    subplot(2,3,3)
    mesh(solDeP.x,z,flip(solDeP.y(3*n+1:4*n,1:end)))
    xlabel('t')
    ylabel('bed length')
    zlabel('Temp')
    title('temp of system')
    
    subplot(2,3,6)
    mesh(solDeP.x,z,flip(solDeP.y(4*n+1:5*n,1:end)))
    xlabel('time')
    ylabel('bed length')
    zlabel('Pressure')
    title('Pressure of system')
%%%%% Purge plot 
    figure(4)
    title('Purge')
    subplot(2,3,1)
    mesh(solPur.x,z,flip(solPur.y(1:n,1:end)))
    xlabel('t')
    ylabel('bed length')
    zlabel('mole fraction y1')
    title('mole fraction of methane')
    
    subplot(2,3,2)
    mesh(solPur.x,z,flip(MoleyPur2))
    xlabel('t')
    ylabel('bed length')
    zlabel('mole fraction y2')
    title('mole fraction of hydrogen')
    
    subplot(2,3,4)
    mesh(solPur.x,z,flip(solPur.y(n+1:2*n,1:end)))
    xlabel('t')
    ylabel('bed length')
    zlabel('loading q1 mol/kg')
    title('loading of methane')
    
    subplot(2,3,5)
    mesh(solPur.x,z,flip(solPur.y(2*n+1:3*n,1:end)))
    xlabel('t')
    ylabel('bed length')
    zlabel('loading q2 mol/kg')
    title('loading of hydrogen')
    
    subplot(2,3,3)
    mesh(solPur.x,z,flip(solPur.y(3*n+1:4*n,1:end)))
    xlabel('t')
    ylabel('bed length')
    zlabel('Temp')
    title('temp of system')
    
    subplot(2,3,6)
    mesh(solPur.x,z,flip(solPur.y(4*n+1:5*n,1:end)))
    xlabel('time')
    ylabel('bed length')
    zlabel('Pressure')
    title('Pressure of system')

%%%%% Velocity plots 
    figure(6)
    title('Velocity plot')
    subplot(2,2,1)
    mesh(solPr.x,zhalf,Prvhalf)
    xlabel('t')
    ylabel('bed length')
    zlabel('Velocity')
    title('Pressurization')
    subplot(2,2,2)
    mesh(solAd.x,zhalf,Advhalf)
    xlabel('t')
    ylabel('bed length')
    zlabel('Velocity')
    title('Adsorption')
    subplot(2,2,3)
    mesh(solDeP.x,zhalf,flip(-1.*DePvhalf))
    xlabel('t')
    ylabel('bed length')
    zlabel('Velocity')
    title('Depressurization')
    subplot(2,2,4)
    mesh(solPur.x,zhalf,flip(-1.*Purvhalf))
    xlabel('t')
    ylabel('bed length')
    zlabel('Velocity')
    title('Purge')

end  
end
%% Solver function
function [Prout, out, DePout, Purout, timeOut, YAdOut, Phigh, Plow] = adsorptionSolver(~,z,yFeed,TPFeed,steptimes,stepend,cycleend)
n = numel(z);      % Size of mesh grid
h = diff(z);
h = h(1);
% start up condidtions
YAdOut = 0;
timeOut = 0;
highP = TPFeed(3);
lowP = TPFeed(4);
Tw = 300;                   % Ambient Tempurature K
Pw = 1e5;                % Ambient Pressure Pa
y1w = 0.80;
cycleEnd = cycleend;
stepEnd = stepend;
signumber = 5;
Quitcycle = 0;

% initialse checkers
PrevPr = zeros(5*n,1);
Pr = zeros(5*n,1);
PrevAd = zeros(5*n,1);
Ad = zeros(5*n,1);
PrevDeP = zeros(5*n,1);
DeP = zeros(5*n,1);
PrevPur = zeros(5*n,1);
Pur = zeros(5*n,1);

for cycleNo = 1:cycleEnd
    
    for stepNo = 1:stepEnd
        t0 = 0;
        %tf = 60;
        tf = steptimes(stepNo);
        tspan = [t0 tf];
        % Start up conditions Conditions
        if cycleNo ==1 && stepNo == 1
            y1_0       = zeros(n,1) + y1w;     

            q1_0        = zeros(n,1);
            q2_0        = zeros(n,1);

            T_0         = zeros(n,1) + Tw;

            P_0         = zeros(n,1) + Pw;
        % Previous conditions Conditions
        else
            y1_0       = zeros(n,1) + Yinit;   

            q1_0        = zeros(n,1) + q1init;
            q2_0        = zeros(n,1) + q2init;

            T_0         = zeros(n,1) + Tinit;

            P_0         = zeros(n,1) + Pinit;
        end

        % y is an array size n*4 of c1 = 1:n, q1 = n+1:2*n, c2 = 2*n+1:3*n, q2 = 3*n+1:4*n
        y0      = [y1_0; q1_0; q2_0; T_0; P_0];
        
        
        %out = ode15s(@(t,y) adsorptionModel(t,y,h,n,yFeed,TPFeed),tspan,y0);
        adparams = ode15s(@(t,y) adsorptionModel(t,y,h,n,yFeed,TPFeed,stepNo,YAdOut,timeOut,highP,lowP),tspan,y0);

        % Step bed conditions
        Yinit = adparams.y(1:n,end);
        q1init = adparams.y(n+1:2*n,end);
        q2init = adparams.y(2*n+1:3*n,end);
        Tinit = adparams.y(3*n+1:4*n,end);
        Pinit = adparams.y(4*n+1:5*n,end);

 

        % Bed end (boundary) conditions
        % Showing condensed values at end of times
        if stepNo == 1
        % for checking css
        Pr(1:n) = round(Yinit,signumber,"significant");
        Pr(n+1:2*n) = round(q1init,signumber,"significant");
        Pr(2*n+1:3*n) = round(q2init,signumber,"significant");
        Pr(3*n+1:4*n) = round(Tinit,signumber,"significant");
        Pr(4*n+1:5*n) = round(Pinit/1e5,signumber,"significant");
        end 
        % Reverse steps bed conditions
        if stepNo == 3

        DeP(1:n) = round(Yinit,signumber,"significant");
        DeP(n+1:2*n) = round(q1init,signumber,"significant");
        DeP(2*n+1:3*n) = round(q2init,signumber,"significant");
        DeP(3*n+1:4*n) = round(Tinit,signumber,"significant");
        DeP(4*n+1:5*n) = round(Pinit/1e5,signumber,"significant");
        end
        
         % Reverse steps bed conditions
        if stepNo == 2 || stepNo == 4
        Yinit = flip(Yinit);
        q1init = flip(q1init);
        q2init = flip(q2init);
        Tinit = flip(Tinit);
        Pinit = flip(Pinit);
        end
        
        

        
        if stepNo == 4
        % set press of boundary conditions
        lowP = Pinit(1);
        
        Pur(1:n) = round(Yinit,signumber,"significant");
        Pur(n+1:2*n) = round(q1init,signumber,"significant");
        Pur(2*n+1:3*n) = round(q2init,signumber,"significant");
        Pur(3*n+1:4*n) = round(Tinit,signumber,"significant");
        Pur(4*n+1:5*n) = round(Pinit/1e5,signumber,"significant");
        end
        
        if stepNo == 2
        % set press of boundary conditions
        Ad(1:n) = round(Yinit,signumber,"significant");
        Ad(n+1:2*n) = round(q1init,signumber,"significant");
        Ad(2*n+1:3*n) = round(q2init,signumber,"significant");
        Ad(3*n+1:4*n) = round(Tinit,signumber,"significant");
        Ad(4*n+1:5*n) = round(Pinit/1e5,signumber,"significant");
        
        highP = Pinit(n);
        
        YAdOut = adparams.y(n,1:end);
        timeOut = adparams.x;
            % purity checker
            timecheck = (cycleNo-1)*sum(steptimes)+ sum(steptimes(1:2));
            Moleyy2 = 1-adparams.y(1:n,1:end);
            format long
            purityh = 100*sum(Moleyy2(n,1:end)) / (sum(adparams.y(n,1:end)) + sum(Moleyy2(n,1:end)));
            fprintf('Purity of hydrogen is %f%% after AD. Time: %4.2fs\n', purityh, timecheck)
        end
       
        if (cycleNo == cycleEnd || Quitcycle == 1) && stepNo == 1
        Prout   = adparams;
        Plow = lowP;
        elseif (cycleNo == cycleEnd || Quitcycle == 1) && stepNo == 2
        out = adparams;
        elseif (cycleNo == cycleEnd || Quitcycle == 1) && stepNo == 3
        DePout   = adparams;
        Phigh = highP;
        elseif (cycleNo == cycleEnd || Quitcycle == 1) && stepNo == 4 
        Purout  = adparams;
        end
    
    end
    % check is prev = this one
    if Quitcycle == 1
        break
    end
        
        statesPC = [PrevPr;PrevAd;PrevDeP;PrevPur];
        statesFC = [Pr;Ad;DeP;Pur];
        
        CSS_perdiff = abs((statesPC-statesFC)./statesFC);
        Approx_CSS = mean(CSS_perdiff);
        disp(Approx_CSS)
         all(CSS_perdiff <= 2.5e-3)

        if all(CSS_perdiff <= 2.5e-3) 
            Quitcycle = 1;
        elseif cycleNo ~= cycleEnd
            PrevPr = Pr;
            PrevAd = Ad;
            PrevDeP = DeP;
            PrevPur = Pur;
        end
    % if not then save prev to current
    % if yes set cycle end to next and pull out of loop
end

end
%% Adsorption model
function [dydt, Vhalf, Pspan, Phalfspan] = adsorptionModel(t,y,h,n,yiFeed,TPFeed,stepNo,yAdOut,TimeoutAD,HighP,LowP)
%stepNo=2;
% Variables allocation
% y is an array size n*5 of y1 = 1:n, q1 = n+1:2*n,  
% q2 = 2*n+1:3*n, T = 3*n+1:4*n, P = 4*n+1:5*n
y1 = max(y(1:n),0);
y1 = min(y1,1);
y2 = 1 - y1;
y2 = max(y2,0);
y2 = min(y2,1);

q1 = max(y(n+1:2*n),0);
q2 = max(y(2*n+1:3*n),0);


T  = y(3*n+1:4*n);
P  = y(4*n+1:5*n);


% Z half points
yhalf = zeros(n+1,1);
Thalf = zeros(n+1,1);
Phalf = zeros(n+1,1);
uhalf = zeros(n+1,1);
deltPhalf = zeros(n+1,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constant physical properties and basic simulation parameters
R = 8.3145;                 % ideal gas constant J/mol/K 
D   = 1.3e-5;               % Axial Dispersion coefficient m2/s
epl = 0.37;               % Void fraction of bed
eplp = 0.546;               % Void fraction of particle
eplt = epl + (1-epl)*eplp;  % Total void fraction 
k   = [0.136;0.259];        % (lumped) Mass Transfer Coefficient 1/s
% k   = [0.195;0.7];        % (lumped) Mass Transfer Coefficient 1/s
% rhop = 716.3;               % Particle density kg/m3
rhop = 850;               % Particle density kg/m3
rhob = rhop*(1-epl);
Cps = 1046.7;               % heat capacity of solid J/kg/K 
Kl = 0.09;                % Thermal diffusivity J/m/s/K
% Kl = 1.2e-5;                % Thermal diffusivity J/m/s/K
deltaH = [24124; 8420];     % heat of adsorption J/mol
ufeed = TPFeed(1);
Tfeed = TPFeed(2);
PH = TPFeed(3); 
PL = TPFeed(4);
ufeed_purge = TPFeed(5);
lambdaP = 1;               % rate of pressure change (1/s)
lambdaD = 2;               % rate of pressure change (1/s)


% Ergun equation parameters
visc = 1.1e-5;     % gas viscosity kg/m/s
Rp = 1.15e-3;     % particle radius m

Acoeff = 1.75*(1-epl)/(2*Rp*epl^3);
Bcoeff = 150*(1-epl)^2/(4*Rp^2*epl^3);

% Constants for CH4
c1a = -0.703029; c1b = 108.4773; c1c = -42.52157; c1d = 5.862788; c1e = 0.678565;
% Constants for H2
c2a = 33.066178; c2b = -11.363417; c2c = 11.432816; c2d = -2.772874; c2e = -0.158558;

% Cpg = heat capacity (J/mol*K) 
Cpg1 = c1a + c1b*(T/1e3) + c1c*(T/1e3).^2 + c1d*(T/1e3).^3 + c1e./((T/1e3).^2);
Cpg2 = c2a + c2b*(T/1e3) + c2c*(T/1e3).^2 + c2d*(T/1e3).^3 + c2e./((T/1e3).^2);
Cpg = (y1.*Cpg1 + y2.*Cpg2);
% adsorbed phase heat capacity
Cpa1 = q1./(q1+q2 + 10e-10);
Cpa2 = q2./(q1+q2 + 10e-10);
Cpa = (Cpa1.*Cpg1 + Cpa2.*Cpg2);
% gas density
rhog = (P .* (y1*16.04e-3 + y2*2.02e-3)) ./ (R*T);  % ρ = PM/RT, M: H2 = 2.02 g/mol CH4 = 16.04 g/mol Air= 28.97 g/mol 50/50 H2 and CH4 = 9.03

%%%%%%%% Boundary conditions: j = 0.5 and n + 0.5 %%%%%%%%%%%%%%%%%%%%%%%%    
    %% Pressurization
if stepNo == 1
   
    Phalf(1) = PH - (PH - LowP)*exp(-lambdaP*t);
 
    Phalf(n+1) = P(n);
    
    yhalf(1) = (y1(1) + yiFeed*uhalf(1)*h/D/2) / (1 + uhalf(1)*h/D/2);
    yhalf(n+1) = y1(n);
    
    rhogInlet = ((Phalf(1)) .* (yhalf(1)*16.04e-3 + (1-yhalf(1))*2.02e-3)) ./ (R*Tfeed);
    Cpg1in = c1a + c1b*(Tfeed/1e3) + c1c*(Tfeed/1e3).^2 + c1d*(Tfeed/1e3).^3 + c1e./((Tfeed/1e3).^2);
    Cpg2in = c2a + c2b*(Tfeed/1e3) + c2c*(Tfeed/1e3).^2 + c2d*(Tfeed/1e3).^3 + c2e./((Tfeed/1e3).^2);
    CpgInlet = (yhalf(1).*Cpg1in + (1-yhalf(1)).*Cpg2in);
    
    Thalf(n+1) = T(n);
    Thalf(1) = (T(1) + Tfeed*uhalf(1)*epl*rhogInlet*CpgInlet*h/Kl/2) / (1 + uhalf(1)*epl*rhogInlet*CpgInlet*h/Kl/2);


    deltPhalf(1) = (P(1)-Phalf(1));
    if deltPhalf(1) < 0
        uhalf(1) = (-(Bcoeff.*visc) + sqrt((Bcoeff.*visc)^2 - 4 * 2/h*deltPhalf(1) .* (Acoeff.*rhogInlet))) ./ (2*(Acoeff.*rhogInlet));
    elseif deltPhalf(1) > 0
        deltPhalf(1) = -deltPhalf(1);
        uhalf(1) = -(-(Bcoeff.*visc) + sqrt((Bcoeff.*visc)^2 - 4 * 2/h*deltPhalf(1) .* (Acoeff.*rhogInlet))) ./ (2*(Acoeff.*rhogInlet));
    else
        uhalf(1) = 0;
    end

    uhalf(n+1) = 0;
        
  
   %% Adsorption
elseif stepNo == 2
    
     uhalf(1) = ufeed - ufeed*exp(-0.5*t);
    Phalf(1) = P(1) + h/2*(Bcoeff*uhalf(1)*visc + Acoeff.*rhog(1)*uhalf(1)^2);
    Phalf(n+1) = PH;
    
    yhalf(1) = (y1(1) + yiFeed*uhalf(1)*h/D/2) / (1 + uhalf(1)*h/D/2);
    yhalf(n+1) = y1(n);
    
    rhogInlet = ((Phalf(1)) .* (yhalf(1)*16.04e-3 + (1-yhalf(1))*2.02e-3)) ./ (R*Tfeed);
    Cpg1in = c1a + c1b*(Tfeed/1e3) + c1c*(Tfeed/1e3).^2 + c1d*(Tfeed/1e3).^3 + c1e./((Tfeed/1e3).^2);
    Cpg2in = c2a + c2b*(Tfeed/1e3) + c2c*(Tfeed/1e3).^2 + c2d*(Tfeed/1e3).^3 + c2e./((Tfeed/1e3).^2);
    CpgInlet = (yhalf(1).*Cpg1in + (1-yhalf(1)).*Cpg2in);
    
    Thalf(1) = (T(1) + Tfeed*uhalf(1)*epl*rhogInlet*CpgInlet*h/Kl/2) / (1 + uhalf(1)*epl*rhogInlet*CpgInlet*h/Kl/2);
    Thalf(n+1) = T(n);
    
    
    rhogOutlet = ((Phalf(n+1)) .* (yhalf(n+1)*16.04e-3 + (1-yhalf(n+1))*2.02e-3)) ./ (R*Thalf(n+1));
    
    deltPhalf(n+1) = (Phalf(n+1)-P(n));
    if deltPhalf(n+1) < 0
        uhalf(n+1) = (-(Bcoeff.*visc) + sqrt((Bcoeff.*visc)^2 - 4 * 2/h*deltPhalf(n+1) .* (Acoeff.*rhogOutlet))) ./ (2*(Acoeff.*rhogOutlet));
    elseif deltPhalf(n+1) > 0
        deltPhalf(n+1) = -deltPhalf(n+1);
        uhalf(n+1) = -(-(Bcoeff.*visc) + sqrt((Bcoeff.*visc)^2 - 4 * 2/h*deltPhalf(n+1) .* (Acoeff.*rhogOutlet))) ./ (2*(Acoeff.*rhogOutlet));
    else
        uhalf(n+1) = 0;
    end
    
%% Forward Depressurization
% can do reverse depress by switching i ititial conditions 
elseif stepNo == 3

    Phalf(n+1) = PL + (HighP - PL)*exp(-lambdaD*t);
    
    Phalf(1) = P(1);
    
    yhalf(1) = y1(1);
    yhalf(n+1) = y1(n);
    
    Thalf(1) = T(1);
    Thalf(n+1) = T(n);
    
    rhogOutlet = ((Phalf(n+1)) .* (yhalf(n+1)*16.04e-3 + (1-yhalf(n+1))*2.02e-3)) ./ (R*Thalf(n+1));
    
    deltPhalf(n+1) = (Phalf(n+1)-P(n));
    if deltPhalf(n+1) < 0
        uhalf(n+1) = (-(Bcoeff.*visc) + sqrt((Bcoeff.*visc)^2 - 4 * 2/h*deltPhalf(n+1) .* (Acoeff.*rhogOutlet))) ./ (2*(Acoeff.*rhogOutlet));
    elseif deltPhalf(n+1) > 0
        deltPhalf(n+1) = -deltPhalf(n+1);
        uhalf(n+1) = -(-(Bcoeff.*visc) + sqrt((Bcoeff.*visc)^2 - 4 * 2/h*deltPhalf(n+1) .* (Acoeff.*rhogOutlet))) ./ (2*(Acoeff.*rhogOutlet));
    else
        uhalf(n+1) = 0;
    end
    
    uhalf(1) = 0;
  %% Forward Purge
elseif stepNo == 4
    
    uhalf(1) = ufeed_purge - ufeed_purge*exp(-0.5*t);
    Phalf(1) = P(1) + h/2*(Bcoeff*uhalf(1)*visc + Acoeff.*rhog(1)*uhalf(1)^2);
    Phalf(n+1) = PL;
    yFeedfromAD2 = sum(trapz(TimeoutAD,yAdOut))/TimeoutAD(end);
    yhalf(1) = (y1(1) + yFeedfromAD2*uhalf(1)*h/D/2) / (1 + uhalf(1)*h/D/2);
    yhalf(n+1) = y1(n);

  
    rhogInlet = ((Phalf(1)) .* (yhalf(1)*16.04e-3 + (1-yhalf(1))*2.02e-3)) ./ (R*Tfeed);
    Cpg1in = c1a + c1b*(Tfeed/1e3) + c1c*(Tfeed/1e3).^2 + c1d*(Tfeed/1e3).^3 + c1e./((Tfeed/1e3).^2);
    Cpg2in = c2a + c2b*(Tfeed/1e3) + c2c*(Tfeed/1e3).^2 + c2d*(Tfeed/1e3).^3 + c2e./((Tfeed/1e3).^2);
    CpgInlet = (yhalf(1).*Cpg1in + (1-yhalf(1)).*Cpg2in);
    
    Thalf(1) = (T(1) + Tfeed*uhalf(1)*epl*rhogInlet*CpgInlet*h/Kl/2) / (1 + uhalf(1)*epl*rhogInlet*CpgInlet*h/Kl/2);
    Thalf(n+1) = T(n);
    
    rhogOutlet = ((Phalf(n+1)) .* (yhalf(n+1)*16.04e-3 + (1-yhalf(n+1))*2.02e-3)) ./ (R*Thalf(n+1));

    deltPhalf(n+1) = (Phalf(n+1)-P(n));
    if deltPhalf(n+1) < 0
        uhalf(n+1) = (-(Bcoeff.*visc) + sqrt((Bcoeff.*visc)^2 - 4 * 2/h*deltPhalf(n+1) .* (Acoeff.*rhogOutlet))) ./ (2*(Acoeff.*rhogOutlet));
    elseif deltPhalf(n+1) > 0
        deltPhalf(n+1) = -deltPhalf(n+1);
        uhalf(n+1) = -(-(Bcoeff.*visc) + sqrt((Bcoeff.*visc)^2 - 4 * 2/h*deltPhalf(n+1) .* (Acoeff.*rhogOutlet))) ./ (2*(Acoeff.*rhogOutlet));
    else
        uhalf(n+1) = 0;
    end
    
end

%% %%%%%% Wall values: j = 1.5 %%%%%%%%%%%%%%%%%%%%%%%%
alpha0P = (2/3) / (P(2)-P(1) + 1e-10)^4;
alpha1P = (1/3) / (2*(P(1)-Phalf(1)) + 1e-10)^4;
Phalf(2) = alpha0P/(alpha0P+alpha1P)*(0.5*(P(1)+P(2))) + alpha1P/(alpha0P+alpha1P)*(2*P(1)-Phalf(1));

alpha0y = (2/3) / (y1(2)-y1(1) + 1e-10)^4;
alpha1y = (1/3) / (2*(y1(1)-yhalf(1)) + 1e-10)^4;
yhalf(2) = alpha0y/(alpha0y+alpha1y)*(0.5*(y1(1)+y1(2))) + alpha1y/(alpha0y+alpha1y)*(2*y1(1)-yhalf(1));

alpha0T = (2/3) / (T(2)-T(1) + 1e-10)^4;
alpha1T = (1/3) / (2*(T(1)-Thalf(1)) + 1e-10)^4;
Thalf(2) = alpha0T/(alpha0T+alpha1T)*(0.5*(T(1)+T(2))) + alpha1T/(alpha0T+alpha1T)*(2*T(1)-Thalf(1));

%%%%%%%% Wall values: j = 2.5 to n - 0.5 %%%%%%%%%%%%%%%%%%%%%%%%
alpha0P = (2/3) ./ (P(3:n)-P(2:n-1) + 1e-10).^4;
alpha1P = (1/3) ./ (P(2:n-1)-P(1:n-2) + 1e-10).^4;
Phalf(3:n) = alpha0P./(alpha0P+alpha1P).*(0.5.*(P(2:n-1)+P(3:n))) + alpha1P./(alpha0P+alpha1P).*(1.5.*P(2:n-1)- 0.5* P(1:n-2));

alpha0y = (2/3) ./ (y1(3:n)-y1(2:n-1) + 1e-10).^4;
alpha1y = (1/3) ./ (y1(2:n-1)-y1(1:n-2) + 1e-10).^4;
yhalf(3:n) = alpha0y./(alpha0y+alpha1y).*(0.5*(y1(2:n-1)+y1(3:n))) + alpha1y./(alpha0y+alpha1y).*(1.5.*y1(2:n-1)- 0.5* y1(1:n-2));

alpha0T = (2/3) ./ (T(3:n)-T(2:n-1) + 1e-10).^4;
alpha1T = (1/3) ./ (T(2:n-1)-T(1:n-2)+ 1e-10).^4;
Thalf(3:n) = alpha0T./(alpha0T+alpha1T).*(0.5*(T(2:n-1)+T(3:n))) + alpha1T./(alpha0T+alpha1T).*(1.5.*T(2:n-1)- 0.5* T(1:n-2));



%% %%%%%% Velocity wall values: j 1.5 to n - 0.5 %%%%%%%%%%%%%%%%%%%%%%%%
rhoghalf = (Phalf .* (yhalf*16.04e-3 + (1-yhalf)*2.02e-3)) ./ (R*Thalf); % ρ = PM/RT, M: H2 = 2.02 g/mol CH4 = 16.04 g/mol Air= 28.97 g/mol 50/50 H2 and CH4 = 9.03
deltPhalf(2:n) = (P(2:n) - P(1:n-1));
%deltPhalf(2:n) = (-Phalf(3:n+1)+8*P(2:n)-8*P(1:n-1)+Phalf(1:n-1))/6;

% Velocity calcs: uhalf(2:n) needs P(2:n) - P(1:n-1))
for i = 2:n
    if deltPhalf(i) < 0
        uhalf(i) = (-(Bcoeff.*visc) + sqrt((Bcoeff.*visc)^2 - 4 * deltPhalf(i)/h .* (Acoeff.*rhoghalf(i)))) ./ (2*(Acoeff.*rhoghalf(i)));
    elseif deltPhalf(i) > 0
        deltPhalf(i) = -deltPhalf(i);
        uhalf(i) = -(-(Bcoeff.*visc) + sqrt((Bcoeff.*visc)^2 - 4 * deltPhalf(i)/h .* (Acoeff.*rhoghalf(i)))) ./ (2*(Acoeff.*rhoghalf(i)));
    else
        uhalf(i) = 0;
    end
end

u = zeros(n,1);
deltP = (Phalf(2:n+1) - Phalf(1:n));
for i = 1:n
    if deltP(i) < 0
        u(i) = (-(Bcoeff.*visc) + sqrt((Bcoeff.*visc)^2 - 4 * deltP(i)/h .* (Acoeff.*rhog(i)))) ./ (2*(Acoeff.*rhog(i)));
    elseif deltP(i) > 0
        deltP(i) = -deltP(i);
        u(i) = -(-(Bcoeff.*visc) + sqrt((Bcoeff.*visc)^2 - 4 * deltP(i)/h .* (Acoeff.*rhog(i)))) ./ (2*(Acoeff.*rhog(i)));
    else
        u(i) = 0;
    end
end

%%%%%%%%%%%%%%% Langmuir Equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % CH4
a11 = -9.5592;	a21 = 4638.5;	b11 = 3.725e-4/1e5;	b21 = 1.443e4;
% % H2
a12 = -23.131;	a22 = 8069.1;	b12 = 2.248/1e5;	b22 = -1.435e4;
qs1 = a11 + (a21./T);
qs2 = a12 + (a22./T);
B1 = b11.*exp(b21./(8.3145.*T));
B2 = b12.*exp(b22./(8.3145*T));
qstar1 = (qs1.*B1.*P.*y1) ./ (1+((P.*y1.*B1)+(P.*(1-y1).*B2))); %mols/kg
qstar2 = (qs2.*B2.*P.*(1-y1)) ./ (1+((P.*y1.*B1)+(P.*(1-y1).*B2)));

% %% LRC method
% Plrc = P/101325;
% % CH4
% k1c = 23.86e-3;	k2c = -0.0562e-3;	k3c = 2.81093e-3;	k4c = 1220; k5c = 1.628; k6c = -248.9;
% 
% % % H2
% k1h = 7.34345e-3;	k2h = -0.013e-3;	k3h = 0.932e-3;	k4h = 506.306; k5h = 0.586972; k6h = 154.455;
% 
% qsi1 = k1c + k2c*T;
% qsi2 = k1h + k2h*T;
% 
% B1 = k3c.*exp(k4c./T);
% B2 = k3h.*exp(k4h./T);
% 
% n1 = (k5c+k6c./T);
% n2 = (k5h+k6h./T); 
% 
% qstar1 = (qsi1.*B1.*(Plrc.*y1).^n1) ./ (1+(((Plrc.*y1).^n1.*B1)+((Plrc.*y2).^n2.*B2)))*1e3; %mols/kg
% qstar2 = (qsi2.*B2.*(Plrc.*y2).^n2) ./ (1+(((Plrc.*y1).^n1.*B1)+((Plrc.*y2).^n2.*B2)))*1e3; %mols/kg

LT1 = qstar1-q1;
LT2 = qstar2-q2;

%%%%%%%%%%%%%%% Eqautions for loading of component i %%%%%%%%%%%%%%
dq1dt = k(1)*LT1;
dq2dt = k(2)*LT2;

%%%%%%%%%%%%%% Balance equations

dPdt = zeros(n,1);
dy1dt = zeros(n,1);
dTdt = zeros(n,1);
dTwdt = zeros(n,1);

dTdz1 = zeros(n,1);
dTdz2 = zeros(n,1);
dTdz3 = zeros(n,1);
dTdz4 = zeros(n,1);
dPdz1 = zeros(n,1);
dPdz2 = zeros(n,1);
dPdz3 = zeros(n,1);
dydz1 = zeros(n,1);
dydz2 = zeros(n,1);
dydz3 = zeros(n,1);
dydz4 = zeros(n,1);
dydz5 = zeros(n,1);

% % Wall balance, not used
% Tatm = 292.5 ;                   % ambient temp
% hi = 38.4928;                   % heat transfer coefficient J/m2/s/K
% ho = 14.2256;
% 
% rhow = 7830;
% Cpw = 502.416;
% Rbi = 2.20447e-2;          % Bed inner radius m 
% Rbo = 2.2073e-2;        % Bed outer raduis m
% Aw = pi*(Rbo^2-Rbi^2);
% 
% dTwdt(1:n) =  0;% 1./(rhow*Cpw*Aw) .* ((2*pi*Rbi*hi*(T(1:n) - Tw(1:n)) - 2*pi*Rbo*ho*(T(1:n) - Tatm)));

%% Energy balance
sink =1./(eplt.*rhog.*Cpg + ((1-epl).*Cpa.*(q1+q2) + rhob.*Cps));

dTdz1(1) = Kl.*( (T(2)-T(1))./h - 2.*(T(1)-Thalf(1))./h )./h;
dTdz1(2:n-1) = Kl.*( (T(3:n)-T(2:n-1))./h - (T(2:n-1)-T(1:n-2))./h )./h;
dTdz1(n) = Kl.*( 2.*(Thalf(n+1)-T(n))./h - (T(n)-T(n-1))./h )./h;

dTdz2(1:n) = -Cpg(1:n).*epl./R.*( (Phalf(2:n+1).*uhalf(2:n+1)-Phalf(1:n).*uhalf(1:n)) - T(1:n).*(Phalf(2:n+1).*uhalf(2:n+1)./Thalf(2:n+1)-Phalf(1:n).*uhalf(1:n)./Thalf(1:n)) )/h;

dTdz3(1:n) = rhob.*(deltaH(1).*dq1dt(1:n) + deltaH(2).*dq2dt(1:n));

% Temporal EB
dTdt(1:n) = sink.* (dTdz1(1:n) + dTdz2(1:n) + dTdz3(1:n));% + dTdz4(1:n));

%% total mass balance

dPdz1(1:n) = - T(1:n) .* (uhalf(2:n+1).*Phalf(2:n+1)./Thalf(2:n+1) - uhalf(1:n).*Phalf(1:n)./Thalf(1:n))./h;

% no rhop.*
dPdz2(1:n) = - (R.*rhob.*T(1:n)./epl) .* (dq1dt(1:n)+dq2dt(1:n));

dPdz3(1:n) = P(1:n)./T(1:n).*dTdt(1:n);

% Temporal EB
dPdt(1:n) = dPdz1(1:n) + dPdz2(1:n) + dPdz3(1:n);




%% component mass balance
dy2dz = zeros(n,1);
dy2dz(1) = ( (y1(2)-y1(1))./h - 2.*(y1(1)-yhalf(1))./h )./h;
dy2dz(2:n-1) = ( y1(3:n)- 2.*y1(2:n-1) + y1(1:n-2) )./h./h;
dy2dz(n) = ( 2.*(yhalf(n+1)-y1(n))/h - (y1(n)-y1(n-1))./h )./h;


dydz1(1:n) = D.*( dy2dz(1:n) + ( (Phalf(2:n+1)-Phalf(1:n))./h .* (yhalf(2:n+1) - yhalf(1:n))./h )./P(1:n) - ( (Thalf(2:n+1)-Thalf(1:n))./h .* (yhalf(2:n+1) - yhalf(1:n))./h )./T(1:n));

dydz2(1:n) = -(T(1:n)./P(1:n)).*( (Phalf(2:n+1).*uhalf(2:n+1).*yhalf(2:n+1)./Thalf(2:n+1)-Phalf(1:n).*uhalf(1:n).*yhalf(1:n)./Thalf(1:n)) - y1(1:n).*(Phalf(2:n+1).*uhalf(2:n+1)./Thalf(2:n+1)-Phalf(1:n).*uhalf(1:n)./Thalf(1:n)) )./h;

dydz3(1:n) = - ((rhob*R*T(1:n))/epl./P(1:n)).*(dq1dt(1:n)-y1(1:n).*(dq1dt(1:n)+dq2dt(1:n)));

dydz4(1:n) = - y1(1:n)./P(1:n).*dPdt(1:n);

dydz5(1:n) = y1(1:n)./T(1:n).*dTdt(1:n);

dy1dt(1:n) = dydz1(1:n) + dydz2(1:n) + dydz3(1:n);% + dydz4(1:n) + dydz5(1:n);


%% caclutaions for graphs
Vhalf=uhalf;
Pspan = Phalf;

Phalfspan = zeros(n,1);
Phalfspan(1) = 2*(-Phalf(1)+P(1))/h;
Phalfspan(n+1) = 2*(Phalf(n+1)-P(n))/h;
Phalfspan(2:n) = (P(2:n) - P(1:n-1))./h;

%%%%%%%%%%%%%%%% Concatenate vectors of time derivatives
dydt = [dy1dt;dq1dt;dq2dt;dTdt;dPdt];
end

function condensedArray = condenseArrayToFiveMeans(originalArray)
    m = 10; % new array size 
    n = length(originalArray); % Length of the original array
    segmentSize = floor(n / m); % Size of each segment, floor to handle non-divisible lengths
    condensedArray = zeros(1, m); % Initialize the array for 5 averaged values
    
    for i = 1:m
        if i < m
            % For the first 4 segments
            startIndex = (i - 1) * segmentSize + 1;
            endIndex = i * segmentSize;
        else
            % For the last segment, include any remaining elements
            startIndex = (i - 1) * segmentSize + 1;
            endIndex = n; % Go to the end of the array
        end
        % Calculate the mean of the segment
        condensedArray(i) = mean(originalArray(startIndex:endIndex));
    end
end
