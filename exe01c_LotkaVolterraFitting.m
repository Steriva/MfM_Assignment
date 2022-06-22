clear; clc; close all

addpath('matlab')

colours = struct2cell(colori());
markersize = 10;
%% Load data
opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";
opts.VariableNames = ["Year", "SnowshoeHarePopulationthousands", "CanadaLynxPeltsthousands"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
populationData = readtable("data/populationData.txt", opts);
populationData = table2array(populationData);
clear opts

initialTime = populationData(1,1);

% Resize the time
time = (populationData(:,1)-initialTime);

IC = [populationData(1,2); populationData(1,3)];

%% Minimization Problem

% Set nondefault solver options
options2 = optimoptions("fmincon","Algorithm","interior-point","Display","iter",...
    "PlotFcn","optimplotfvalconstr",'MaxFunctionEvaluations',1000);

% options = optimoptions('lsqnonlin','OptimalityTolerance',1e-12,'Algorithm','Levenberg-Marquardt',...
%    'FunctionTolerance',1e-12,'StepTolerance',1e-12,'UseParallel',false,'Display','iter');

objfun = @(mu) computingOF(mu, populationData(:,2:3), time, IC);

mu0 = [0.7; 0.1; 0.1; 0.7];

% [muSolution,~] = lsqnonlin(objfun,mu0,[0 0 0 0], [+Inf +Inf +Inf +Inf],options);

[muSolution,~] = fmincon(objfun,mu0,[],[],[],[],[0 0 0 0], [+Inf 0.6 0.6 +Inf],...
   [],options2);
 
bSol = muSolution(1);
pSol = muSolution(2);
rSol = muSolution(3);
dSol = muSolution(4);

[tSolution, solution] = solvingLotkaVolterra(bSol,pSol,rSol,dSol, [min(time), max(time)], IC);

figure(1)
subplot(2,1,1), plot(populationData(:,1),populationData(:,2),'o','MarkerSize',markersize,'MarkerFaceColor','r','MarkerEdgeColor','k')
grid on; grid minor; hold on
plot(tSolution+initialTime, solution(:,1),'b')
xlabel('Time [years]','Interpreter','latex','FontSize',30)
ylabel('Snowshoe Hare Pelts','Interpreter','latex','FontSize',30,'Color','r')

subplot(2,1,2), plot(populationData(:,1), populationData(:,3),'o','MarkerSize',markersize,'MarkerFaceColor','b','MarkerEdgeColor','k')
grid on; grid minor; hold on
plot(tSolution+initialTime, solution(:,2),'b')
xlabel('Time [years]','Interpreter','latex','FontSize',30)
ylabel('Canada Lynx Pelts','Interpreter','latex','FontSize',30,'Color','b')
%sgtitle('Population Data','Interpreter','latex','FontSize',30)
