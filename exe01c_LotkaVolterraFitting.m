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

option = optimset("PlotFcn","optimplotfvalconstr",'Maxiter',200);
objfun = @(mu) computingOF(mu, populationData(:,2:3), time, IC);

mu0 = [0.7; 0.01; 0.01; 0.7];

% [0.7; 0.03; 0.02; 0.8];

% [muSolution,~] = fmincon(objfun,mu0,[],[],[],[],[0 0 0 0.5], [2 +Inf +Inf +Inf],...
%    [],options2);
 
muSolution = fminsearch(objfun,mu0, option);

bSol = muSolution(1);
pSol = muSolution(2);
rSol = muSolution(3);
dSol = muSolution(4);

[tSolution, solution] = solvingLotkaVolterra(bSol,pSol,rSol,dSol, [min(time), max(time)], IC);

figure(1)
subplot(2,1,1), plot(populationData(:,1),populationData(:,2),'o','MarkerSize',markersize,'MarkerFaceColor',colours{1},'MarkerEdgeColor','k')
grid on; grid minor; hold on
plot(tSolution+initialTime, solution(:,1),'Color',colours{1},'LineStyle','--','LineWidth',2)
xlabel('Time [years]','Interpreter','latex','FontSize',30)
ylabel('Snowshoe Hare Pelts','Interpreter','latex','FontSize',30,'Color',colours{1})
legend('Data','Fitting model','Interpreter','latex','FontSize',30)

subplot(2,1,2), plot(populationData(:,1),populationData(:,3),'o','MarkerSize',markersize,'MarkerFaceColor',colours{2},'MarkerEdgeColor','k')
grid on; grid minor; hold on
plot(tSolution+initialTime, solution(:,2),'Color',colours{2},'LineStyle','--','LineWidth',2)
xlabel('Time [years]','Interpreter','latex','FontSize',30)
ylabel('Snowshoe Hare Pelts','Interpreter','latex','FontSize',30,'Color',colours{2})
legend('Data','Fitting model','Interpreter','latex','FontSize',30)
%sgtitle('Population Data','Interpreter','latex','FontSize',30)
