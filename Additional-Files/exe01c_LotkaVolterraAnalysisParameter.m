clear; clc; close all

addpath('matlab')

colours = struct2cell(colori());
markersize = 10;

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

bSol = 1;
pSol = 0.01;
rSol = 0.01;
dSol = 1;

% [tSolution, solution] = solvingLotkaVolterra(bSol,pSol,rSol,dSol, [min(time), max(time)], IC);
% 
% figure(1)
% subplot(2,1,1), plot(populationData(:,1),populationData(:,2),'o','MarkerSize',markersize,'MarkerFaceColor','r','MarkerEdgeColor','k')
% grid on; grid minor; hold on
% plot(tSolution+initialTime, solution(:,1),'b')
% xlabel('Time [years]','Interpreter','latex','FontSize',30)
% ylabel('Snowshoe Hare Pelts','Interpreter','latex','FontSize',30,'Color','r')
% 
% subplot(2,1,2), plot(populationData(:,1), populationData(:,3),'o','MarkerSize',markersize,'MarkerFaceColor','b','MarkerEdgeColor','k')
% grid on; grid minor; hold on
% plot(tSolution+initialTime, solution(:,2),'b')
% xlabel('Time [years]','Interpreter','latex','FontSize',30)
% ylabel('Canada Lynx Pelts','Interpreter','latex','FontSize',30,'Color','b')
% %sgtitle('Population Data','Interpreter','latex','FontSize',30)

%% Varying b

IC = [populationData(1,2); populationData(1,3)];

b = linspace(0, 10,1e2);
t = linspace(min(time), max(time),200);
[xPlot, yPlot] = meshgrid(b,t);
p = 0.01;
r = 0.01;
d = 1;

prey_b = zeros(length(t),length(b));
predator_b = zeros(length(t),length(b));

for ii = 1 : length(b)
   [~, solution] = solvingLotkaVolterra(b(ii),p,r,d, t, IC);
   prey_b(:,ii) = solution(:, 1);
   predator_b(:,ii) = solution(:, 2);
end

figure
subplot(2,1,1), surf(xPlot,yPlot,prey_b), shading interp, colormap(turbo), colorbar
ylabel('time')
xlabel('b')

subplot(2,1,2), surf(xPlot,yPlot,predator_b), shading interp, colormap(turbo), colorbar
ylabel('time')
xlabel('b')

%% Varying p

IC = [populationData(1,2); populationData(1,3)];

p = linspace(1e-4, 1,1e2);
t = linspace(min(time), max(time),200);
[xPlot, yPlot] = meshgrid(p,t);
b = 1;
r = 0.01;
d = 1;

prey_b = zeros(length(t),length(p));
predator_b = zeros(length(t),length(p));

for ii = 1 : length(p)
   [~, solution] = solvingLotkaVolterra(b,p(ii),r,d, t, IC);
   prey_b(:,ii) = solution(:, 1);
   predator_b(:,ii) = solution(:, 2);
end

figure
subplot(2,1,1), surf(xPlot,yPlot,prey_b), shading interp, colormap(turbo), colorbar
ylabel('time')
xlabel('p')

subplot(2,1,2), surf(xPlot,yPlot,predator_b), shading interp, colormap(turbo), colorbar
ylabel('time')
xlabel('p')

%% Varying r

IC = [populationData(1,2); populationData(1,3)];

r = linspace(1e-4, 1,1e2);
t = linspace(min(time), max(time),200);
[xPlot, yPlot] = meshgrid(p,t);
b = 1;
p = 0.01;
d = 1;

prey_b = zeros(length(t),length(r));
predator_b = zeros(length(t),length(r));

for ii = 1 : length(r)
   [~, solution] = solvingLotkaVolterra(b,p,r(ii),d, t, IC);
   prey_b(:,ii) = solution(:, 1);
   predator_b(:,ii) = solution(:, 2);
end

figure
subplot(2,1,1), surf(xPlot,yPlot,prey_b), shading interp, colormap(turbo), colorbar
ylabel('time')
xlabel('r')

subplot(2,1,2), surf(xPlot,yPlot,predator_b), shading interp, colormap(turbo), colorbar
ylabel('time')
xlabel('r')

%% Varying d

IC = [populationData(1,2); populationData(1,3)];

d = linspace(0, 1,1e2);
t = linspace(min(time), max(time),200);
[xPlot, yPlot] = meshgrid(d,t);
b = 1;
p = 0.01;
r = 0.01;

prey_b = zeros(length(t),length(d));
predator_b = zeros(length(t),length(d));

for ii = 1 : length(d)
   [~, solution] = solvingLotkaVolterra(b,p,r,d(ii), t, IC);
   prey_b(:,ii) = solution(:, 1);
   predator_b(:,ii) = solution(:, 2);
end

figure
subplot(2,1,1), surf(xPlot,yPlot,prey_b), shading interp, colormap(turbo), colorbar
ylabel('time')
xlabel('r')

subplot(2,1,2), surf(xPlot,yPlot,predator_b), shading interp, colormap(turbo), colorbar
ylabel('time')
xlabel('r')