clear; clc; close all

markersize = 10;

% Add path for the optimized DMD
addpath('matlab')
addpath('matlab/optimizedDMD')
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
time = (populationData(:,1)-initialTime)';

figure(1)
subplot(2,1,1), plot(populationData(:,1),populationData(:,2),'o','MarkerSize',markersize,'MarkerFaceColor','r','MarkerEdgeColor','k')
grid on; grid minor
xlabel('Time [years]','Interpreter','latex','FontSize',30)
ylabel('Snowshoe Hare Pelts','Interpreter','latex','FontSize',30,'Color','r')

subplot(2,1,2), plot(populationData(:,1), populationData(:,3),'o','MarkerSize',markersize,'MarkerFaceColor','b','MarkerEdgeColor','k')
grid on; grid minor
xlabel('Time [years]','Interpreter','latex','FontSize',30)
ylabel('Canada Lynx Pelts','Interpreter','latex','FontSize',30,'Color','b')
%sgtitle('Population Data','Interpreter','latex','FontSize',30)

%% Perform SVD
X = (populationData(:,2:3))';
[u,s,v] = svd(X,'econ');

% Plot the eigenvalues sigma_j
figure(2)
plot(diag(s)/sum(diag(s))*100,'b-o','MarkerFaceColor','b','Markersize',markersize,'MarkerEdgeColor','k','Linewidth',1)
grid on; grid minor;
xlabel('Dimension $j$','Interpreter','latex','FontSize',30)
ylabel('$\frac{\sigma_j}{\sum_k \sigma_k}\;\;[\%]$','Interpreter','latex','FontSize',50,'Color','k')

% Plot the v-modes to see the time behaviour
figure(3)
plot(populationData(:,1), v(:,1),'r^-','MarkerFacecolor','r','MarkerEdgeColor','k','Markersize',markersize,'Linewidth', 2)
hold on; grid on; grid minor;
plot(populationData(:,1), v(:,2),'b^-','MarkerFacecolor','b','MarkerEdgeColor','k','Markersize',markersize,'Linewidth', 2)
xlim([min(populationData(:,1)),max(populationData(:,1))])
legend('First', 'Second','Interpreter','latex','FontSize',30,'Location','Best')
xlabel('Time [years]','Interpreter','latex','FontSize',30)
title('$V-$modes','Interpreter','latex','FontSize',50,'Color','k')

%% Perform exact DMD
r = 2;

X1 = X(:,1:end-1); X2 = X(:,2:end);
[U, Sigma, V] = svd(X1,'econ');
U = U(:,1:r);
Sigma = Sigma(1:r,1:r);
V = V(:,1:r);

Atilde = U'*X2*V*diag(1./diag(Sigma)); % similarity transformation A to Atilde

[eV, D] = eig(Atilde); % eigen-decomposition: eV vector, D eigenvalues
mu = diag(D); % mu: diagonal of eigenvalue matrix D
dt = populationData(2,1)-populationData(1,1);
omega = log(mu)/dt;
Phi = X2*(V/Sigma)*eV; % DMD modes
alpha1 = Sigma*V(1,:)';
bj = (eV * D)\alpha1;

u_modes=zeros(size(V,2), length(time));
for iter = 1:length(time)
    u_modes(:,iter) = bj.*exp(omega*time(iter));
end
u_dmd = Phi*u_modes;

% Plot the DMD reconstruction vs true solution
figure(4)
subplot(2,1,1), plot(populationData(:,1),populationData(:,2),'o','MarkerSize',markersize,'MarkerFaceColor','r','MarkerEdgeColor','k')
grid on; grid minor; hold on
plot(populationData(:,1),u_dmd(1,:),'o','MarkerSize',markersize,'MarkerFaceColor','b','MarkerEdgeColor','k')
legend('FOM', 'DMD','Interpreter','latex','FontSize',30,'Location','Best')
xlabel('Time [years]','Interpreter','latex','FontSize',30)
ylabel('Snowshoe Hare Pelts','Interpreter','latex','FontSize',30,'Color','r')

subplot(2,1,2), plot(populationData(:,1),populationData(:,3),'o','MarkerSize',markersize,'MarkerFaceColor','r','MarkerEdgeColor','k')
grid on; grid minor; hold on
plot(populationData(:,1), u_dmd(2,:),'o','MarkerSize',markersize,'MarkerFaceColor','b','MarkerEdgeColor','k')
legend('FOM', 'DMD','Interpreter','latex','FontSize',30,'Location','Best')
xlabel('Time [years]','Interpreter','latex','FontSize',30)
ylabel('Canada Lynx Pelts','Interpreter','latex','FontSize',30,'Color','b')
%sgtitle('Population Data','Interpreter','latex','FontSize',30)

%% Optimized DMD

[w,e1,b] = optdmd(X,time,2,2);

% reconstructed values
x1 = w*diag(b)*exp(e1*time);
relerr_r = norm(x1-X,'fro')/norm(X,'fro');


if max(max(abs(imag(x1))))<1e-8
    x1 = real(x1);
else
    disp('Check the imaginary part: it may be too high')
end

figure(5)
subplot(2,1,1), plot(populationData(:,1),populationData(:,2),'o','MarkerSize',markersize,'MarkerFaceColor','r','MarkerEdgeColor','k')
grid on; grid minor; hold on
plot(time(1:end)+initialTime,x1(1,:),'b-o','Linewidth',1.5,'MarkerSize',markersize,'MarkerFaceColor','b','MarkerEdgeColor','k')
legend('FOM', 'opt-DMD','Interpreter','latex','FontSize',30,'Location','Best')
xlabel('Time [years]','Interpreter','latex','FontSize',30)
ylabel('Snowshoe Hare Pelts','Interpreter','latex','FontSize',30,'Color','k')

subplot(2,1,2), plot(populationData(:,1),populationData(:,3),'o','MarkerSize',markersize,'MarkerFaceColor','r','MarkerEdgeColor','k')
grid on; grid minor; hold on
plot(time(1:end)+initialTime,x1(2,:),'b-o','Linewidth',1.5,'MarkerSize',markersize,'MarkerFaceColor','b','MarkerEdgeColor','k')
legend('FOM', 'opt-DMD','Interpreter','latex','FontSize',30,'Location','Best')
xlabel('Time [years]','Interpreter','latex','FontSize',30)
ylabel('Canada Lynx Pelts','Interpreter','latex','FontSize',30,'Color','k')
sgtitle(strcat('r=',num2str(r)),'Interpreter','latex','FontSize',30)