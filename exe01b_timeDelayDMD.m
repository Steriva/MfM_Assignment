clear; clc; close all

markersize = 10;

% Add path for auxiliary code
addpath('matlab')
addpath('matlab/optimizedDMD')
colours = struct2cell(colori());

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
sgtitle('Population Data','Interpreter','latex','FontSize',30)

%% Perform SVD at different values of p

p = 5:5:25;
X = (populationData(:,2:3))';

H = cell(length(p),1);
Legend = cell(length(p),1);
figure(2)

for jj = 1:length(p)

    H{jj} = buildHankelMatrix(X, p(jj));
    
    [u,s,v] = svd(H{jj},'econ');
    
    % Plot the eigenvalues sigma_j
    subplot(1,2,1), plot(diag(s)/sum(diag(s))*100,'-o','Color',colours{jj},...
        'MarkerFaceColor',colours{jj},'Markersize',markersize,'MarkerEdgeColor','k','Linewidth',1)
    hold on
     
    subplot(1,2,2), semilogy(diag(s)/sum(diag(s)),'-^','Color',...
        colours{jj},'MarkerFaceColor',colours{jj},'Markersize',markersize,'MarkerEdgeColor','k','Linewidth',1);
    hold on
    
    Legend{jj} = strcat( '$p= ',num2str(p(jj)),'$' );
end

subplot(1,2,1), grid on; grid minor; 
xlabel('Dimension $j$','Interpreter','latex','FontSize',30)
ylabel('$100\cdot\frac{\sigma_j}{\sum_k \sigma_k}\;\;[\%]$','Interpreter','latex','FontSize',50,'Color','k')
title('Linear Scale','Interpreter','latex','FontSize',30)
legend(Legend,'Interpreter','latex','FontSize',50)

subplot(1,2,2), grid on; grid minor; 
xlabel('Dimension $j$','Interpreter','latex','FontSize',30)
ylabel('$\frac{\sigma_j}{\sum_k \sigma_k}\;\;[-]$','Interpreter','latex','FontSize',50,'Color','k')    
title('Log Scale','Interpreter','latex','FontSize',30)
legend(Legend,'Interpreter','latex','FontSize',50)

%% Perform exact DMD

r = 10;
pOpt = size(X,2)/(size(X,1)+1);
H1 = H{p==pOpt}(:,1:end-1);
H2 = H{p==pOpt}(:,2:end);

[U, Sigma, V] = svd(H1,'econ');
U = U(:,1:r);
Sigma = Sigma(1:r,1:r);
V = V(:,1:r);

Atilde = U'*H2*V*diag(1./diag(Sigma)); % similarity transformation A to Atilde

[eV, D] = eig(Atilde); % eigen-decomposition: eV vector, D eigenvalues
mu = diag(D); % mu: diagonal of eigenvalue matrix D
dt = populationData(2,1)-populationData(1,1);
omega = log(mu)/dt;
Phi = H2*(V/Sigma)*eV; % DMD modes
alpha1 = Sigma*V(1,:)';
bj = (eV * D)\alpha1;

figure(3)
plot(mu,'o','MarkerSize',markersize,'MarkerFaceColor','r','MarkerEdgeColor','k')
grid on; grid minor; hold on
plot([0 0], [-1 1],[-1 1], [0 0],'Linewidth',2,'Color','k')
xlabel('Re','Interpreter','latex','FontSize',30)
ylabel('Im','Interpreter','latex','FontSize',30)


u_modes=zeros(size(V,2), length(time));
for iter = 1:length(time)
    u_modes(:,iter) = bj.*exp(omega*time(iter));
end
u_dmd = Phi*u_modes;

if max(max(abs(imag(u_dmd))))<1e-12
    u_dmd = real(u_dmd);
else
    disp('Check the imaginary part: it may be too high')
end

sampleTime = linspace(min(time),max(time),1e3);
u_modesSample=zeros(size(V,2), length(sampleTime));
for iter = 1:length(sampleTime)
    u_modesSample(:,iter) = bj.*exp(omega*sampleTime(iter));
end
u_dmdSample = Phi*u_modesSample;

if max(max(abs(imag(u_dmdSample))))<1e-12
    u_dmdSample = real(u_dmdSample);
else
    disp('Check the imaginary part: it may be too high')
end

% Plot the DMD reconstruction vs true solution
figure(4)
subplot(2,1,1), plot(populationData(:,1),populationData(:,2),'o','MarkerSize',markersize,'MarkerFaceColor','r','MarkerEdgeColor','k')
grid on; grid minor; hold on
plot(populationData(:,1),u_dmd(1,:),'o','MarkerSize',markersize,'MarkerFaceColor','b','MarkerEdgeColor','k')
plot(sampleTime+initialTime,u_dmdSample(1,:),'Color','b','Linewidth',1.5)
legend('FOM', 'DMD','Interpreter','latex','FontSize',30,'Location','Best')
xlabel('Time [years]','Interpreter','latex','FontSize',30)
ylabel('Snowshoe Hare Pelts','Interpreter','latex','FontSize',30,'Color','k')

subplot(2,1,2), plot(populationData(:,1),populationData(:,3),'o','MarkerSize',markersize,'MarkerFaceColor','r','MarkerEdgeColor','k')
grid on; grid minor; hold on
plot(populationData(:,1), u_dmd(2,:),'o','MarkerSize',markersize,'MarkerFaceColor','b','MarkerEdgeColor','k')
plot(sampleTime+initialTime,u_dmdSample(2,:),'Color','b','Linewidth',1.5)
legend('FOM', 'DMD','Interpreter','latex','FontSize',30,'Location','Best')
xlabel('Time [years]','Interpreter','latex','FontSize',30)
ylabel('Canada Lynx Pelts','Interpreter','latex','FontSize',30,'Color','k')
sgtitle(strcat('r=',num2str(r)),'Interpreter','latex','FontSize',30)


%% Optimized DMD

[w,e1,b] = optdmd(H{p==pOpt},time(1:end-pOpt),r,1);

% reconstructed values
x1 = w*diag(b)*exp(e1*time(1:end-pOpt));
relerr_r = norm(x1-H{p==pOpt},'fro')/norm(H{p==pOpt},'fro');

if max(max(abs(imag(x1))))<1e-8
    x1 = real(x1);
else
    disp('Check the imaginary part: it may be too high')
end

u_optDMD = [x1(1,:), x1(end-1,end-pOpt+1:end);...
            x1(2,:), x1(end,  end-pOpt+1:end)];

figure(5)
subplot(2,1,1), plot(populationData(:,1),populationData(:,2),'o','MarkerSize',markersize,'MarkerFaceColor','r','MarkerEdgeColor','k')
grid on; grid minor; hold on
plot(time(1:end)+initialTime,u_optDMD(1,:),'b-o','Linewidth',1.5,'MarkerSize',markersize,'MarkerFaceColor','b','MarkerEdgeColor','k')
legend('FOM', 'opt-DMD','Interpreter','latex','FontSize',30,'Location','Best')
xlabel('Time [years]','Interpreter','latex','FontSize',30)
ylabel('Snowshoe Hare Pelts','Interpreter','latex','FontSize',30,'Color','k')

subplot(2,1,2), plot(populationData(:,1),populationData(:,3),'o','MarkerSize',markersize,'MarkerFaceColor','r','MarkerEdgeColor','k')
grid on; grid minor; hold on
plot(time(1:end)+initialTime,u_optDMD(2,:),'b-o','Linewidth',1.5,'MarkerSize',markersize,'MarkerFaceColor','b','MarkerEdgeColor','k')
legend('FOM', 'opt-DMD','Interpreter','latex','FontSize',30,'Location','Best')
xlabel('Time [years]','Interpreter','latex','FontSize',30)
ylabel('Canada Lynx Pelts','Interpreter','latex','FontSize',30,'Color','k')
sgtitle(strcat('r=',num2str(r)),'Interpreter','latex','FontSize',30)

