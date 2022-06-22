clear; clc; close all

addpath('matlab')

if isfile('data/reaction_diffusion_big.mat')
    load('data/reaction_diffusion_big.mat')
else
    check=input('The Reaction diffusion equation will be solved: do you want to continue? [Y/N]\n','s');
    if strcmp(check,'Y')
        [t,x,y,u,v] = reaction_diffusion();
        save('data/reaction_diffusion_big.mat','t','x','y','u','v')
    else
        return
    end
end
    
markersize = 10;

%% Reshaping the data

n = length(x);
m = length(t);
dt = t(2)-t(1);
dx = x(2)-x(1);
dy = y(2)-y(1);

% The X matrix will be (2n^2 x m) in which the spatial behaviour of u and v
% has been rearranged as column vector
Xu = zeros(n*n, m);
Xv = zeros(n*n, m);

for jj = 1:m
textwaitbar(jj, m, 'Reshaping full order data');
Xu(:,jj) = reshape(u(:,:,jj), [n*n, 1]);
Xv(:,jj) = reshape(v(:,:,jj), [n*n, 1]);
end



%% Performing SVD

disp('Performing the SVD')
[Uu, Su, Vu] = svd(Xu,'econ');
[Uv, Sv, Vv] = svd(Xv,'econ');

% Plot the eigenvalues sigma_j
figure(1)
subplot(2,2,1)
semilogy(0:length(diag(Su))-1, diag(Su)/sum(diag(Su)),'b-o','MarkerFaceColor','b','Markersize',markersize,'MarkerEdgeColor','k','Linewidth',1)
grid on; grid minor;
xlim([0 100])
%xlabel('Dimension $j$','Interpreter','latex','FontSize',30)
ylabel('$\frac{\sigma^{(u)}_j}{\sum_k \sigma_k^{(u)}}\;\;[-]$','Interpreter','latex','FontSize',50,'Color','k')

subplot(2,2,2)
plot(0:length(diag(Su))-1, 1-diag(Su)/sum(diag(Su)),'b-o','MarkerFaceColor','b','Markersize',markersize,'MarkerEdgeColor','k','Linewidth',1)
grid on; grid minor;
xlim([0 20])
%xlabel('Dimension $j$','Interpreter','latex','FontSize',30)
ylabel('$1-\frac{\sigma_j^{(u)}}{\sum_k \sigma_k^{(u)}}\;\;[-]$','Interpreter','latex','FontSize',50,'Color','k')

subplot(2,2,3)
semilogy(0:length(diag(Sv))-1, diag(Sv)/sum(diag(Sv)),'b-o','MarkerFaceColor','b','Markersize',markersize,'MarkerEdgeColor','k','Linewidth',1)
grid on; grid minor;
xlim([0 100])
xlabel('Dimension $j$','Interpreter','latex','FontSize',30)
ylabel('$\frac{\sigma^{(v)}_j}{\sum_k \sigma_k^{(v)}}\;\;[-]$','Interpreter','latex','FontSize',50,'Color','k')

subplot(2,2,4)
plot(0:length(diag(Sv))-1, 1-diag(Sv)/sum(diag(Sv)),'b-o','MarkerFaceColor','b','Markersize',markersize,'MarkerEdgeColor','k','Linewidth',1)
grid on; grid minor;
xlim([0 20])
xlabel('Dimension $j$','Interpreter','latex','FontSize',30)
ylabel('$1-\frac{\sigma_j^{(v)}}{\sum_k \sigma_k^{(v)}}\;\;[-]$','Interpreter','latex','FontSize',50,'Color','k')

%% L2-error on the training set - for u - 
rVec = [50, 40, 30, 25, 20, 15, 10, 5, 1];
globalErrorL2norm = zeros(length(rVec),2);

disp(' ')

for kk = 1:length(rVec)

    r = rVec(kk);
    disp(strcat('Computing L2-error with r=',num2str(r)))
    Utmp = Uu(:,1:r); Stmp = Su(1:r, 1:r); Vtmp = Vu(:,1:r);
    
    lowRankData = Utmp*Stmp*Vtmp';
    
    uLowRank = zeros(size(u));
%   vLowRank = zeros(size(v));
    for jj = 1:m
        uLowRank(:,:,jj) = reshape(lowRankData(1:end,jj), [n, n, 1]);
        % vLowRank(:,:,jj) = reshape(lowRankData(1:end,jj), [n, n, 1]);
    end
    
    errorL2norm = zeros(m,1);
        
    for jj = 1:m
        errorL2norm(jj,1) = L2normOnRectangle(x,y, u(:,:,jj)-uLowRank(:,:,jj)) ...
                          / L2normOnRectangle(x,y, u(:,:,jj) );
%         errorL2norm(jj,2) = L2normOnRectangle(x,y, v(:,:,jj)-vLowRank(:,:,jj)) ...
%                           / L2normOnRectangle(x,y, v(:,:,jj) );
    end
    
    globalErrorL2norm(kk,:) = mean(errorL2norm);

end

figure(2)
semilogy(rVec, globalErrorL2norm(:,1),'r-^','MarkerFaceColor','r','Markersize',markersize,'MarkerEdgeColor','k','Linewidth',1.5)
grid on; grid minor; 
xlabel('Rank $r$','Interpreter','latex','FontSize',30)
ylabel('Relative error: $ \left\langle\frac{|| u-\mathcal{P}_r(u)||_{L^2}}{|| u||_{L^2}}\right\rangle_t $','Interpreter','latex','FontSize',30)

%% Training the NN in low rank variables

r = 5;

disp(' ')
disp('Starting the training of the NNs')

% train the NN for u
input  = (Vu(1:end-1,1:r))'; 
output = (Vu(2:end,1:r))'; 
net_u = feedforwardnet([10 10 5]);
net_u.layers{1}.transferFcn = 'purelin';
net_u.layers{2}.transferFcn = 'purelin';
net_u.layers{3}.transferFcn = 'purelin';
net_u = train(net_u,input,output);
 
% % train the NN for v
% input  = (Vv(1:end-1,1:r))'; 
% output = (Vv(2:end,1:r))'; 
% net_v = feedforwardnet([5 5 5]);
% net_v.layers{1}.transferFcn = 'logsig';
% net_v.layers{2}.transferFcn = 'radbas';
% net_v.layers{3}.transferFcn = 'purelin';
% net_v = train(net_v,input,output);

%% Testing the NN prediction - only for u 

V_NN = zeros(m,r);
V_NN(1,:) = Vu(1,1:r);

for jj=2:m
    textwaitbar(jj, m, 'Predicting V - for u - through NN');
    tmpIN = V_NN(jj-1,:)';
    tmpOUT=net_u(tmpIN);
    V_NN(jj,:)=tmpOUT.'; 
    tmpIN=tmpOUT;
end


V_L2error = zeros(size(V_NN,2),1);

figure(3)
for jj = 1:4
    subplot(2,2,jj), plot(t,V_NN(:,jj),'r')
    hold on; grid on; grid minor
    plot(t, Vu(:,jj),'b')
    xlabel('Time $t$','Interpreter','latex','FontSize',20)
    ylabel(strcat('$V_{',num2str(jj),'}^{(u)}$ mode'),'Interpreter','latex','FontSize',20)
    legend('Predicted','True','Interpreter','latex','FontSize',20)

    V_L2error(jj) = traprule(t, (V_NN(:,jj) -Vu(:,jj)).^2 ) ...
                  / traprule(t,Vu(:,jj).^2);
    title(strcat('Relative error: $\frac{|| V^{(u)}_{',num2str(jj),'}-\mathcal{P}_r(V^{(u)}_{',num2str(jj),'})||_{L^2}}{|| V^{(u)}_{',num2str(jj),'}||_{L^2}} =$',num2str(V_L2error(jj))),...
        'Interpreter','latex','FontSize',20)
end


Xpred = Uu(:,1:r) * Su(1:r,1:r) * V_NN';

u_NN = zeros(size(u));

errorL2norm = zeros(m,1);

for jj = 1:m
    u_NN(:,:,jj) = reshape(Xpred(:,jj), [n, n, 1]);

    errorL2norm(jj,1) = L2normOnRectangle(x,y, u(:,:,jj)-u_NN(:,:,jj)) ...
                      / L2normOnRectangle(x,y, u(:,:,jj) );
end

ave_errorL2norm = mean(errorL2norm);

% figure(4)
% semilogy(t, errorL2norm(:,1),'r-','Linewidth',1.5)
% grid on; grid minor; 
% xlabel('Time $t$','Interpreter','latex','FontSize',30)
% ylabel('Relative error: $ \frac{|| u-\mathcal{P}_r(u)||_{L^2}}{|| u||_{L^2}} $','Interpreter','latex','FontSize',30)


[xPlot, yPlot] = meshgrid(x,y);

% figure(5)
% for jj = 1:m
%     subplot(1,2,1), contourf(xPlot, yPlot, u(:,:,jj))
%     title(strcat('$t=', num2str(dt*jj),'$'),'Interpreter','latex','FontSize',20)
%     xlabel('$x$','Interpreter','latex','FontSize',20)
%     ylabel('$y$','Interpreter','latex','FontSize',20)
%     colorbar; clim([-1 1])
%     colormap("turbo")
% 
%     subplot(1,2,2), contourf(xPlot, yPlot, u_NN(:,:,jj))
%     title(strcat('$t=', num2str(dt*jj),'$'),'Interpreter','latex','FontSize',20)
%     xlabel('$x$','Interpreter','latex','FontSize',20)
%     ylabel('$y$','Interpreter','latex','FontSize',20)
%     colorbar; clim([-1 1])
%     colormap("turbo")
%     pause(0.01)
% end

figure(5)
for jj = 1:4
    timePlot = dt * jj * 50;
    subplot(2,4,2*(jj-1)+1), contourf(xPlot, yPlot, u(:,:,timePlot/dt))
    if 2*(jj-1)+1 == 1 || 2*(jj-1)+1 == 5
        ylabel('$y$','Interpreter','latex','FontSize',20)
    end
    title(strcat('True - $t=', num2str(timePlot),'$'),'Interpreter','latex','FontSize',20)
    colorbar; clim([-1.5 1.5])
    colormap("turbo")

    subplot(2,4,2*(jj)), contourf(xPlot, yPlot, u_NN(:,:,timePlot/dt))
    if jj>2
        xlabel('$x$','Interpreter','latex','FontSize',20)
    end
    title(strcat('NN - $t=', num2str(timePlot),'$'),'Interpreter','latex','FontSize',20)
    colorbar; clim([-1.5 1.5])
    colormap("turbo")
end
