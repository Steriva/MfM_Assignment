clear; clc; close all

addpath('matlab')

T = 4;      % final time
N = 2^9;    % spatial discretization
dt = 5e-3;  % temporal discretization

%% Creating the training data from the KS PDE


nu = 0.04; % IC parameter

[x_train, t_train, u_train] =  solveKS(nu, T, dt, N);

disp(strcat('The KS equation has been solved with nu=',num2str(nu),' and its data are stored into data folder.'))
disp(' ')

% The following lines are used to generate the input and output data to
% train the NN.
% The following matrices are built X_n and X_n+1, so that the time advance
% can be modelled in the NN

input  = u_train(:,1:end-1);
output = u_train(:,2:end);

save('data/KS_trainingData.mat', 'input', 'output', 't_train',"x_train", "u_train", 'nu')

%% Creating the validation data from the KS PDE


nu = 0.1; % IC parameter

[x_test, t_test, u_test] =  solveKS(nu, T, dt, N);


disp(strcat('The KS equation has been solved with nu=',num2str(nu),' and its data are stored into data folder.'))
disp(' ')

% The following lines are used to generate the input and output data to
% train the NN.
% The following matrices are built X_n and X_n+1, so that the time advance
% can be modelled in the NN

input  = u_test(:,1:end-1);
output = u_test(:,2:end);

save('data/KS_testingData.mat', 'input', 'output', 't_test',"x_test", "u_test", 'nu')


