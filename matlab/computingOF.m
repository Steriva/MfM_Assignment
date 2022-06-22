function J = computingOF (mu, data, timeSave, IC)

% This auxiliary function is used to compute the Objective Function of the
% minimization problem, given the data n x 2 (the first column are the prey
% and the second the predator). timeSave is a vector with the time at which
% the data are collected and it is used to interpolate the results of the
% mathematical model.

% The equations are
% dx/dt = (b - p * y) * x
% dy/dt = (r * x - d) * y
% the parameters are all positive.
    bTest = mu(1);
    pTest = mu(2);
    rTest = mu(3);
    dTest = mu(4);

    [~, solution] = solvingLotkaVolterra(bTest,pTest,rTest,dTest, timeSave, IC);

    J = sum( (data(:,1)-solution(:,1)).^2 + (data(:,2)-solution(:,2)).^2 );
 
end