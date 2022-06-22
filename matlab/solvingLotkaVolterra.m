function [tSolution, solution] = solvingLotkaVolterra(b,p,r,d, timeSave, IC)

% This simple function is able to solve the LotkaVolterra equations, given
% the necessary coefficents and the IC.
% The solution will be output at the time given by timeSave

ode_options = odeset('RelTol',1e-12, 'AbsTol',1e-11);
[tSolution,solution] = ode23(@(t,y) LotkaVolterra(y,b,p,r,d),timeSave,IC,ode_options);

    function dydt = LotkaVolterra(y,b,p,r,d)
    
    dydt = [(b-p*y(2))*y(1); (r*y(1)-d)*y(2)];
    end

end