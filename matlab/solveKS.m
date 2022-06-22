
function [x, tt, uu] = solveKS(nu, tmax, h, N)

% This function is called to solve the  Kuramoto-Sivashinsky equation 
% (from Trefethen)
% u_t = -u*u_x - u_xx - nu*u_xxxx,  periodic BCs
% 
% The solution is output along with the spatial and temporal variable 
% (discretized) and the solution itself. 
% 'nu' is a parameter that is able to change the IC (default = 87).
% 'tmax' is the ending time
% 'h' is the temporale discretization (NOTE: 10^-2 or lower are preferred)
% 'N' is the spatial discretization (NOTE: use a power of 2)
    
    x = 2*pi*(1:N)'/N;
    u = -sin(x)+2*cos(2*x)+3*cos(3*x)-4*sin(4*x);
    v = fft(u);
    
    % % % % % %
    %Spatial grid and initial condition:
    k = [0:N/2-1 0 -N/2+1:-1]';
    L = k.^2 - nu*k.^4; %% added the 4* for my specif
    E = exp(h*L); E2 = exp(h*L/2);
    M = 16;
    r = exp(1i*pi*((1:M)-.5)/M);
    LR = h*L(:,ones(M,1)) + r(ones(N,1),:);
    Q = h*real(mean( (exp(LR/2)-1)./LR ,2));
    f1 = h*real(mean( (-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3 , 2) );
    f2 = h*real(mean( (2+LR+exp(LR).*(-2+LR))./LR.^3 ,2));
    f3 = h*real(mean( (-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3 ,2));
    
    % Main time-stepping loop:
    nmax = round(tmax/h); g = -0.5i*k;
    tt = zeros(1,nmax);
    uu = zeros(N,nmax);

    for n = 1:nmax
        textwaitbar(n, nmax, 'Solving the PDE');
        t = n*h;
        Nv = g.*fft(real(ifft(v)).^2);
        a = E2.*v + Q.*Nv;
        Na = g.*fft(real(ifft(a)).^2); 
        b = E2.*v + Q.*Na;
        Nb = g.*fft(real(ifft(b)).^2);
        c = E2.*a + Q.*(2*Nb-Nv);
        Nc = g.*fft(real(ifft(c)).^2);
        v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3;     
        
        % Saving the solution
        u = real(ifft(v));
        uu(:,n) = u; 
        tt(n) = t;
    end

end