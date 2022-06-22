function norm = L2normOnRectangle(x,y,fun)

% This function can be used to approximate the L2-norm of a function fun 
% (whose format must be a matrix with the same dimension of x and y) on a
% rectangle which is uniformly discretized. 
% The integral is approximated with the constant rule:
%
% int_a^b fun dx = sum_i dx_i * fun(x_i)
% 


dx = x(2)-x(1);
dy = y(2)-x(1);
area = dx*dy;

nx = length(x);
ny = length(y);

norm = 0.;
for ii = 1:nx
    for jj = 1:ny
        norm = norm + area * fun(ii,jj).^2;
    end
end
norm = sqrt(norm);

end