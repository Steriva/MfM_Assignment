function int=traprule(x,y)
% This function allows to perform an integration on a function using the
% trapezoidal rule . 
%               int = traprule ( x , y )
% x and y must be vectors of the same dimensions.

    A=zeros(length(x),1);
    for k=1:length(x)-1
        dx=x(k+1)-x(k);
        A(k+1)=dx/2*(y(k+1)+y(k))+A(k);
    end
    int=A(end);

end