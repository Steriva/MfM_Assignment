function H = buildHankelMatrix(X, p)

% This function can be used to build the Hankel matrix starting from the
% measuremnt matrix X (n x M).
% The parameter p can be used to defined the desired number of columns for
% H, which is (p*n) x (M-p).

    [n,M] = size(X);
    H = zeros(p*n, M-p);
    
    for ii = 1:p
        partialStateSpaceVector = [X(:,ii:M-p+(ii-1))];
        H(2*(ii-1)+1,:) = partialStateSpaceVector(1,:);
        H(2*ii,:) = partialStateSpaceVector(2,:);
    end

end