function X1 = steinMean(X)


n = size(X,3);

X0 = 0;
for i=1:n
    X0 = X0 + X(:,:,i);
end
iY = X0./n;

iY0 = X(:,:,1);
% while norm(iY-iY0)>1e-15
while JBLD(iY,iY0)>1e-15
    iY0 = iY;
    Y = 0;
    for i=1:n
        Y = Y + inv((iY0+X(:,:,i))/2);
    end
    Y = Y./n;
    iY = inv(Y);
end

X1 = iY;

end

function T_inv = inverse(T)

[V, D] = schur(T);
d = diag(D);
T_inv = V * diag(1./d) * V';
T_inv = (T_inv + T_inv') / 2;

end