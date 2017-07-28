function [val, num, den, eigenValues] = minCutValue(W, label)

l = unique(label);
K = length(l);
num = 0;

for i = 1:K
    idx = (label==l(i));
    Wg = W(~idx,idx);
    cutA = sum(Wg(:));
    num = num+cutA;
end


Lambda = diag(1./sqrt(sum(W)));
L = Lambda * W * Lambda;
L = (L + L') / 2;
L = eye(size(W)) - L;

eigenValues = eig(L);
eigenValues = sort(eigenValues);
den = abs(eigenValues(K+1)-eigenValues(K));

if abs(den) < 1e-3
    val = inf;
else
    val = num / den;
end

end