function [Xm] = logEucMean(varargin)

n = nargin;
S = 0;
for i = 1:n
    [V, D] = eig(varargin{i});
    X = V * diag(log(diag(D))) * V';
    X = (X + X') / 2;
    S = S + X;
end
S = S / n;
[VS, DS] = eig(S);
Xm = VS * diag(exp(diag(DS))) * VS';
Xm = (Xm + Xm') / 2;