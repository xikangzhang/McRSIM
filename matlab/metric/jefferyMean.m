function Xm = jefferyMean(varargin)

n = nargin;
P = 0;
Q = 0;
for i = 1:n
    P = P + inv(varargin{i});
    Q = Q + varargin{i};
end
[V, D] = eig(P);
d = diag(D).^0.5;
M = V * diag(d) * V';
M = (M + M') / 2;
N = V * diag(1./d) * V';
N = (N + N') / 2;
S = M * Q * M;
S = (S + S') /2;
[VS, DS] = eig(S);
ds = diag(DS).^0.5;
T = VS * diag(ds) * VS';
T = (T + T') / 2;
Xm = N * T * N;
Xm = (Xm + Xm') / 2;

end