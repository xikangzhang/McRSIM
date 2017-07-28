function d = LERM(X,Y)

[VX, DX] = eig(X);
X2 = VX * diag(log(diag(DX))) * VX';
X2 = (X2 + X2') / 2;

[VY, DY] = eig(Y);
Y2 = VY * diag(log(diag(DY))) * VY';
Y2 = (Y2 + Y2') / 2;

d = norm(X2 - Y2, 'fro');

% d = norm(logm(X) - logm(Y), 'fro');

end
