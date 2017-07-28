function d = AIRM(X, Y)

% [U,S,V] = svd(X);
% s = diag(S);
% q = 1./sqrt(s);
% Q = diag(q);
% Z = U*Q*V';
% Z = (Z+Z')/2;
% d = norm(logm(Z*Y*Z),'fro');

d = dist(X, Y);

end