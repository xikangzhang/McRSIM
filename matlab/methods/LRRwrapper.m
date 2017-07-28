function [missrate, label, gt] = LRRwrapper(X, s, verbose)

if nargin < 3
    X(3:3:end-1, :) = []; % retain the last row of ones
end

lambda = 4;
K = max(s);

%run lrr
Z = solve_lrr(X,lambda);
%post processing
[U,S,V] = svd(Z,'econ');
S = diag(S);
r = sum(S>1e-4*S(1));
U = U(:,1:r);S = S(1:r);
U = U*diag(sqrt(S));
U = normr(U);
L = (U*U').^4;
% spectral clustering
D = diag(1./sqrt(sum(L,2)));
L = D*L*D;
[U,S,V] = svd(L);
V = U(:,1:K);
V = D*V;
grps = kmeans(V,K,'emptyaction','singleton','replicates',20,'display','off');
[miss, index] = missclassGroups(grps,s,K);
missrate =  miss/length(grps);
[~, labelIndex] = sort(index);
label = labelIndex(grps)';
gt = s;

end