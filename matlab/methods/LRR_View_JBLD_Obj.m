function [missrate, grp, index] = LRR_View_JBLD_Obj(X, s, camID)

K = max(s);
X1 = X(:, camID==1);
X2 = X(:, camID==2);

opt.metric = 'JBLD';
opt.sigma = 10^-4;
opt.H_structure = 'HtH';
opt.H_rows = 10;

[grp1, HHcenter1] = lrrjbld(X1, K, opt);
[grp2, HHcenter2] = lrrjbld(X2, K, opt);
D = HHdist(HHcenter1, HHcenter2, opt);
[ind] = lapjv(D);
grp(camID==1, :) = grp1;
grp(camID==2, :) = grp2(:, ind);

ind = (s~=0);
grp = (grp(ind,:));
s = s(ind);

[missrate, index] = ErrorRate2(grp, s); % calculate the error rate

end

function [grp, HHcenter] = lrrjbld(X, K, opt)

X(3:3:end-1, :) = []; % retain the last row of ones
lambda = 4;

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

% JBLD
X(end, :) = [];
features = cell(1, size(X, 2));
for j=1:size(X, 2)
    t = reshape(X(:, j), 2, []);
    v = diff(t,1,2);
    features{j} = v;
end

opt.metric = 'JBLD';
opt.sigma = 10^-4;
opt.H_structure = 'HtH';
opt.H_rows = 10;

HH  = getHH(features,opt);
D = HHdist(HH, [], opt);
D = D / max(D(:));
Wj = exp(-D / 1);

L = L .* Wj;

% spectral clustering
D = diag(1./sqrt(sum(L,2)));
L = D*L*D;
[U,S,V] = svd(L);
V = U(:,1:K);
V = D*V;
grps = kmeans(V,K,'emptyaction','singleton','replicates',20,'display','off');
grp = zeros(length(grps), K);
for i = 1:K
    grp(grps==i, i) = 1;
end

HHcenter = cell(1, K);
for i = 1:K
    HHcenter{i} = steinMean(cat(3,HH{grp(:,i)==1}));
end

end