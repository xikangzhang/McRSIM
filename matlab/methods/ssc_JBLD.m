function [missrate, grp, W, index] = ssc_JBLD(X, s, verbose)

if nargin < 3
    X(3:3:end, :) = [];
end
r = 0; affine = true; alpha = 800; outlier = false; rho = 0.7;

n = max(s);
Xp = DataProjection(X,r);

if (~outlier)
    CMat = admmLasso_mat_func(Xp,affine,alpha);
    C = CMat;
else
    CMat = admmOutlier_mat_func(Xp,affine,alpha);
    N = size(Xp,2);
    C = CMat(1:N,:);
end

W = BuildAdjacency(thrC(C,rho));


% JBLD
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

HHt  = getHH(features,opt);
D = HHdist(HHt, [], opt);
D = D / max(D(:));
Wj = exp(-D / 1);

W = W .* Wj;


grp = SpectralClustering(W, n);
grp = bestMap(s,grp);
missrate = sum(s(:) ~= grp(:)) / length(s);
index = 1:max(s);

% [grp,~,~] = ncutW(W,n);
% ind = (s~=0);
% grp = (grp(ind,:));
% s = s(ind);
% [missrate, index] = ErrorRate2(grp, s); % calculate the error rate

end