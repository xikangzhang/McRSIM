function [missrate, grp, index] = SSC_MDD_LA(X, s, camID)
% Inputs:
% X: data matrix
% s: groundtruth labels
% camID: camera ID of each trajectory

K = max(s);
X1 = X(:, camID==1);
X2 = X(:, camID==2);

opt.metric = 'JBLD';
opt.sigma = 10^-4;
opt.H_structure = 'HtH';
opt.H_rows = 10;

[grp1, HHcenter1] = sscjbld(X1, K, opt);
[grp2, HHcenter2] = sscjbld(X2, K, opt);
D = HHdist(HHcenter1, HHcenter2, opt);
[ind] = lapjv(D);
grp(camID==1, :) = grp1;
grp(camID==2, :) = grp2(:, ind);

ind = (s~=0);
grp = (grp(ind,:));
s = s(ind);

[missrate, index] = ErrorRate2(grp, s); % calculate the error rate

end

function [grp, HHcenter] = sscjbld(X, K, opt)

X(3:3:end, :) = [];
r = 0; affine = true; alpha = 800; outlier = false; rho = 0.7;

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

% MDD
features = cell(1, size(X, 2));
for j=1:size(X, 2)
    t = reshape(X(:, j), 2, []);
    v = diff(t,1,2);
    features{j} = v;
end

HH  = getHH(features,opt);
D = HHdist(HH, [], opt);
D = D / max(D(:));
Wj = exp(-D / 1);

W = W .* Wj;
grp = ncutW(W, K);

HHcenter = cell(1, K);
for i = 1:K
    HHcenter{i} = steinMean(cat(3,HH{grp(:,i)==1}));
end

end