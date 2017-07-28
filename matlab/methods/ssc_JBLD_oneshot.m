function [missrate, grp, W, index] = ssc_JBLD_oneshot(X, s, camID)

X(3:3:end, :) = [];
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

% KNN on W12 and W21
kNN = 5;
W1 = Wj(camID==1, camID==2);
for j = 1:size(W1, 1)
    [~,ind] = sort(W1(j,:),'descend');
    W1(j, ind(kNN+1:end)) = 0;
end
W2 = Wj(camID==1, camID==2);
for j = 1:size(W2, 2)
    [~,ind] = sort(W2(:,j),'descend');
    W2(ind(kNN+1:end),j) = 0;
end
W3 = (W1 + W2) / 2;
%     W3 = min(W1, W2);
%     W3 = W3.^30;
W(camID==1, camID==2) = W3;
W(camID==2, camID==1) = W3';

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