function [missrate, grp, index] = RSIM_MDD_LA(X, s, UpperD, LowerD, camID)

K = max(s);
X1 = X(:, camID==1);
X2 = X(:, camID==2);

opt.metric = 'JBLD';
opt.sigma = 10^-4;
opt.H_structure = 'HtH';
opt.H_rows = 10;

[grp1, HHcenter1] = rsimjbld(X1, K, UpperD, LowerD, opt);
[grp2, HHcenter2] = rsimjbld(X2, K, UpperD, LowerD, opt);
D = HHdist(HHcenter1, HHcenter2, opt);
[ind] = lapjv(D);
grp(camID==1, :) = grp1;
grp(camID==2, :) = grp2(:, ind);

ind = (s~=0);
grp = (grp(ind,:));
s = s(ind);

[missrate, index] = ErrorRate2(grp, s); % calculate the error rate

end

function [grp, HHcenter] = rsimjbld(X, K, UpperD, LowerD, opt)

if(nargin<4)
	LowerD = 1;
end
if(nargin<3)
	UpperD = 4;
end

r = LowerD*K:UpperD*K; % rank from lower bound K to upper bound 4K
[~,~,VR] = svd(X,'econ'); % take the right singular vector of X
clusterLabel = {};
approxBound = [];
Aff = {};
eigenValues = [];

% MDD
features = cell(1, size(X, 2));
for j=1:size(X, 2)
    t = reshape(X(:, j), 3, []);
    t(3, :) = [];
    v = diff(t,1,2);
    features{j} = v;
end

HH  = getHH(features,opt);
D = HHdist(HH, [], opt);
D = D / max(D(:));
Wj = exp(-D / 1);

for ii = 1:length(r)
    rnk = r(ii);
	V = VR(:,1:rnk);
	V = normr(V);
	Z = V*V';
	W = real(Z.^3.5); 
    W = W .* Wj;  
	Aff{ii} = W;
	[clusterLabel{ii},~,~] = ncutW(W,K);
	D = diag(1./sum(W));
	L = D*W;	
	eigenValues = eigs(L,K+1);
	approxBound(ii) = ComputeNcutValue(W,clusterLabel{ii})/(eigenValues(K)-eigenValues(K+1));
end

[minNcutValue, idx] = min(approxBound);
grp = clusterLabel{idx};

HHcenter = cell(1, K);
for i = 1:K
    HHcenter{i} = steinMean(cat(3,HH{grp(:,i)==1}));
end

end