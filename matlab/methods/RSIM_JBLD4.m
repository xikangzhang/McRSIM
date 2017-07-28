function [missrate, grp, bestRank, minNcutValue,W] = RSIM_JBLD4(X, s, UpperD, LowerD, camID)
% Input: X --  data matrix, s -- groundtruth label, UpperD -- largest rank
% Pan Ji, pan.ji@anu.edu.au
if(nargin<4)
	LowerD = 1;
end
if(nargin<3)
	UpperD = 4;
end
K = max(s);
r = LowerD*K:UpperD*K; % rank from lower bound K to upper bound 4K
% [~,~,VR] = svd(X,'econ'); % take the right singular vector of X
[~,~,VR11] = svd(X(:,camID==1),'econ'); % take the right singular vector of X
[~,~,VR22] = svd(X(:,camID==2),'econ'); % take the right singular vector of X
clusterLabel = {};
approxBound = [];
Aff = {};
eigenValues = [];

% JBLD
features = cell(1, size(X, 2));
for j=1:size(X, 2)
    t = reshape(X(:, j), 3, []);
    t(3, :) = [];
    v = diff(t,1,2);
    %         v = diff(v,1,2);
    features{j} = v;
end

opt.metric = 'JBLD';
opt.sigma = 10^-4;
opt.H_structure = 'HtH';
opt.H_rows = 10;

HHt  = getHH(features,opt);
D = HHdist(HHt, [], opt);
% Wj2 = HHdist_ms(HHt,HHt,opt);
% Wj2 = Wj2.^(-1/log(min(Wj2(:))));
% D = D / norm(D);
% D = D .* D;
D = D / max(D(:));

Wj = exp(-D / 1);
% NcutDiscrete = ncutW(W, nCluster);
% label = sortLabel_order(NcutDiscrete, 1:size(D,1));

for ii = 1:length(r)
    rnk = r(ii);
	V11 = VR11(:,1:rnk);	
    V11 = normr(V11); % normalize each row
    Z11 = V11*V11'; % new shape interaction matrix
    W11 = real(Z11.^3.5); % enhance block-diagonal structure;
    V22 = VR22(:,1:rnk);
    V22 = normr(V22); % normalize each row
    Z22 = V22*V22'; % new shape interaction matrix
    W22 = real(Z22.^3.5); % enhance block-diagonal structure;

%     W = W .* Wj;
    W12 = Wj(camID==1, camID==2);
    W21 = Wj(camID==2, camID==1);
    W = [W11, W12; W21, W22];
                      
	Aff{ii} = W;
	[clusterLabel{ii},~,~] = ncutW(W,K); % install your ncutW function 
	                                     % from https://www.cis.upenn.edu/~jshi/software/ 
		
	D = diag(1./sum(W));
	L = D*W;	
	eigenValues = eigs(L,K+1);	% you can also easily modify the ncutW function and 
	                            % let it output the eignvalues to save the above three steps
	approxBound(ii) = ComputeNcutValue(W,clusterLabel{ii})/(eigenValues(K)-eigenValues(K+1));
end

[minNcutValue, idx] = min(approxBound);
W = Aff{idx};
grp = clusterLabel{idx};
bestRank = r(idx);
missrate = ErrorRate(grp, s); % calculate the error rate

end