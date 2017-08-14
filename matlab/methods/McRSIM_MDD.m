function [missrate, grp, bestRank, W, index] = McRSIM_MDD(X, s, UpperD, LowerD, camID, verbose)
% Inputs:
% X: data matrix
% s: groundtruth labels
% UpperD: rank upper bound divided by number of motions
% LowerD: rank upper bound divided by number of motions
% camID: camera ID of each trajectory
% verbose: indicator of whether to use homogeneous coordinates, if it exist
% use x-y, if not, use x-y-1

if(nargin<4)
	LowerD = 1;
end
if(nargin<3)
	UpperD = 4;
end
K = max(s);
r = LowerD*K:UpperD*K;
[~,~,VR1] = svd(X(:,camID==1),'econ'); % take the right singular vector of X
[~,~,VR2] = svd(X(:,camID==2),'econ'); % take the right singular vector of X
clusterLabel = {};
approxBound = [];
Aff = {};
eigenValues = [];

% MDD
features = cell(1, size(X, 2));
for j=1:size(X, 2)
    if ~exist('verbose','var')
        t = reshape(X(:, j), 3, []);
        t(3, :) = [];
    else
        t = reshape(X(:, j), 2, []);
    end
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

for ii = 1:length(r)
    rnk = r(ii);
	V1 = VR1(:,1:rnk); V1 = normr(V1);
    V2 = VR2(:,1:rnk); V2 = normr(V2);
    n1 = size(V1, 1);
    n2 = size(V2, 1);
    if n1 > n2
        M = V2' * V1(1:n2, :);
        [Um, Sm, Vm] = svd(M);
        R = Um * Vm';
        V2t = V2 * R;
        W = [V1; V2t] * [V1', V2t'];
    else
        M = V1' * V2(1:n1, :);
        [Um, Sm, Vm] = svd(M);
        R = Um * Vm';
        V1t = V1 * R;
        W = [V1t; V2] * [V1t', V2'];
    end
    W = real(W.^3.5);
    
    W = W .* Wj;
                      
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

ind = (s~=0);
grp = (grp(ind,:));
s = s(ind);

% missrate = ErrorRate(grp, s); % calculate the error rate
[missrate, index] = ErrorRate2(grp, s); % calculate the error rate

end