function [missrate, grp, bestRank,W, index] = McRSIM_MDD_Incomplete(X, mask, s, UpperD, LowerD, camID)
% Inputs:
% X: data matrix
% mask: indicator matrix of missing data
% s: groundtruth label
% UpperD: rank upper bound divided by number of motions
% LowerD: rank upper bound divided by number of motions
% camID: camera ID of each trajectory

if(nargin<4)
	LowerD = 1;
end
if(nargin<3)
	UpperD = 4;
end
K = max(s);
LowerR = LowerD*K;
UpperR = UpperD*K;
r = LowerR:UpperR;
X1 = X(camID==1, :);
X2 = X(camID==2, :);
M1 = mask(camID==1, :);
M2 = mask(camID==2, :);
VR1 = findVR(X1, M1, UpperR);
VR2 = findVR(X2, M2, UpperR);

clusterLabel = {};
approxBound = [];
Aff = {};
eigenValues = [];

% MDD
opt.metric = 'JBLD';
opt.sigma = 10^-4;
opt.H_structure = 'HtH';
opt.H_rows = 5;

HH  = getHH_missing(X', mask', opt);
D = HHdist(HH, [], opt);
D = D / (max(D(:))+1e-10);
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