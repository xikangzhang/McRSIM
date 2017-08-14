function [missrate, grp, bestRank, minNcutValue,W, index] = RSIM_MDD(X, s, UpperD, LowerD, verbose)
% Inputs:
% X: data matrix
% s: groundtruth labels
% UpperD: rank upper bound divided by number of motions
% LowerD: rank upper bound divided by number of motions
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
[~,~,VR] = svd(X,'econ'); % take the right singular vector of X
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
    v = diff(t, 1, 2);
    features{j} = v;
end

opt.metric = 'JBLD';
opt.sigma = 10^-4;
opt.H_structure = 'HtH';
opt.H_rows = 10;

HHt  = getHH(features, opt);
D = HHdist(HHt, [], opt);
D = D / max(D(:));
Wj = exp(-D / 1);

for ii = 1:length(r)
    rnk = r(ii);
	V = VR(:,1:rnk);
	
	V = normr(V); % normalize each row
	
	Z = V*V'; % new shape interaction matrix	
		
	W = real(Z.^3.5); % enhance block-diagonal structure; 
	                  %	On hopkins155, average err = 0.79% with powering value gamma = 3.8.
					  % You can also try other powering values in [3,4].

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