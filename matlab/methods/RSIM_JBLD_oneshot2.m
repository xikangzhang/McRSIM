function [missrate, grp, bestRank, minNcutValue,W, index] = RSIM_JBLD_oneshot2(X, s, UpperD, LowerD, camID)
% Input: X --  data matrix, s -- groundtruth label, UpperD -- largest rank
% Pan Ji, pan.ji@anu.edu.au
if(nargin<4)
	LowerD = 1;
end
if(nargin<3)
	UpperD = 4;
end
% K = max(s);
K = max(s);
r = LowerD*K:UpperD*K; % rank from lower bound K to upper bound 4K
[~,~,VR] = svd(X,'econ'); % take the right singular vector of X
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
% opt.metric = 'JBLD_denoise'; opt.epsilon = 0.01;
opt.sigma = 10^-4;
opt.H_structure = 'HtH';
opt.H_rows = 10;

HHt  = getHH(features,opt);
D = HHdist(HHt, [], opt);
D = D / max(D(:));
Wj = exp(-D / 1);
% NcutDiscrete = ncutW(W, nCluster);
% label = sortLabel_order(NcutDiscrete, 1:size(D,1));

for ii = 1:length(r)
    rnk = r(ii);
	V = VR(:,1:rnk);
	
	V = normr(V); % normalize each row
	
	Z = V*V'; % new shape interaction matrix	
		
	W = real(Z.^3.5); % enhance block-diagonal structure; 
	                  %	On hopkins155, average err = 0.79% with powering value gamma = 3.8.
					  % You can also try other powering values in [3,4].
           
	Aff{ii} = W;
	[clusterLabel{ii},~,~] = ncutW(W,2*K); % install your ncutW function 
	                                     % from https://www.cis.upenn.edu/~jshi/software/ 
		
	Lambda = diag(1./sum(W));
	L = Lambda*W;	
    opts.tol = 1e-6;
	eigenValues = eigs(L,2*K+1, 'LM', opts);	% you can also easily modify the ncutW function and 
	                            % let it output the eignvalues to save the above three steps
	approxBound(ii) = ComputeNcutValue(W,clusterLabel{ii})/(eigenValues(2*K)-eigenValues(2*K+1));
end

[minNcutValue, idx] = min(approxBound);
W = Aff{idx};
grp = clusterLabel{idx};
bestRank = r(idx);

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
    
%     % KNN on W
%     kNN = 50;
%     W(camID==1, camID==2) = Wj(camID==1, camID==2).^300;
%     W(camID==2, camID==1) = Wj(camID==2, camID==1).^300;
%     W1 = W;
%         for j = 1:size(W1, 1)
%         [~,ind] = sort(W1(j,:),'descend');
%         W1(j, ind(kNN+1:end)) = 0;
%         end
%     W = (W1 + W1') / 2;

% % scale Wj to Wr
% W11 = W(camID==1, camID==1);
% W22 = W(camID==2, camID==2);
% Wrm = (norm(W11, 'fro') + norm(W22, 'fro')) / (numel(W11) + numel(W22));
% W12 = Wj(camID==1, camID==2).^1000;
% Wjm = norm(W12, 'fro') / numel(W12);
% scale = Wrm / Wjm;
% W12 = W12 * scale;
% W(camID==1, camID==2) = W12;
% W(camID==2, camID==1) = W12';

[grp,~,~] = ncutW(W,K);

ind = (s~=0);
grp = (grp(ind,:));
s = s(ind);

% missrate = ErrorRate(grp, s); % calculate the error rate
[missrate, index] = ErrorRate2(grp, s); % calculate the error rate

end