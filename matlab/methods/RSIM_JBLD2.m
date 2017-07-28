function [missrate, grp, bestRank, minNcutValue,W] = RSIM_JBLD2(X, s, UpperD, LowerD)
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
opt.sigma = 10^-4;
opt.H_structure = 'HtH';
opt.H_rows = 10;

kNN_list = [40];
scale_list = [10];
rankThres_list = (0.5:0.1:0.9)';

% HHt  = getHH(features,opt);
% D = HHdist(HHt, [], opt);
% % Wj2 = HHdist_ms(HHt,HHt,opt);
% % Wj2 = Wj2.^(-1/log(min(Wj2(:))));
% D = D / max(D(:));
% Wj = exp(-D / 1);
% NcutDiscrete = ncutW(W, nCluster);
% label = sortLabel_order(NcutDiscrete, 1:size(D,1));


nr = length(r);
nk = length(kNN_list);
ns = length(scale_list);
nt = length(rankThres_list);
approxBound = zeros(nr, nk, ns);
clusterLabel = cell(nr, nk, ns);
Aff = cell(nr, nk, ns);
for ii = 1:length(r)
    for kk = 1:nk
        for ss = 1:ns
            for tt = 1:nt
            
%             kNN = min(kNN_list(kk), length(s));
%             scale_sig = scale_list(ss);
%             D1 = D - min(D(:));
%             D2 = D1;
%             for j=1:size(D1,1)
%                 [~,ind] = sort(D1(:,j));
%                 D2(ind(kNN+1:end),j) = Inf;
%                 D2(ind(1:kNN),j) = D1(ind(1:kNN),j) / max(D1(ind(1:kNN),j)) * 0.5;
%             end
%             % D = min(D2,D2');%(B+B')/2;
%             D1 = (D2+D2')/2;
            % D = max(D2,D2');
%             Wj = exp(-D1.^2/(2*scale_sig^2));
%             Wj = exp(-D1/scale_sig);
%             Wj = exp(-D/max(D(:)));
            
            rnk = r(ii);
            V = VR(:,1:rnk);
            
            V = normr(V); % normalize each row
            
            Z = V*V'; % new shape interaction matrix
            
            W = real(Z.^3.5); % enhance block-diagonal structure;
            %	On hopkins155, average err = 0.79% with powering value gamma = 3.8.
            % You can also try other powering values in [3,4].
            
            opt.thres = rankThres_list(tt);
            HH  = getHH(features, opt);
            D = HHdist(HH, [], opt);
            % Wj2 = HHdist_ms(HHt,HHt,opt);
            % Wj2 = Wj2.^(-1/log(min(Wj2(:))));
            D = D / max(D(:));
            Wj = exp(-D / 1);
            
            W = W .* Wj;
%             W = Wj;
            
            Aff{ii, kk, ss} = W;
            [clusterLabel{ii, kk, ss},~,~] = ncutW(W,K); % install your ncutW function
            % from https://www.cis.upenn.edu/~jshi/software/
            
            Lambda = diag(1./sum(W));
            L = Lambda*W;
            eigenValues = eigs(L,K+1);	% you can also easily modify the ncutW function and
%             opts.tol = 1e-3; eigenValues = eigs(L,K+1,'LR',opts);
%             ev = eig(L);
%             eigenValues = ev(end:-1:1);
            % let it output the eignvalues to save the above three steps
            approxBound(ii, kk, ss) = ComputeNcutValue(W,clusterLabel{ii})/(eigenValues(K)-eigenValues(K+1));
            
            end
        end
    end
end

[minNcutValue, idx] = min(approxBound(:));
W = Aff{idx};
grp = clusterLabel{idx};
% bestRank = r(idx);
bestRank = 1;
missrate = ErrorRate(grp, s); % calculate the error rate

end