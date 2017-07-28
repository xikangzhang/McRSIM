function HH_hat = pcaClean(HH)

w = size(HH,2);
h = size(HH,1);
d = min([w,h]);
% assert(rank(HH)==d);    % check that HH is full rank
% [coef,score,~,~,~,mu] = pca(HH','Centered',false);
% HH_hat = bsxfun(@plus,coef(:,1:d-1)*score(:,1:d-1)',mu'); % reduce rank by one
[U,S,V] = svd(HH);
S(end,end) = 0;
HH_hat = U*diag(diag(S))*U';
HH_hat = (HH_hat+HH_hat')/2;

end