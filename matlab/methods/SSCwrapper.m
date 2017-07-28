function [missrate,C1, label, gt] = SSCwrapper(X, s, verbose)

if nargin < 3
    X(3:3:end, :) = [];
end
r = 0; affine = true; alpha = 800; outlier = false; rho = 0.7;
[missrate,C1, label, gt] = SSC(X,r,affine,alpha,outlier,rho,s);

end