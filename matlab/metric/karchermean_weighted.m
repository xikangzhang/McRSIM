function G = karchermean_weighted(Gi, w)

N = length(Gi);
if nargin<2, w = 1/N*ones(1,N); end
assert(length(w)==N);
G = Gi{1};
Gpre = G;
for iter = 1:1000
    [V,D] = eig(G);
    sqrtG = V*diag(sqrt(diag(D)))*V';
    sqrtGinv = V*diag(1./sqrt(diag(D)))*V';
    M = cellfun(@(x) logm(sqrtGinv*x*sqrtGinv),Gi,'UniformOutput',false);
    w = reshape(w, [1 1 N]);
    M = sum(bsxfun(@times,cat(3,M{:}),w),3);
    G = sqrtG*expm(M)*sqrtG;
    if norm(G-Gpre,'fro')<=1e-3
        break;
    end
    Gpre = G;
end