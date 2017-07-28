function val = quickNcutValue(W, K)

[clusterLabel,~,~] = ncutW(W,K);
Lambda = diag(1./sum(W));
L = Lambda*W;
eigenValues = eigs(L,K+1);	% you can also easily modify the ncutW function and
% let it output the eignvalues to save the above three steps
val = ComputeNcutValue(W,clusterLabel)/(eigenValues(K)-eigenValues(K+1));

end