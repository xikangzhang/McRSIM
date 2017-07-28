function y = JLDistance(seg1,seg2,nc)
% Jasen-Bregman Divergence Distance

D = size(seg1,2);
assert(size(seg2,2)==D);

% nr = 4;
if ~exist('nc','var')
    nc = 4;
end

% % take out mean
% seg1 = bsxfun(@minus, seg1, mean(seg1));
% seg2 = bsxfun(@minus, seg2, mean(seg2));

% H1 = blockHankel(seg1',[nr size(seg1,1)-nr/D+1]);
% H2 = blockHankel(seg2',[nr size(seg2,1)-nr/D+1]);
H1 = blockHankel(seg1',[size(seg1,1)-nc+1, nc]);
H2 = blockHankel(seg2',[size(seg2,1)-nc+1, nc]);

H1_p = H1 / (norm(H1*H1','fro')^0.5);
H2_p = H2 / (norm(H2*H2','fro')^0.5);

% y = 2 - norm(H1_p*H1_p' + H2_p*H2_p', 'fro');
% y = 2 - norm(H1_p'*H1_p + H2_p'*H2_p, 'fro');

HH1 = H1_p' * H1_p;
HH2 = H2_p' * H2_p;

HH1 = HH1 + 1e-6 * eye(nc);
HH2 = HH2 + 1e-6 * eye(nc);

y = log(det((HH1+HH2)/2)) - 0.5*log(det(HH1*HH2));
% y = log(det((HH1+HH2)/2));

end