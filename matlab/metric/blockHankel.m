function H = blockHankel(L,nrnc)
%
%   H = blockHankel(X);
%
%   H = blockHankel(X,nrnc);
%
%   X: Data DIMS x N 
%   nrnc: size of hankel matrix nr x nc
%
% treat every column as a observation
% form a multi output hankel e.g.
%       | c1 c2 c3 |
% H =   | c2 c3 c4 |
%       | c3 c4 c5 |
%
% original hankel_mo.m by Caglayan Dicle


[dim, N] = size(L);

if dim>N
    warning('DIMS>N. Make sure X is DIMSxN (row vector)!')
end

if nargin<2
    nr = ceil(N/(dim+1))*dim;
    nc = N - ceil(N/(dim+1))+1;
else
    nr = nrnc(1);
    nc = nrnc(2);
end

cidx = (0 : nc-1);
ridx = (1 : nr)';

H = ridx(:,ones(nc,1)) + dim*cidx(ones(nr,1),:);  % Hankel subscripts 
t = L(:);

temp.type = '()';
temp.subs = {H(:)};
H = reshape( subsref( t, temp ), size( H ) );