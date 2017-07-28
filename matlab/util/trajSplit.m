function seg = trajSplit(S, v)
% split contour according to indices
% Input:
% S: 2-by-n matrix, the contour coordinates
% v: 1-by-n vector, indication vector, 1 means keeping, 0 means
% getting rid of
% Output:
% seg: 1-by-k cell array, segments of contours

n = size(S, 2);
assert(length(v)==n);

% if all points are valid return S
if all(v), seg{1} = S; return; end

%  split trajectory
dv = diff(v);
indS = find(dv==1)+1;
indE = find(dv==-1);
if v(1), indS = [1, indS]; end
if v(end), indE = [indE, n]; end
k = length(indS);
assert(length(indE)==k);
seg = cell(1, k);
for i = 1:k
    seg{i} = S(:, indS(i):indE(i));
end

end