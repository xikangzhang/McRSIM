% The clustering method of Normalized cuts 
% Input:
%    D: the distance matrix
%    k: the number of clustering
%    order: the order information
% Output:
%    label: the unique index for each cluster


function label = Ncuts(D, k, order)

W = exp(-D);     % the similarity matrix
NcutDiscrete = ncutW(W, k);

% label the data from low order to high order
index = cell(k, 1);
order_index = zeros(k, 1);

for i = 1:k
    index{i} = find(NcutDiscrete(:, i));
    order_index(i) = sum(order(index{i}))/numel(index{i});
end

[~, I] = sort(order_index);

for i = 1:k
    label(index{I(i)}) = i;
end

end