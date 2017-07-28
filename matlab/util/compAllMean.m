% compuate mean of all combinations of a vector

function [mn, index] = compAllMean(y)

n = length(y);
mn = zeros(1, 2^n-1);
index = cell(1, 2^n-1);
count = 1;
for i = 1:n
    ind = nchoosek((1:n), i);
    for j = 1:size(ind, 1)
        mn(count) = mean(y(ind(j, :)));
        index{count} = ind(j, :);
        count = count + 1;
    end
end
mn(count:end) = [];
index(count:end) = [];

end