function segment = fLabel2sLabel(label)
% change frame label format to segment format

dif = diff(label);
if nnz(dif) == 0
    segment = [label(1), 1, length(label)];
    return;
end

ind = find(dif);
index = [1; reshape(ind, [], 1) + 1];

interval = zeros(length(index), 2);
interval(1:end-1, 1) = index(1:end-1);
interval(1:end-1, 2) = index(2:end)-1;
interval(end, :) = [index(end), length(label)];

segment = [reshape(label(index), [], 1), interval];

end