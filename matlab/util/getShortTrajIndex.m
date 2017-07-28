function indShort = getShortTrajIndex(m, minLen)

N = size(m, 1);
indShort = [];
for kk = 1:N
    if nnz(m(kk, :)) < minLen
        indShort = [ indShort, kk ];
        continue;
    end
    Temp = fLabel2sLabel(m(kk, :));
    idx = (Temp(:,1)==1);
    l = Temp(idx,3) - Temp(idx,2) + 1;
    if all(l < minLen) || isempty(l)
        indShort = [ indShort, kk ];
    end
end

end