function seg = removeShortTraj(seg, minLen)

l = cellfun(@(x)size(x, 2), seg);
seg(l < minLen) = [];

end