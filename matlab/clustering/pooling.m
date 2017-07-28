function [Y] = pooling(X_all, poolMaxSize)

if nargin < 2
    poolMaxSize = 10000;
end

rng('default');
nX = length(X_all);
ri = randperm(nX);
Y(1:poolMaxSize) = struct('dsca',[], 'H', [], 'HH',[], 'sigma',[],'loc',[]);
counter = 0;
for i = 1:nX
    X = X_all{ri(i)};
    nd = length(X);
    counterEnd = counter + nd;
    if counterEnd > poolMaxSize,
        break;
    end
    Y(counter+1:counterEnd) = X;
    counter = counterEnd;
end
Y(counter+1:end) = [];

end