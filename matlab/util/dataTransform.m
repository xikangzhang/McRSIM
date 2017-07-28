function [x,y,camID] = dataTransform(x, y, K, theta, T, ratio)

if nargin < 4
    theta = 0;
end
if nargin < 5
    T = [0; 0];
end
if nargin < 6
    ratio = 0.5;
end

y3 = permute(y, [1 3 2]);
x3 = permute(x, [1 3 2]);
N = size(y3, 3);

rng(0);
ind = (rand(1, N) > ratio);
camID = ones(N, 1);
for k = 1:N
    if ind(k) == true
        y3(1:2, :, k) = rotation(y3(1:2, :, k), theta); % rotation
        y3(1:2, :, k) = bsxfun(@plus, y3(1:2, :, k), T); % translation
        x3(:, :, k) = K \ squeeze(y3(:, :, k));
        camID(k) = 2;
    end
end

y = permute(y3, [1 3 2]);
x = permute(x3, [1 3 2]);


end