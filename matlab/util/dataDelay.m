function [x, y, camID] = dataDelay(x, y, delay, ratio)

D = size(x, 1);
N = size(x, 2);
F = size(x, 3);
temp = zeros(D, N, F-delay);
rng(0);
ind = (rand(1, N) > ratio);
camID = zeros(N, 1);
temp(:, ind, :) = x(:, ind, 1:F-delay); camID(ind) = 1;
temp(:, ~ind, :) = x(:, ~ind, delay+1:F); camID(~ind) = 2;
x = temp;
temp(:, ind, :) = y(:, ind, 1:F-delay); camID(ind) = 1;
temp(:, ~ind, :) = y(:, ~ind, delay+1:F); camID(~ind) = 2;
y = temp;

end