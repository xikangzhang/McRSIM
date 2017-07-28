function [x, y, camID] = dataRotation(x, y, K, theta, ratio)

T = [0; 0];
[x,y,camID] = dataTransform(x, y, K, theta, T, ratio);

end