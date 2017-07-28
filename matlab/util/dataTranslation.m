function [x, y, camID] = dataTranslation(x, y, K, T, ratio)

theta = 0;
[x,y,camID] = dataTransform(x, y, K, theta, T, ratio);

end