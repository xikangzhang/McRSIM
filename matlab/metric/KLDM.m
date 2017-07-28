function d = KLDM(X, Y)

I = eye(size(X));
d = trace(X \ Y + Y \ X - 2 * I) / 2;

end