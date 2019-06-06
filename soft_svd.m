function M = soft_svd(M, tau)
s1 = size(M, 1);
[U, S, V] = svd(double(M));
for i = 1:s1
    S(i, i) = abs(S(i, i) - tau);
end
M = U*S*V';
end