function vec = Dinv(X, v)
sx1 = size(X, 1);
sx2 = size(X, 2);
sv2 = size(v, 2);
mat = zeros(sx1, sx2);
onv = ones(sx1, 1);
for i = 1:sx1
    for j = 1:sx2
        if X(i, j) ~= 0
            mat(i, j) = 1/X(i, j);
        else
            mat(i, j) = 0;
        end
    end
end
vec = (mat .* v)'*onv;
end