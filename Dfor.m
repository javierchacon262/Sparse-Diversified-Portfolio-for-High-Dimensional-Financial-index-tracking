function mat = Dfor(X, w)
sw = size(X, 1);
mat = X .* (ones(sw, 1)*w');
end