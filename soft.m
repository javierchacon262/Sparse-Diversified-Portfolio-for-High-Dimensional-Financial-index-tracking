function vec = soft(vec, alpha)
s = size(vec, 1);
for i = 1:s
    if vec(i) - alpha <= 0
        vec(i) = 0;
    else
        vec(i) = vec(i) - alpha;
    end
end
end