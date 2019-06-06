function v3 = sum_v3(v3, w, y2)
sm = sum(v3);
while sm > 1
    v3 = min(1, max(0, w+y2));
    v3 = v3/norm(v3, 1);
end
end