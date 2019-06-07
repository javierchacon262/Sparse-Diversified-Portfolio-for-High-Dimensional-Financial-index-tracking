function v3 = sum_v3(v3, w, y3)
sm = sum(v3);
while sm > 1.000000000001 || sm < 0.99999999999
    v3 = min(1, max(0, w+y3));
    v3 = v3/norm(v3, 1);
    sm = sum(v3);
end
end