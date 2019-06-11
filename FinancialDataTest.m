clc
clear
close all
dbstop if error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Market data lecture
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data    = load('MarketSP500.mat');
stocks  = Data.DataMat;
gspc    = Data.gspc;
symbols = Data.symbols;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Price-based returns calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1    = size(stocks, 1);
s2    = size(stocks, 2);
X     = zeros(s1, s2);
rb    = zeros(s1, 1);
dates = [];
for i = 1:s1
    rb(i) = (gspc(i).Close - gspc(i).Open)/gspc(i).Open;
    dates = [dates; gspc(i).Date];
    for j = 1:s2
        X(i, j) = (stocks(i, j).Close - stocks(i, j).Open)/stocks(i, j).Open;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Variable initialization random
v1     = rand(s2, 1);
v2     = rand(s1, s2);
v3     = rand(s2, 1);
y1     = rand(s2, 1);
y2     = rand(s1, s2);
y3     = rand(s2, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameter initialization
iter   = 100;
ro     = 0.001;
mu     = 0.1;
lambda = 0.1;
alpha  = 0.001;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ciclo de optimizacion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X2       = X'*X;
sx2      = size(X2, 1);
norm1vec = [];
norm2vec = [];
normNvec = [];
objfunc  = [];
for i = 1:iter
    disp(i)
    w        = ((2*X2) + (3*ro*eye(sx2, sx2)))\((2*X'*rb) + ro*(v1 + Dinv(X, v2) + v3 - y1 - Dinv(X, y2) - y3));
    v1       = soft(w+y1, alpha);
    v2       = soft_svd(Dfor(X, w)+y2, -ro);
    v3       = sum_v3(v3, w, y3);
    y1       = y1 + (ro*(v1 - w));
    y2       = y2 + (ro*(v2 - Dfor(X, w)));
    y3       = y3 + (ro*(v3 - w));
    w(w<0) = 0;
    w = w/norm(w, 1);
    norm1vec = [norm1vec, norm(w, 1)];
    norm2vec = [norm2vec, norm((X*w) - rb)];
    %normNvec = [normNvec, norm(svd(Dfor(X, w)),1)];
    normNvec = [normNvec, norm(svd(X*w),1)];
    obj      = norm2vec(i) + lambda*norm1vec(i) - normNvec(i);
    objfunc  = [objfunc, obj];
end
xiter       = 1:iter;
[ws, idxws] = sort(w, 'descend');
ws(20:end)  = 0;
wss         = zeros(sx2, 1);
Port_s      = [];
Port_f      = [];
for i = 1:sx2
    wss(idxws(i)) = ws(i);
end
wss = wss/norm(wss, 1);
for i = 1:sx2
    if wss(i) ~= 0
        Port_s = [Port_s, symbols(i)];
        Port_f = [Port_f, wss(i)];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Restore price from returns for every portfolio and comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pw    = [gspc(1).Close];
pwss  = [gspc(1).Close];
rxw   = X*w;
rxwss = X*wss;
for i = 1:s1-1
    pw   = [pw, (pw(i) + pw(i)*rxw(i))];
    pwss = [pwss, (pwss(i) + pwss(i)*rxwss(i))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First graph, norms behaviour
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name', 'Norms Behaviour', 'WindowState', 'maximized')
subplot(1, 4, 1)
scatter(xiter, norm2vec)
title('||Xw - rb||^2_2')
xlabel('Iterations')
ylabel('L2-Norm')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1, 4, 2)
scatter(xiter, norm1vec)
title('||w||_1')
xlabel('Iterations')
ylabel('L1-Norm')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1, 4, 3)
scatter(xiter, normNvec)
title('||D_x(w)||_*')
xlabel('Iterations')
ylabel('L*-Norm')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1, 4, 4)
scatter(xiter, objfunc)
title('Objective function evaluation')
xlabel('Iterations')
ylabel('OBJ-Func')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(1, 2, 1)
xr2 = 1:s1;
plot(xr2, X*w)
hold on
plot(xr2, X*wss)
plot(xr2, rb)
title('Returns Comparison')
legend('X*w','X*wss','rb')
hold off
subplot(1, 2, 2)
spc = zeros(s1);
for i = 1:s1
    spc(i) = gspc(i).Close;
end
plot(xr2, spc)
hold on
plot(xr2, pw)
plot(xr2, pwss)
legend('Benchmark', 'PW', 'PWSS')
hold off
xlabel('Time Periods: Days')
ylabel('Index Price: USD$')
ylim([2000, 3000])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('')
disp('L2-Norm between the selected 20-asset portfolio and the benchmark:')
disp(norm(X*wss - rb))
disp('')
disp('L2-Norm between the selected 20-asset portfolio and the +0- (220-230)-asset portfolio:')
disp(norm(X*wss - X*w))
disp('')
disp('L2-Norm between the +0- (220-230)-asset portfolio and the benchmark:')
disp(norm(X*w - rb))
disp('')
disp('20-asset portfolio found:')
disp(Port_s)
disp(Port_f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%