clc
clear
close all
dbstop if error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Market data lecture
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data   = load('MarketSP500.mat');
stocks = Data.DataMat;
gspc   = Data.gspc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Price based returns calculation
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
%Returns correlation matrix calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rCov =          cov(returns_M); %Covariance matrix
% [rCor, sigma] = corrcov(rCov);  %Correlation matrix
%figure
%subplot(1, 2, 1)
%imagesc(rCor)
%colorbar
%title('Returns Correlation Matrix:')
%sumvec = zeros(1, s2);
%for i = 1:s2
%    sumvec(i) = sum(rCor(:, i));
%end
%sumvec = (1/max(sumvec)).*sumvec;
%x = 1:s2;
%subplot(1, 2, 2)
%scatter(x, sumvec)
%title('Normalized Stocks Correlation Sum:')
%xlabel('stock number')
%ylabel('normalized correlation sum')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ADMM implementation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Variable initialization
w      = rand(s2, 1);
v1     = zeros(s2, 1);
v2     = zeros(s1, s2);
v3     = zeros(s2, 1);
y1     = zeros(s2, 1);
y2     = zeros(s1, s2);
y3     = zeros(s2, 1);
iter   = 10000;
ro     = 0.0000001;
mu     = 0.0000001;
lambda = 0.0000001;
alpha  = lambda/ro;
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
    v1       = soft(v1, alpha);
    v2       = soft_svd(v2, -ro);
    v3       = sum_v3(v3);
    y1       = y1 + (ro*(v1 - w));
    y2       = y2 + (ro*(v2 - Dfor(X, w)));
    y3       = y3 + (ro*(v3 - w));
    norm1vec = [norm1vec, norm(w, 1)];
    norm2vec = [norm2vec, norm((X*w) - rb)];
    normNvec = [normNvec, norm(svd(Dfor(X, w)),1)];
    obj      = norm2vec(i) + lambda*norm1vec(i) - normNvec(i);
    objfunc  = [objfunc, obj];
end
xiter = 1:iter;
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
xr2 = 1:340;
plot(xr2, X*w)
hold on
plot(xr2, rb)

