function [idt, V] = L2E(X, Y, thresh,verbose)

% Authors: Jiayi Ma (jyma2010@gmail.com)
% Date:    04/17/2012

beta = 0.8;
lambda = 0;
is_grad = 1;

n_ker = 15;
[N, D] = size(X); 
tmp_X = unique(X, 'rows'); idx = randperm(size(tmp_X,1)); idx = idx(1:min(n_ker,size(tmp_X,1))); ctrl_X=tmp_X(idx,:);
K=con_K(ctrl_X, ctrl_X, beta);
U = con_K(X, ctrl_X, beta);

x0 = zeros(n_ker*D, 1);
sigma2 = 0.05;

%%
options = optimset( 'display',verbose, 'MaxIter', 50);
if is_grad
    options = optimset(options, 'GradObj', 'on');
end
param = fminunc(@(x)costfun(x, X, Y, K, U, lambda, sigma2, is_grad), x0, options);
for ii = 1:3
    sigma2 = sigma2*0.5;
    param = fminunc(@(x)costfun(x, X, Y, K, U, lambda, sigma2, is_grad), param, options);
end

C = param(1:end);
C = reshape(C, [n_ker D]);
V=U*C;
Pb =  exp(-sum((Y-V).^2, 2) / (2*sigma2)) ;

idt = find(Pb > thresh);
