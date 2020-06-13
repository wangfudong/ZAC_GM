function K=con_K(x,y,beta)
% CON_K constructs the kernel K, 
%   where K(i, j) = k(x, y) = exp(-beta*||x-y||^2).

% Authors: Jiayi Ma (jyma2010@gmail.com)
% Date:    04/17/2012

[n, d]=size(x); [m, d]=size(y);

K=repmat(x,[1 1 m])-permute(repmat(y,[1 1 n]),[3 2 1]);
K=squeeze(sum(K.^2,2));
K=-beta * K;
K=exp(K);