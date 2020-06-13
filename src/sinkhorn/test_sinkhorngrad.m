% test for SinkhornGrad
n = 30;
m = 30;

dx = 1/n:1/n:1;
dy = 1/m:1/m:1;

dxy = abs(bsxfun(@minus,dx',dy)).^1;
dxy = dxy/max(dxy(:));

a = exp(-(dx'-0.3).^2/0.1);
a = a/sum(a);
b = exp(-(dy'-0.7).^2/0.05);
b = b/sum(b);

% a = 1/n*ones(n,1);
% b = 1/m*ones(m,1);
% 
% a = 1*(dx'-0.3).^2+1;
% a = a/sum(a);
% b = 2*(dy'-0.7).^2+1;
% b = b/sum(b);

a = rand(n,1);
a = a/sum(a);
b = rand(m,1);
b = b./sum(b,1);

figure,plot(dx,a,'r');
figure,plot(dy,b,'b')
figure,imagesc(dxy);
%%
lambda = 1/200;
tol = 1.0e-7;
maxiter = 5000;
verbose = 0;
tic;
[map,Dis,u,v,iter,err_hold] = sinkhorn_OT(a,b,dxy,lambda,tol,maxiter,verbose); % running with VERBOSE
 [Dis1] = Sinkhorn_1vsN(a,b,dxy,lambda,tol,maxiter);
toc;

figure,imagesc(map),title('map T');
figure,surf(map),title('map T');

[map_lp] = LP_OT(a,b,dxy,1);
figure,imagesc(map_lp);
figure,surf(map_lp),title('map lp');
%%
K = exp(-1/lambda*dxy);
tic
[d,grad_a,grad_b] = SinkhornGrad(a,b,K,dxy,500);
toc
figure,imagesc(grad_a)
figure,imagesc(grad_b)























