function [d,grad_a,grad_b] = SinkhornGrad(a,b,K,M,niter)
% to compute the gradient of entropy-regularized wasserstein diatance w.r.t
% a and b

%% operator handles
u_update = @(v,a) (repmat(a,1,length(b(1,:)))./(K*v));
v_update = @(u,b) (b./(K'*u));

DuDat = @(x,dvdat,a,v) ( bsxfun(@rdivide,x,K*v)...
    -dvdat(K'*(bsxfun(@times,x,(repmat(a,1,length(b(1,:)))./((K*v).^2))))) );

DvDat = @(x,dudat,b,u) ( -dudat(K*(bsxfun(@times,x,(b./((K'*u).^2))))) );

DuDbt = @(x,dvdbt,a,v) (-dvdbt(K'*(bsxfun(@times,x,repmat(a,1,length(b(1,:)))./((K*v).^2)))));

DvDbt = @(x,dudbt,b,u) ( bsxfun(@rdivide,x,K'*u)...
    -dudbt(K*(bsxfun(@times,x,b./((K'*u).^2)))));

%% settings

n = size(a,1);
m = size(b,1);

DVDAT = @(EPS) zeros(n,size(EPS,2));
DVDBT = @(EPS) zeros(m,size(EPS,2));

v = ones(m,size(b,2));

%% iterations
for j = 1:niter
    u = u_update(v,a);
    DUDAT = @(x) ( DuDat(x,DVDAT,a,v) );
    DUDBT = @(x) ( DuDbt(x,DVDBT,a,v) );
    
    v = v_update(u,b);
    DVDAT = @(x) ( DvDat(x,DUDAT,b,u) );
    DVDBT = @(x) ( DvDbt(x,DUDBT,b,u) );
end

%% output

U = K.*M;

d = sum(u.*(U*v));

grad_a = ( DUDAT(U*v) + DVDAT(U'*u) );
grad_b = ( DUDBT(U*v) + DVDBT(U'*u) );




% % operator handles
% u_update = @(v,a) (a./(K*v));
% v_update = @(u,b) (b./(K'*u));
% 
% DuDat = @(x,dvdat,a,v) ( bsxfun(@rdivide,x,K*v)...
%     -dvdat(K'*(bsxfun(@times,x,(a./((K*v).^2))))) );
% 
% DvDat = @(x,dudat,b,u) ( -dudat(K*(bsxfun(@times,x,(b./((K'*u).^2))))) );
% 
% DuDbt = @(x,dvdbt,a,v) (-dvdbt(K'*(bsxfun(@times,x,a./((K*v).^2)))));
% 
% DvDbt = @(x,dudbt,b,u) ( bsxfun(@rdivide,x,K'*u)...
%     -dudbt(K*(bsxfun(@times,x,b./((K'*u).^2)))));
% 
% % settings
% 
% n = size(a,1);
% m = size(b,1);
% 
% DVDAT = @(EPS) zeros(n,size(EPS,2));
% DVDBT = @(EPS) zeros(m,size(EPS,2));
% 
% v = ones(m,size(b,2));
% 
% % iterations
% for j = 1:niter
%     u = u_update(v,a);
%     DUDAT = @(x) ( DuDat(x,DVDAT,a,v) );
%     DUDBT = @(x) ( DuDbt(x,DVDBT,a,v) );
%     
%     v = v_update(u,b);
%     DVDAT = @(x) ( DvDat(x,DUDAT,b,u) );
%     DVDBT = @(x) ( DvDbt(x,DUDBT,b,u) );
% end
% 
% % output
% 
% U = K.*M;
% 
% d = sum(u.*(U*v));
% 
% grad_a = ( DUDAT(U*v) + DVDAT(U'*u) );
% grad_b = ( DUDBT(U*v) + DVDBT(U'*u) );





































