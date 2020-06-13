function [OT_dis] = Sinkhorn_1vsN(a,b,dis,lambda,tol,maxiter)

if nargin < 3
    error('not enough inputs');% check marginal constraints by default
end
if nargin < 4 || isempty(lambda),
    lambda = 1/100;
end
if nargin < 5 || isempty(tol),
    tol = 1.0e-5;
end
if nargin < 6 || isempty(maxiter),
    maxiter = 5000;
end

%%
I = (a > 0);
a = a(I);
dis = dis(I,:);
K = exp(-1/lambda*dis);
U = K.*dis;
ainvK=bsxfun(@rdivide,K,a); % precomputation of this matrix saves a d1 x N Schur product at each iteration.

%% Fixed point counter
compt = 1;
Dold = ones(1,size(b,2)); %initialization of vector of distances.
%% Initialization of Left scaling Factors, N column vectors.
u = ones(size(a,1),size(b,2))/size(a,1);

while compt <= maxiter
    
    u = 1./(ainvK*(b./(K'*u)));
    compt = compt+1;
    
    % check the stopping criterion every 20 fixed point iterations
    % or, if that's the case, before the final iteration to store the most
    % recent value for the matrix of right scaling factors v.
    if mod(compt,20) == 0 || compt == maxiter,
        v = b./(K'*u); % main iteration of Sinkhorn's algorithm
        u = 1./(ainvK*v);
        
        D = sum(u.*(U*v));
        Criterion = max(abs(D./Dold-1));
        if Criterion < tol || isnan(Criterion)
            break;
        end
        Dold = D;
        
        if any(isnan(Criterion)), % stop all computation if a computation of one of the pairs goes wrong.
            error('NaN values have appeared, lambda is too large,');
        end
        compt = compt+ 1;
    end
    
end
% disp(['Last_Sinkhorniter :',num2str(compt-1),' Criterion : ',num2str(Criterion)]);
compt = compt - 1;
OT_dis = sum(u.*(U*v));
end