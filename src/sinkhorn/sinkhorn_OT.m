function [Map,Dis,u,v,iter,err_hold] = sinkhorn_OT(a,b,dis,lambda,tolerance,maxiter,verbose)
% compute entropy regularized OT distance and lower bounds on the EMD.
% Inputs:
%  a: m*1, marginal measure
%  b: n*1, marginal measure
%  M: m*n, distance between site_a and site_b
%  tolerence: stop creteria
%  iter: max iteration
%  verbose: whether display the iteration
%
% Outputs:
%  Map: optimal transportation map
%  Dis: optimal transport distance and lower bounds on the EMD
%  u: m*1, left scaling
%  v: n*1, right scaling
%  iter: iteration times
%  err_hold: errors at iteration steps
%% processing inputs
if nargin < 3
    error('no enough inputs');
end
if nargin < 4 || isempty(lambda)
    lambda = 1/100;
end
if(nargin < 5 || isempty(tolerance))
    tolerance = 1.0e-5;
end
if (nargin < 6 || isempty(maxiter))
    maxiter = 5000;
end
if (nargin < 7 || isempty(verbose))
    verbose = 0;
end

%% preprocessing: check the dimensions of a,b,and dis
if size(a,1) ~= size(dis,1) || size(b,1) ~= size(dis,2)
    error('check the inputs');
end

%% preprocessing: remove zero values of a(i) in 1-vs-n case

index = (a>0);
a = a(index);
dis = dis(index,:);
K = exp(-1/lambda*dis);
%K = K + repmat((1.0e-200)*(sum(K,1) == 0),length(K(:,1)),1);
U = K.*dis;
ainvK = bsxfun(@rdivide,K,a); % K./[a]

%% iterator, initialization
cnt = 1;
u = ones(size(a,1),size(b,2))/length(index);% апоРа©
dis_hold = 1;

err_hold = zeros(maxiter,3);
%% iterator
while cnt <= maxiter
    
    u = 1./(ainvK*(b./(K'*u)));
    cnt = cnt+1;
    
    if (mod(cnt,20) == 0 || cnt == maxiter)
        v = b./(K'*u);
        u = 1./(ainvK*v);
        
        
        dis = sum(u.*(U*v));% = v*K*u*M = T.*M
        Criterion = abs(dis./dis_hold-1);
        
        if Criterion < tolerance || isnan(Criterion)
            break;
        end
        
        dis_hold = dis;
        err_hold(cnt,1) = cnt;
        err_hold(cnt,2) = Criterion;
        err_hold(cnt,3) = dis_hold;
        
        if verbose > 0
            disp(['    Sinkhorn_Iteration :',num2str(cnt),' Criterion: ',num2str(Criterion)]);
        end
        if any(isnan(Criterion)), % stop all computation if a computation of one of the pairs goes wrong.
            error('NaN values appeared during the fixed point iteration, may be caused by high value of lambda.');
        end
        
        cnt = cnt+1;
        
    end
end

%disp(['Last Sinkhorn_Iteration :',num2str(cnt-1),' Criterion: ',num2str(Criterion)]);

iter = cnt -1;
Dis.disOT = sum(u.*(U*v));

alpha = log(u);
beta = log(v);
beta(beta == -inf) = 0; % zero values of v (corresponding to zero values in b) generate inf numbers.
Dis.lowEMD = (a'* alpha + sum(b.*beta))*lambda;

Map = bsxfun(@times,v',(bsxfun(@times,u,K)));




















